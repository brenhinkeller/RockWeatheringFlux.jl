## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem
    using ProgressMeter
    using StatsBase
    using DelimitedFiles

    # Local utilities
    include("Utilities.jl")


## --- Load Earthchem bulk geochemical data
    @info "Loading EarthChem bulk data"
    bulk_raw = matopen("data/bulk.mat")
    bulk_dict = read(bulk_raw, "bulk")
    close(bulk_raw)
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

    # Filter ages younger than 0 or greater than the age of the earth
    # TO DO: see if I should redo AgeEst... otherwise I'm not using Age
    # invalid_age = vcat(findall(>(4000), bulk.Age), findall(<(0), bulk.Age))
    # bulk.Age[invalid_age] .= NaN

    # Get rock types
    bulk_cats = match_earthchem(bulk.Type, major=false)

    # Calculate mean and standard deviation of major element oxides for each rock type
    # TO DO: quality control on bulk rock names: why are there granites with 50% silica?
    # TO DO: separate igneous rocks by silica content
    @time geochem = major_elements(bulk, bulk_cats)

    # Reduce bulk to only the data we need
    # I don't know that this actually makes a difference though? Since we pass variables by reference...
    # Although unreduced is ~1.6 GB so maybe it matters for overall memory space
    reduced_bulk = Array{Array{Float64}}(undef, length(first(geochem)) + 4, 1)
    geochemkeys = keys(first(geochem))
    for i in 1:length(geochemkeys)
        reduced_bulk[i] = bulk[geochemkeys[i]]
    end
    reduced_bulk[end-3] = bulk.P2O5
    reduced_bulk[end-2] = bulk.Latitude
    reduced_bulk[end-1] = bulk.Longitude
    reduced_bulk[end] = bulk.AgeEst     # TO DO: replace or recalculate AgeEst?
    bulk = NamedTuple{(geochemkeys..., :P2O5, :Latitude, :Longitude, :Age)}(reduced_bulk)    

## --- Load Earthchem metadata
    @info "Loading EarthChem sample metadata"
    bulktext_raw = matopen("data/bulktext.mat")
    bulktext_dict = read(bulktext_raw, "bulktext")
    close(bulktext_raw)
    bulktext = NamedTuple{Tuple(Symbol.(keys(bulktext_dict)))}(values(bulktext_dict))

    # Ignoring sparse arrays for now because they make me sad
    # Correct indices for 0-index offset
    composition_idx = Int.(bulktext.index["Composition"] .+ 1)
    reference_idx = Int.(bulktext.index["Reference"] .+ 1)
    rockname_idx = Int.(bulktext.index["Rock_Name"] .+ 1)
    source_idx = Int.(bulktext.index["Source"] .+ 1)
    type_idx = Int.(bulktext.index["Type"] .+ 1)
    material_idx = Int.(bulktext.index["Material"] .+ 1)

    # Parse numeric codes in index into arrays
    # TO DO: NamedTuple?
    bulk_composition = lowercase.(string.(bulktext.Composition[composition_idx]))
    bulk_reference = lowercase.(string.(bulktext.Reference[reference_idx]))
    bulk_rockname = lowercase.(string.(bulktext.Rock_Name[rockname_idx]))
    bulk_source = lowercase.(string.(bulktext.Source[source_idx]))
    bulk_type = lowercase.(string.(bulktext.Type[type_idx]))
    bulk_material = lowercase.(string.(bulktext.Material[material_idx]))

    # All relevant rock type identifiers
    bulk_lith = bulk_type .* " " .* bulk_material .* " " .* bulk_rockname


## --- Load Macrostrat data
    @info "Loading Macrostrat lithologic data"
    macrostrat = importdataset("data/toy_responses.tsv", '\t', importas=:Tuple)     # Reduced size file
    # macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)

    # Match data to rock types
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip, major=false)


## --- For each Macrostrat sample, estimate the log likelihood that each EarthChem sample
#      accurately represents that sample in space, time, and geochemistry
# TO DO: NaN -> for missing age and location isn't good, but some equivilent?
    # age - pick negative sample age

# views: 19.881672 seconds (1.82 M allocations: 35.402 GiB, 27.09% gc time, 2.11% compilation time)
# not:   24.634141 seconds (6.79 M allocations: 35.645 GiB, 18.63% gc time, 32.48% compilation time) - FIRST RUN

# first run no views: 30.999972 seconds (7.16 M allocations: 35.664 GiB, 27.01% gc time, 26.28% compilation time)
# first run w/ views: 33.637817 seconds (7.17 M allocations: 35.665 GiB, 25.97% gc time, 24.22% compilation time)


# NOTE: random sample selection means this is hitting kill again... but it doesn't allocate?
    # full data: 1518.174307 seconds (9.26 M allocations: 3.886 TiB, 24.87% gc time, 0.02% compilation time)
    # toy data: 
	
    @info "Matching samples"
    bulk_idxs = collect(1:length(bulk.SiO2))
    # TO DO: NamedTuple instead of a Dict?
    matches = Dict{Symbol, Array{Int64}}(
    	:sed => Array{Int64}(undef, count(macro_cats.sed), 1),
    	:ign => Array{Int64}(undef, count(macro_cats.ign), 1),
    	:met => Array{Int64}(undef, count(macro_cats.met), 1)
    )
    #@time @showprogress "Matching lithographic / geochemical samples..." for type in eachindex(macro_cats)
    @time for type in eachindex(matches)
        @info "Matching samples for $type"

        # NOTE: doing int vars in a function is (slightly) worse than this
        # Intermediate Earthchem variables (unaffected by chunking)
        bulksamples = bulk_cats[type]                      # EarthChem BitVector
        @views bulklat = bulk.Latitude[bulksamples]        # EarthChem latitudes
		@views bulklon = bulk.Longitude[bulksamples]       # EarthChem longitudes
        @views bulkage = bulk.Age[bulksamples]             # EarthChem age
	    @views sampleidx = bulk_idxs[bulksamples]          # Indices of EarthChem samples
	    
	    geochemdata = geochem[type]                        # Major element compositions
	    
	    # Earthchem samples for only major elements for this rock type
		    bulkgeochem = Array{Array}(undef, length(geochemdata), 1)
		    for i in 1:length(geochemkeys)
		        bulkgeochem[i] = bulk[geochemkeys[i]][bulksamples] 
		    end
		    bulkgeochem = NamedTuple{geochemkeys}(bulkgeochem)
		    
		# Get start and end coordinates for sample chunks
		len = count(macro_cats[type]) 
		if len != 0
			chunks = vcat(collect(1:25:len), len+1)
		else
			continue
		end
        
        # Find most likely Earthchem sample in chunks of at most 100 Macrostrat samples
		for i in 1:length(chunks[1:end-1])
			# Get chunks
			j = chunks[i]
        	k = chunks[i+1]-1

		    # Intermediate Macrostrat variables for this chunk
		    @views lat = macrostrat.rocklat[macro_cats[type]]         # Macrostrat latitude
		    @views lon = macrostrat.rocklon[macro_cats[type]]         # Macrostrat longitude
		    @views sampleage = macrostrat.age[macro_cats[type]]       # Macrostrat age

            # for n in j:k 
            #     matched_sample = likelihood(lat[n], lon[n], bulklat, bulklon, sampleidx[n],
        	# 	sampleage, bulkage, geochemdata, bulkgeochem)
            # end

		    # Find most likely sample
		    # TO DO: pass argument as a tuple?
		    # TO DO: is there any way to do fewer than ~several M computations per sample?
        	matched_sample = likelihood(lat[j:k], lon[j:k], bulklat, bulklon, sampleidx,
        		sampleage, bulkage, geochemdata, bulkgeochem)
        	matches[type][j:k] .= matched_sample
        	GC.gc()		# Needed to stop SigKill
        end
    end


## --- Separate data by rock type
    # Create one long array of all indices
    # Note that because not all rocks in Macrostrat were matched, not all will have data
    allmatches = zeros(Int64, length(macro_cats.ign))
    allmatches[macro_cats.sed] .= matches[:sed]
    allmatches[macro_cats.ign] .= matches[:ign]
    allmatches[macro_cats.met] .= matches[:met]

    # Write data to a file
    writedlm("output/matched_bulkidx.tsv", vcat("bulkidx", allmatches),"\t")

    # Separate into subtypes
    allmatches = (
        siliciclast = allmatches[macro_cats.siliciclast],
        shale = allmatches[macro_cats.shale],
        carb = allmatches[macro_cats.carb],
        chert = allmatches[macro_cats.chert],
        evaporite = allmatches[macro_cats.evaporite],
        coal = allmatches[macro_cats.coal],
        sed = allmatches[macro_cats.sed],

        volc = allmatches[macro_cats.volc],
        plut = allmatches[macro_cats.plut],
        ign = allmatches[macro_cats.ign],

        metased = allmatches[macro_cats.metased],
        metaign = allmatches[macro_cats.metaign],
        met = allmatches[macro_cats.met],
    )

## --- EOF
