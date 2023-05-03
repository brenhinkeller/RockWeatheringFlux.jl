## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem
    using ProgressMeter
    using StatsBase

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
    reduced_bulk = Array{Array}(undef, length(first(geochem)) + 3, 1)
    geochemkeys = keys(first(geochem))
    for i in 1:length(geochemkeys)
        reduced_bulk[i] = bulk[geochemkeys[i]]
    end
    reduced_bulk[end-2] = bulk.Latitude
    reduced_bulk[end-1] = bulk.Longitude
    reduced_bulk[end] = bulk.AgeEst
    bulk = NamedTuple{(geochemkeys..., :Latitude, :Longitude, :AgeEst)}(reduced_bulk)    

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

    # Match data to rock types
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip, major=false)


## --- For each Macrostrat sample, estimate the log likelihood that each EarthChem sample
#      accurately represents that sample in space, time, and geochemistry
# TO DO: likelihood vs log likelihood? Is this already log likelihood?
# TO DO: NaN -> for missing age and location isn't good, but some equivilent?
# TO DO: send samples through in chunks isntead of by rock type
	
    @info "Matching samples"
    bulk_idxs = collect(1:length(bulk.SiO2))
    matches = Dict{Symbol, Matrix{Int64}}()
    #@time @showprogress "Matching lithographic / geochemical samples..." for type in eachindex(macro_cats)
    @time for type in eachindex(macro_cats)        
        # Cover is not in EarthChem data; skip it
        if type==:cover continue end
        
        # If going to crash: don't
        if type==:sed
            @warn "Skipping $type to avoid SigKill"
            continue
        end
        
        # Temporary output
        println("$type")
        
        # NOTE: doing int vars in a function is (slightly) worse than this
        # Intermediate Earthchem variables (unaffected by chunking)
        bulksamples = bulk_cats[type]                      # EarthChem BitVector
        bulklat = bulk.Latitude[bulksamples]               # EarthChem latitudes
		bulklon = bulk.Longitude[bulksamples]              # EarthChem longitudes
        bulkage = bulk.AgeEst[bulksamples]                 # EarthChem age
	    sampleidx = bulk_idxs[bulksamples]                 # Indices of EarthChem samples
	    
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
			chunks = vcat(collect(1:50:len), len)
		else
			continue
		end
        
		for i in 1:length(chunks[1:end-1])
			# Get chunks
			j = chunks[i]
        	k = chunks[i+1]-1

		    # Get intermediate variables for the rock type we're looking at
		    lat = macrostrat.rocklat[macro_cats[type]]         # Macrostrat latitude
		    lon = macrostrat.rocklon[macro_cats[type]]         # Macrostrat longitude
		    sampleage = macrostrat.age[macro_cats[type]]       # Macrostrat age

		    # Find most likely sample
		    # TO DO: pass argument as a tuple?
		    # TO DO: is there any way to do fewer than ~several M computations per sample?
		    # matched_sample = likelihood(lat, lon, bulklat, bulklon, sampleidx, sampleage, bulkage, geochemdata, bulkgeochem)
		    # setindex!(matches, matched_sample, type)
		    
		    # Find most likely Earthchem sample in chunks of at most 100 Macrostrat samples
        
        	
        	matched_sample = likelihood(lat[j:k], lon[j:k], bulklat, bulklon, sampleidx,
        		sampleage, bulkage, geochemdata, bulkgeochem)
        end
    end

#=
Test a subset of potential matches

s = matches[:ign][1]
i = findfirst(macro_cats.ign)

println("
    $(rocktype[i]) 
    $(rockdescrip[i]) 
    $(rockname[i]) 
    $(rockstratname[i]) 
    $(rockcomments[i]) 
    $(age[i])"
)
bulk_lith[s]
bulk.AgeEst[s]
=#

## --- EOF
