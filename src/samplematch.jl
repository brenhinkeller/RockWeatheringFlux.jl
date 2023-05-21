## --- Match Macrostrat samples to the most likely EarthChem sample
    # The thing to do here is probably to have nothing happen in global scope
    # Packages
    using MAT
    using JLD
    using StatGeochem
    using ProgressMeter
    using StatsBase
    using DelimitedFiles
    using StaticArrays
    using LoopVectorization

    # Profiling
    using Profile
    using PProf
    using BenchmarkTools
    using Test

    # Local utilities
    include("Utilities.jl")

## --- Load Macrostrat data
    @info "Loading Macrostrat lithologic data"
    # macrostrat = importdataset("data/toy_responses.tsv", '\t', importas=:Tuple)     # Reduced size file
    macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)

    # Match data to rock types
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip, major=false)

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
    geochemkeys = (:SiO2, :Al2O3, :Fe2O3T, :TiO2, :MgO, :CaO, :Na2O, :K2O, :MnO)

    # Reduce bulk to only the data we need
    reduced_bulk = Array{Array{Float64}}(undef, length(geochemkeys) + 4, 1)
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
    rockname_idx = Int.(bulktext.index["Rock_Name"] .+ 1)
    type_idx = Int.(bulktext.index["Type"] .+ 1)
    material_idx = Int.(bulktext.index["Material"] .+ 1)

    # reference_idx = Int.(bulktext.index["Reference"] .+ 1)        # Not needed
    # source_idx = Int.(bulktext.index["Source"] .+ 1)              # Not needed

    # Parse numeric codes in index into arrays
    # TO DO: NamedTuple?
    bulk_composition = lowercase.(string.(bulktext.Composition[composition_idx]))
    bulk_rockname = lowercase.(string.(bulktext.Rock_Name[rockname_idx]))
    bulk_type = lowercase.(string.(bulktext.Type[type_idx]))
    bulk_material = lowercase.(string.(bulktext.Material[material_idx]))

    # bulk_reference = lowercase.(string.(bulktext.Reference[reference_idx]))
    # bulk_source = lowercase.(string.(bulktext.Source[source_idx]))

	
## --- Find matching Earthchem sample for each Macrostrat sample
    @info "Matching EarthChem and Macrostrat samples"
    bulk_idxs = collect(1:length(bulk.SiO2))
    # TO DO: NamedTuple instead of a Dict?
    matches = Dict{Symbol, Array{Int64}}(
    	:sed => Array{Int64}(undef, count(macro_cats.sed), 1),
    	:ign => Array{Int64}(undef, count(macro_cats.ign), 1),
    	:met => Array{Int64}(undef, count(macro_cats.met), 1)
    )
    
    @timev for type in eachindex(matches)
        # Intermediate Earthchem variables
        bulksamples = bulk_cats[type]                   # EarthChem BitVector
        bulklat = bulk.Latitude[bulksamples]            # EarthChem latitudes
		bulklon = bulk.Longitude[bulksamples]           # EarthChem longitudes
        bulkage = bulk.Age[bulksamples]                 # EarthChem age
	    sampleidx = bulk_idxs[bulksamples]              # Indices of EarthChem samples
        bulkname = bulk_rockname[bulksamples]           # Rock names
        bulktype = bulk_type[bulksamples]               # Types--volcanic, plutonic, siliciclastic, etc.
        bulkmaterial = bulk_material[bulksamples]       # Ign, met, sed, xenolith, etc.
	    
	    # Earthchem samples for only major elements for this rock type
        bulkgeochem = Array{Array}(undef, length(geochemkeys), 1)
        for i in 1:length(geochemkeys)
            bulkgeochem[i] = zeronan!(bulk[geochemkeys[i]][bulksamples])
        end
        bulkgeochem = NamedTuple{geochemkeys}(bulkgeochem)

        # Macrostrat samples
        lat = macrostrat.rocklat[macro_cats[type]]                 # Macrostrat latitude
        lon = macrostrat.rocklon[macro_cats[type]]                 # Macrostrat longitude
        sampleage = macrostrat.age[macro_cats[type]]               # Macrostrat age
        rocktype = macrostrat.rocktype[macro_cats[type]]           # Sample rock type
        rockname = macrostrat.rockname[macro_cats[type]]           # Sample name
        rockdescrip = macrostrat.rockdescrip[macro_cats[type]]     # Sample description

        # Progress bar
        p = Progress(length(lat), desc="Matching $type samples...")

        @inbounds for i in eachindex(lat)
            geochemfilter = find_earthchem(rocktype[i], rockname[i], rockdescrip[i], bulkname, bulktype, bulkmaterial)
            geochemdata = major_elements(bulk, bulksamples, geochemfilter)
            matches[type][i] = likelihood(bulkage, sampleage[i], bulklat, bulklon, lat[i], lon[i], bulkgeochem, geochemdata)
            
            # Manual garbage collection or else the code will run out of memory
            if i % 7 == 0
                GC.gc()
            end
            next!(p)
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
