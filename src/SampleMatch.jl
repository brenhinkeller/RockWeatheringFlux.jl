## --- Match Macrostrat samples to the most likely EarthChem sample
    # The thing to do here is probably to have nothing happen in global scope
    # Packages
    using HDF5
    using StatGeochem
    using ProgressMeter
    using StatsBase
    using DelimitedFiles
    using StaticArrays
    using LoopVectorization

    # Local utilities
    include("utilities/Utilities.jl")

## --- Load Macrostrat data
    @info "Loading Macrostrat lithologic data"
    macrostrat = importdataset("$macrostrat_io", '\t', importas=:Tuple)
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip, major=false)


## --- Load Earthchem bulk geochemical data
    @info "Loading EarthChem data"
    bulkfid = h5open("output/bulk.h5", "r")

    # Bulk
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    bulk_cats = match_earthchem(bulk.Type, major=false)

    # Bulk rock name, type, and material
    path = bulkfid["bulktext"]["sampledata"]
    header = read(path["header"])
    index = read(path["index"])

    target = ["Rock_Name", "Type", "Material"]
    targetind = [findall(==(i), header)[1] for i in target]

    bulktext = NamedTuple{Tuple(Symbol.(target))}(
        [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] for i in eachindex(target)]
    )    

    close(bulkfid)


## --- Find matching Earthchem sample for each Macrostrat sample
    geochemkeys, = get_elements()               # Major elements
    bulk_idxs = collect(1:length(bulk.SiO2))
    matches = (
        sed = Array{Int64}(undef, count(macro_cats.sed), 1),
    	ign = Array{Int64}(undef, count(macro_cats.ign), 1),
    	met = Array{Int64}(undef, count(macro_cats.met), 1)
    )

    @timev for type in eachindex(matches)
        bulksamples = bulk_cats[type]                        # EarthChem BitVector
        EC = (
            bulklat = bulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = bulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = bulk.Age[bulksamples],                 # EarthChem age
            sampleidx = bulk_idxs[bulksamples],              # Indices of EarthChem samples
            bulkname = bulktext.Rock_Name[bulksamples],      # Basalt, rhyolite, dacite, etc.
            bulktype = bulktext.Type[bulksamples],           # Volcanic, siliciclastic, etc.
            bulkmaterial = bulktext.Material[bulksamples],   # Ign, met, sed, xenolith, etc.
        )
	    
	    # Earthchem samples for only major elements for this rock type
        bulkgeochem = NamedTuple{Tuple(geochemkeys)}(
            [zeronan!(bulk[i][bulksamples]) for i in geochemkeys]
        )

        # Macrostrat samples
        macrosamples = macro_cats[type]                              # Macrostrat BitVector
        MS = (
            lat = macrostrat.rocklat[macrosamples],                  # Macrostrat latitude
            lon = macrostrat.rocklon[macrosamples],                  # Macrostrat longitude
            sampleage = macrostrat.age[macrosamples],                # Macrostrat age
            rocktype = macrostrat.rocktype[macrosamples],            # Sample rock type
            rockname = macrostrat.rockname[macrosamples],            # Sample name
            rockdescrip = macrostrat.rockdescrip[macrosamples],      # Sample description
        )
        
        # Progress bar
        p = Progress(length(MS.lat), desc="Matching $type samples...")

        @inbounds for i in eachindex(MS.lat)
            geochemfilter = find_earthchem(MS.rocktype[i], MS.rockname[i], MS.rockdescrip[i], 
                EC.bulkname, EC.bulktype, EC.bulkmaterial
            )
            geochemdata = major_elements(bulk, bulksamples, geochemfilter)
            matches[type][i] = likelihood(EC.bulkage, MS.sampleage[i], EC.bulklat, EC.bulklon, 
                MS.lat[i], MS.lon[i], bulkgeochem, geochemdata
            )
            next!(p)
        end
    end


## --- Separate data by rock type
    # Create one long array of all indices
    allmatches = zeros(Int64, length(macro_cats.ign))
    allmatches[macro_cats.sed] .= matches[:sed]
    allmatches[macro_cats.ign] .= matches[:ign]
    allmatches[macro_cats.met] .= matches[:met]

    # Write data to a file
    writedlm("$matchedbulk_io", allmatches,"\t")


## --- End of File
