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
    using Static
    using Measurements

    # Local utilities
    include("utilities/Utilities.jl")

## --- Load Macrostrat data
    @info "Loading Macrostrat lithologic data ($macrostrat_io)"
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocktype = read(macrofid["rocktype"]),
        rockname = read(macrofid["rockname"]),
        rockdescrip = read(macrofid["rockdescrip"]),
        rocklat = read(macrofid["rocklat"]),
        rocklon = read(macrofid["rocklon"]),
        age = read(macrofid["age"]),
    )
    close(macrofid)
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip, major=false)


## --- Load Earthchem bulk geochemical data
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
    # Pre-define
    geochemkeys, = get_elements()               # Major elements
    bulk_idxs = collect(1:length(bulk.SiO2))    # Indices of bulk

    # Zero-NaN version of the major elements in bulk
    bulkzero = deepcopy(bulk)
    bulkzero = NamedTuple{Tuple(geochemkeys)}(
        [zeronan!(bulkzero[i]) for i in geochemkeys]
    )
    
    # Preallocate
    matches = zeros(Int64, length(macro_cats.sed))

    @timev for i in eachindex(matches)
        # Progress bar
        p = Progress(length(matches), desc="Matching samples...")

        # Get the type of the sample
        type = get_type(macro_cats, i)
        (type==:cover || type==nothing) && continue

        # Get EarthChem data for that type
        bulksamples = bulk_cats[type]                        # EarthChem BitVector
        count(bulksamples) < 1 && continue

        EC = (
            bulklat = bulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = bulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = bulk.Age[bulksamples],                 # EarthChem age
            sampleidx = bulk_idxs[bulksamples],              # Indices of EarthChem samples
            bulkname = bulktext.Rock_Name[bulksamples],      # Basalt, rhyolite, dacite, etc.
            bulktype = bulktext.Type[bulksamples],           # Volcanic, siliciclastic, etc.
            bulkmaterial = bulktext.Material[bulksamples],   # Ign, met, sed, xenolith, etc.
        )

        # Get all EarthChem samples for that rock type
        bulkgeochem = NamedTuple{Tuple(geochemkeys)}(
            [bulkzero[i][bulksamples] for i in geochemkeys]
        )

        # Average geochemistry for that type
        geochemdata = major_elements(bulk, bulksamples)

        # Find match
        matches[i] = likelihood(EC.bulkage, macrostrat.age[i], EC.bulklat, EC.bulklon, 
            macrostrat.rocklat[i], macrostrat.rocklon[i], bulkgeochem, geochemdata, 
            EC.sampleidx
        )

        next!(p)
    end


## --- Write data to a file
    writedlm("$matchedbulk_io", matches, "\t")


## --- End of File
