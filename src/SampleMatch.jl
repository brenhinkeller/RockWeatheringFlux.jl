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
        rocklat = read(macrofid["rocklat"]),
        rocklon = read(macrofid["rocklon"]),
        age = read(macrofid["age"]),
        type = read(macrofid["typecategory"]),
        rocktype = read(macrofid["rocktype"]),
        rockname = read(macrofid["rockname"]),
        rockdescrip = read(macrofid["rockdescrip"]),
    )
    close(macrofid)
    macro_cats = match_rocktype(macrostrat.type)

    
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


## --- Create average geochemistry lookup table for each rock name
    # Get rock names for each Macrostrat sample
    name_cats = match_rockname(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    rocknames = string.(keys(name_cats))

    # Get EarthChem samples for each rock name
    typelist = get_rock_class(false, 1)[1]
    nbulk = length(bulktext.Rock_Name)
    bulk_lookup = NamedTuple{keys(name_cats)}([falses(nbulk) for _ in eachindex(name_cats)])

    p = Progress(length(rocknames)+1, desc="Finding EarthChem samples for each rock name")
    next!(p)
    for i in eachindex(rocknames)
        bulk_lookup[i] .= find_earthchem(rocknames[i], bulktext.Rock_Name, bulktext.Type, 
            bulktext.Material
        )

        # If no matches, jump up a class
        if count(bulk_lookup[i]) == 0
            newsearch = class_up(typelist, rocknames[i])
            newsearch==:carb && (newsearch=="carbonate")    # No carbonatites
            bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
                bulktext.Type, bulktext.Material
            )

            # If still no matches, jump up a class again
            if count(bulk_lookup[i]) == 0
                newsearch = class_up(typelist, rocknames[i])
                bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
                    bulktext.Type, bulktext.Material
                )
            end
        end
        next!(p)
    end

    # Get average geochemistry for each rock name
    bulk_lookup = NamedTuple{keys(name_cats)}([major_elements(bulk, bulk_lookup[i]) 
        for i in eachindex(bulk_lookup)]
    )


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

    @timev @showprogress for i in eachindex(matches)
        # Progress bar
        # p = Progress(length(matches)+1, desc="Matching samples...")
        # next!(p)

        # Get the rock type and rock names of the sample
        type = get_type(macro_cats, i)
        name = get_type(name_cats, i, all_keys=true)
        (type==:cover || type==nothing) && continue

        # Get EarthChem data for that type
        bulksamples = bulk_cats[type]                        # EarthChem BitVector
        count(bulksamples) < 1 && continue

        EC = (
            bulklat = bulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = bulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = bulk.Age[bulksamples],                 # EarthChem age
            sampleidx = bulk_idxs[bulksamples],              # Indices of EarthChem samples
        )

        # Get all EarthChem samples for that rock type
        bulkgeochem = NamedTuple{Tuple(geochemkeys)}(
            [bulkzero[i][bulksamples] for i in geochemkeys]
        )

        # Average geochemistry by rock names
        # geochemdata = major_elements(bulk, bulksamples)
        geochemdata = bulk_lookup[name]

        # Find match
        matches[i] = likelihood(EC.bulkage, macrostrat.age[i], EC.bulklat, EC.bulklon, 
            macrostrat.rocklat[i], macrostrat.rocklon[i], bulkgeochem, geochemdata, 
            EC.sampleidx
        )

        # next!(p)
    end

    # Write data to a file
    writedlm("$matchedbulk_io", matches, "\t")


## --- End of File
