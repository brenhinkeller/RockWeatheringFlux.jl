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
        type = read(macrofid["type"]),
        rocktype = read(macrofid["rocktype"]),
        rockname = read(macrofid["rockname"]),
        rockdescrip = read(macrofid["rockdescrip"]),
    )
    close(macrofid)
    # macro_cats = match_rocktype(macrostrat.type)
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
        macrostrat.rockdescrip, unmultimatch=false, inclusive=false, source=:macrostrat
    )
    
## --- Load Earthchem bulk geochemical data
    bulkfid = h5open("output/bulk.h5", "r")

    # Bulk
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    # bulktype = read(bulkfid["bulk"]["type"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

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
    bulk_cats = match_rocktype(bulktext.Rock_Name, bulktext.Type, bulktext.Material; 
        unmultimatch=false, inclusive=false, source=:earthchem
    )


## --- Create average geochemistry lookup table for each rock name
    # Get rock names for each Macrostrat sample
    name_cats = match_rockname(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    rocknames = string.(keys(name_cats))

    # Get EarthChem samples for each rock name
    typelist = get_rock_class(false, true)      # Get subtypes, major types inclusive
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
                newsearch = class_up(typelist, string(newsearch))
                bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
                    bulktext.Type, bulktext.Material
                )
            end
        end
        next!(p)
    end

    # Get average geochemistry for each rock name
    geochem_lookup = NamedTuple{keys(name_cats)}([major_elements(bulk, bulk_lookup[i]) 
        for i in eachindex(bulk_lookup)]
    )

## --- Calculate inverse spatial weight for each rock name
    # This current method means samples without a latitude or longitude will never get picked.
    # In theory, I could undo that by assigning a weight to samples without spatial data;
    # probably the mean of the weights (although this could cause some weird statistical
    # behavior depending on the distribution...). In any case, it may not end up mattering
    # if there's enough samples to choose from.

    # # Preallocate
    # spatial_lookup = NamedTuple{keys(name_cats)}([fill(NaN, nbulk) for _ in eachindex(name_cats)])

    # @info "Calculating inverse spatial weights for each rock name"
    # for n in eachindex(keys(spatial_lookup))
    #     println("$n ($n/$(length(rocknames))) \n")
    #     spatial_lookup[n][bulk_lookup[n]] .= invweight_location(bulk.Latitude[bulk_lookup[n]], 
    #         bulk.Latitude[bulk_lookup[n]]
    #     )
    # end

    # # Save to file
    # A = Array{Float64}(undef, nbulk, length(spatial_lookup))
    # for i in eachindex(keys(spatial_lookup))
    #     A[:,i] = spatial_lookup[i]
    # end

    # fid = h5open("output/invspatial.h5", "w")
    #     fid["header"] = collect(rocknames)
    #     fid["k"] = A
    # close(fid)


## --- Alternatively, load spatial weights from a file
    # fid = h5open("output/invspatial.h5", "r")
    #     header = read(fid["header"])
    #     k = read(fid["k"])
    # close(fid)
    # # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    # p = 1.0 ./ k
    # zeronan!(p)
    # spatial_lookup = NamedTuple{Tuple(Symbol.(header))}(p[:,i] for i in eachindex(header))
    # header = collect(rocknames)
    # spatial_lookup = NamedTuple{Tuple(Symbol.(header))}(ones(nbulk) for _ in eachindex(header))


## --- Find matching Earthchem sample for each Macrostrat sample
    # Preallocate
    matches = zeros(Int64, length(macro_cats.sed))
    geochemkeys, = get_elements()               # Major elements
    bulk_idxs = collect(1:length(bulk.SiO2))    # Indices of bulk

    # Zero-NaN version of the major elements in bulk
    bulkzero = deepcopy(bulk)
    bulkzero = NamedTuple{Tuple(geochemkeys)}(
        [zeronan!(bulkzero[i]) for i in geochemkeys]
    )

    p = Progress(length(matches)+1, desc="Matching samples...")
    next!(p)
    @timev for i in eachindex(matches)
        # Get the rock type and randomly select one sample rock name
        type = get_type(macro_cats, i, all_keys=true)
        if type==(:cover,) || type===nothing
            next!(p)
            continue
        end

        # Get EarthChem data for that type
        bulksamples = falses(length(bulk_cats[1]))           # EarthChem BitVector
        for t in type
            bulksamples .|= bulk_cats[t]
        end

        if count(bulksamples) == 0
            next!(p)
            continue
        end

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

        # Randomly select one EarthChem sample from a randomly selected sample name, 
        # proportional to it's spatial weight. Assign an error sampled from a normal 
        # distribution with mean and standard deviation equal to the mean and standard 
        # deviation of that rock name.
        #
        # This will represent the assumed geochemistry of the Macrostrat sample. The
        # assumption of this method is that there are enough samples of each name that 
        # outliers get ironed out.
        name = rand(get_type(name_cats, i, all_keys=true))
        t = @. bulk_lookup[name] & bulksamples
        @assert count(t) > 0 "$i"

        randsample = rand(bulk_idxs[t])

        geochemdata = geochem_lookup[name]
        errs = NamedTuple{Tuple(geochemkeys)}([abs(randn()*geochemdata[i].e) for i in geochemkeys])

        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((bulkzero[i][randsample]), errs[i])) for i in geochemkeys]
        )

        # Find match
        matches[i] = likelihood(EC.bulkage, macrostrat.age[i], EC.bulklat, EC.bulklon, 
            macrostrat.rocklat[i], macrostrat.rocklon[i], bulkgeochem, geochemdata, 
            EC.sampleidx
        )

        next!(p)
    end

    # Write data to a file
    writedlm("$matchedbulk_io", matches, '\t')


## --- End of File