## --- Match Macrostrat samples to the most likely EarthChem sample
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
    using Dates

    # Local utilities
    include("utilities/Utilities.jl")

    # Start timer
    start = now()
    @info "Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))"


## --- Load matches
    @info "Loading matched Tuples $(Dates.format(now(), "HH:MM"))"
    fid = h5open("output/matches.h5", "r")

    # Macrostrat rock types
    header = read(fid["vars"]["macro_cats_head"])
    data = read(fid["vars"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Macrostrat rock names
    header = read(fid["vars"]["name_cats_head"])
    data = read(fid["vars"]["name_cats"])
    data = @. data > 0
    name_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # EarthChem rock types
    header = read(fid["vars"]["bulk_cats_head"])
    data = read(fid["vars"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # EarthChem rock names
    rocknames = read(fid["vars"]["bulk_lookup_head"])
    data = read(fid["vars"]["bulk_lookup"])
    data = @. data > 0
    bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([data[:,i] for i in eachindex(rocknames)])

    close(fid)


## --- Load Macrostrat data
    @info "Loading Macrostrat data ($macrostrat_io) $(Dates.format(now(), "HH:MM"))"
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["rocklat"]),
        rocklon = read(fid["rocklon"]),
        age = read(fid["age"]),
        type = read(fid["type"]),
        rocktype = read(fid["rocktype"]),
        rockname = read(fid["rockname"]),
        rockdescrip = read(fid["rockdescrip"]),
    )
    close(fid)

    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
    #     macrostrat.rockdescrip, unmultimatch=false, inclusive=false, source=:macrostrat
    # )
    
## --- Load Earthchem bulk geochemical data
    @info "Loading EarthChem data $(Dates.format(now(), "HH:MM"))"
    fid = h5open("output/bulk.h5", "r")

    # Bulk
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Bulk rock name, type, and material
    path = fid["bulktext"]["sampledata"]
    header = read(path["header"])
    index = read(path["index"])

    target = ["Rock_Name", "Type", "Material"]
    targetind = [findall(==(i), header)[1] for i in target]

    bulktext = NamedTuple{Tuple(Symbol.(target))}(
        [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] 
            for i in eachindex(target)
        ]
    )
    close(fid)

    # bulk_cats = match_rocktype(bulktext.Rock_Name, bulktext.Type, bulktext.Material; 
    #     unmultimatch=false, inclusive=false, source=:earthchem
    # )


## --- Create average geochemistry lookup table for each rock name
    @info "Creating EarthChem lookup tables $(Dates.format(now(), "HH:MM"))"

    # # Get rock names for each Macrostrat sample
    # name_cats = match_rockname(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    # rocknames = string.(keys(name_cats))

    # # Get EarthChem samples for each rock name
    # typelist = get_rock_class(false, true)      # Subtypes, major types include minors
    # nbulk = length(bulktext.Rock_Name)
    # bulk_lookup = NamedTuple{keys(name_cats)}([falses(nbulk) for _ in eachindex(name_cats)])

    # p = Progress(length(rocknames)+1, desc="Finding EarthChem samples for each rock name")
    # next!(p)
    # for i in eachindex(rocknames)
    #     bulk_lookup[i] .= find_earthchem(rocknames[i], bulktext.Rock_Name, bulktext.Type, 
    #         bulktext.Material
    #     )

    #     # If no matches, jump up a class
    #     if count(bulk_lookup[i]) == 0
    #         newsearch = class_up(typelist, rocknames[i])
    #         newsearch==:carb && (newsearch=="carbonate")    # No carbonatites!
    #         bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
    #             bulktext.Type, bulktext.Material
    #         )

    #         # If still no matches, jump up a class again
    #         if count(bulk_lookup[i]) == 0
    #             newsearch = class_up(typelist, string(newsearch))
    #             bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
    #                 bulktext.Type, bulktext.Material
    #             )
    #         end
    #     end
    #     next!(p)
    # end

    # Get average geochemistry for each rock name
    geochem_lookup = NamedTuple{keys(name_cats)}([major_elements(bulk, bulk_lookup[i]) 
        for i in eachindex(bulk_lookup)]
    )


## --- Remove all multimatches from major types
    minorsed, minorign, minormet = get_minor_types()
    
    for type in minorsed
        macro_cats.sed .&= .!macro_cats[type]
        bulk_cats.sed .&= .!bulk_cats[type]
    end
    for type in minorign
        macro_cats.ign .&= .!macro_cats[type]
        bulk_cats.ign .&= .!bulk_cats[type]
    end
    for type in minormet
        macro_cats.met .&= .!macro_cats[type]
        bulk_cats.met .&= .!bulk_cats[type]
    end


## --- Get weights for weighted-random selection of rock types and names
    typelist = get_rock_class()                         # Major types exclude minor types
    minortypes = (minorsed..., minorign..., minormet...)

    # Minor rock types
    p_type = (
        sed = float.([count(macro_cats[i]) for i in minorsed]),
        ign = float.([count(macro_cats[i]) for i in minorign]),
        met = float.([count(macro_cats[i]) for i in minormet])
    )
    [p_type[i] ./= sum(p_type[i]) for i in keys(p_type)]

    # Descriptive rock names
    p_name = NamedTuple{minortypes}(
        [[float.(count(name_cats[Symbol(typelist[i][j])])) for j in eachindex(typelist[i])] 
            for i in minortypes
    ])
    [p_name[i] ./= sum(p_name[i]) for i in keys(p_name)]

    # rocktypes = collect(keys(macro_cats))[1:end-1]      # No cover!
    # p_name = NamedTuple{Tuple(rocktypes)}(
    #     [[float.(count(name_cats[Symbol(typelist[i][j])])) for j in eachindex(typelist[i])] 
    #         for i in rocktypes
    # ])

    # ... That are actually descriptive. Zero nondescriptive names
    # Technically now the goal would be to remove the requirement that the algo. picks
    # major types, since there are some names in major types that are actually useful....
    # nondesc = nondescriptive()
    # for k in keys(nondesc)
    #     for i in eachindex(typelist[k])
    #         if typelist[k][i] in nondesc[k]
    #             p_name[k][i] = 0.0
    #         end
    #     end
    # end
    # [p_name[i] ./= sum(p_name[i]) for i in keys(p_name)]


## --- Load spatial weights
    # fid = h5open("output/invspatial.h5", "r")
    #     header = read(fid["header"])
    #     k = read(fid["k"])
    # close(fid)
    # p = 1.0 ./ k
    # zeronan!(p)
    # spatial_lookup = NamedTuple{Tuple(Symbol.(header))}(p[:,i] for i in eachindex(header))

    # rocktypes = keys(macro_cats)
    # spatial_lookup = NamedTuple{rocktypes}([fill(NaN, nbulk) for _ in eachindex(rocktypes)])
    # for k in keys(spatial_lookup)
    #     println("$k\n")

    #     spatial_lookup[k][bulk_cats[k]] .= invweight_location(bulk.Latitude[bulk_cats[k]], 
    #         bulk.Latitude[bulk_cats[k]]
    #     )
    # end


## --- Find matching Earthchem sample for each Macrostrat sample
    # As part of this process, we'll need to assume the geochemistry of the Macrostrat 
    # sample.
    # 
    # To do this and preserve any multi-modal distributions in the data, we'll randomly 
    # pick one rock name matched with the sample, and randomly select one EarthChem sample
    # that was also matched with that rock name.
    # 
    # The error for each major element will be randomly sampled from a normal distribution 
    # with a mean and standard deviation equal to the mean and standard deviation for that
    # major element within the selected rock name.
    # 
    # This method assumes there are enough samples for outliers to get ironed out.

    # Preallocate
    matches = zeros(Int64, length(macro_cats.sed))

    # Definitions
    geochemkeys = get_elements()[1]                             # Major elements
    bulk_idxs = collect(1:length(bulk.SiO2))                    # Indices of bulk samples
    minortypes = (sed = minorsed, ign=minorign, met=minormet)   # Tuple minor types

    # Zero-NaN version of the major elements in bulk
    bulkzero = deepcopy(bulk)
    bulkzero = NamedTuple{Tuple(geochemkeys)}(
        [zeronan!(bulkzero[i]) for i in geochemkeys]
    )

    @info "Starting sample matching $(Dates.format(now(), "HH:MM"))"
    p = Progress(length(matches), desc="Matching samples...")
    @timev for i in eachindex(matches)
        # Get the rock type. Remove cover as a matched type, if it's there
        type = get_type(macro_cats, i, all_keys=true)
        if type===nothing
            next!(p)
            continue
        else
            t = trues(length(type))
            for i in eachindex(type)
                t[i] = type[i] != :cover
            end

            if count(t)==0
                next!(p)
                continue
            else
                type = type[t]
            end
        end
        
        # Pick a random sample to act as the geochemistry for that sample
        samplenames = get_type(name_cats, i, all_keys=true)
        randname, type = get_descriptive_name(samplenames, p_name, type, p_type, typelist, 
            minortypes
        )
        # randsample = bulk_idxs[weighted_rand(spatial_lookup[rand(type)][bulk_lookup[randname]])]
        randsample = rand(bulk_idxs[bulk_lookup[randname]])
        
        geochemdata = geochem_lookup[randname]
        errs = NamedTuple{Tuple(geochemkeys)}([abs(randn()*geochemdata[i].e) 
            for i in geochemkeys]
        )
        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((bulkzero[i][randsample]), errs[i])) for i in geochemkeys]
        )

        # Get EarthChem data for all types
        # majtype = unique([class_up(typelist, string(t)) for t in type])
        if :metased in type || :metaign in type
            type = (type..., :met)
        end

        bulksamples = falses(length(bulk_cats[1]))
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

        # Find match
        matches[i] = likelihood(EC.bulkage, macrostrat.age[i], EC.bulklat, EC.bulklon, 
            macrostrat.rocklat[i], macrostrat.rocklon[i], bulkgeochem, geochemdata, 
            EC.sampleidx
        )

        next!(p)
    end

    # Write data to a file
    writedlm("$matchedbulk_io", matches, '\t')

## --- End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).

    Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """

## --- End of File