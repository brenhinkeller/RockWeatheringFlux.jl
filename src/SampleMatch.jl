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


## --- Load Macrostrat data
    @info "Loading Macrostrat data ($macrostrat_io) $(Dates.format(now(), "HH:MM"))"
    fid = h5open("$macrostrat_io", "r")
    
    # Data
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"]),
        rocklon = read(fid["vars"]["rocklon"]),
        age = read(fid["vars"]["age"]),
        rocktype = read(fid["vars"]["rocktype"]),
        rockname = read(fid["vars"]["rockname"]),
        rockdescrip = read(fid["vars"]["rockdescrip"]),
    )
    
    # Rock type matches
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Rock name matches
    header = read(fid["type"]["name_cats_head"])
    data = read(fid["type"]["name_cats"])
    data = @. data > 0
    name_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    close(fid)
    

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

    # Rock type matches
    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Rock name matches
    rocknames = read(fid["bulktypes"]["bulk_lookup_head"])
    data = read(fid["bulktypes"]["bulk_lookup"])
    data = @. data > 0
    bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([data[:,i] for i in eachindex(rocknames)])

    close(fid)
    

## --- Alternatively, do the matching yourself
    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
    #     macrostrat.rockdescrip, unmultimatch=false, inclusive=false, source=:macrostrat
    # )

    # name_cats = match_rockname(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    # rocknames = string.(keys(name_cats))

    # bulk_cats = match_rocktype(bulktext.Rock_Name, bulktext.Type, bulktext.Material; 
    #     unmultimatch=false, inclusive=false, source=:earthchem
    # )
    
    # Get EarthChem samples for each rock name
    # typelist = get_rock_class(inclusive=true)      # Subtypes, major types include minors
    # classnames = string.(collect(keys(typelist)))
    # bulk_lookup = NamedTuple{keys(name_cats)}([falses(length(bulktext.Rock_Name)) 
    #     for _ in eachindex(name_cats)]
    # )

    # p = Progress(length(rocknames), desc="Finding EarthChem samples for each rock name")
    # for i in eachindex(rocknames)
    #     bulk_lookup[i] .= find_earthchem(rocknames[i], bulktext.Rock_Name, bulktext.Type, 
    #         bulktext.Material
    #     )

    #     # If no matches, jump up a class. Find everything within that class
    #     if count(bulk_lookup[i]) == 0
    #         searchlist = typelist[class_up(typelist, rocknames[i])]

    #         # Search all of those names; each class should at least have something
    #         for j in eachindex(searchlist)
    #             bulk_lookup[i] .|= find_earthchem(searchlist[j], bulktext.Rock_Name, 
    #                 bulktext.Type, bulktext.Material
    #             )
    #         end
    #     end
    #     next!(p)
    # end


## --- Get average geochemistry for each rock name
    geochem_lookup = NamedTuple{keys(name_cats)}([major_elements(bulk, bulk_lookup[i]) 
        for i in eachindex(bulk_lookup)]
    )


## --- Remove all multimatches from major types
    # This means that major types should be ONLY those samples which cannot be matched
    # with any minor types
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
    # Major types exclude minor types
    typelist = get_rock_class(inclusive=false)
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


## --- Remove all multimatches and major matches from Macrostrat rocks
    # Each sample can technically only be one rock type, and samples matched with major
    # types are technically a minor type (e.g. an igneous rock is either volcanic or 
    # plutonic).

    # Preallocate
    types = Array{Tuple{Symbol, Symbol}}(undef, length(macrostrat.age), 1)

    for i in eachindex(types)
        # Unweighted random selection of a rock name
        samplenames = get_type(name_cats, i, all_keys=true)

        
    end
    

    # Assign the sample the corresponding rock type. If the name maps to more than one
    # type, pick randomly. If the name maps to a major type, assign a minor type and 
    # corresponding rock name with a weighted-random selection, where weights are 
    # proportional to the abundance of that type / name in the Macrostrat samples

    
    


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
        
        # Pick a random sample to act as the geochemistry for that sample:
        # Replace all major types with a randomly selected minor type, proportional to that
        # minor type's abundance
        samplenames = get_type(name_cats, i, all_keys=true)
        randname, type = get_descriptive_name(samplenames, p_name, type, p_type, typelist, 
            minortypes
        )
        # If there aren't any matched samples in bulk_lookup, just move on
        # Not perfect but see if it helps?
        if count(bulk_lookup[randname])==0
            next!(p)
            continue
        end

        randsample = rand(bulk_idxs[bulk_lookup[randname]])
        
        geochemdata = geochem_lookup[randname]
        errs = NamedTuple{Tuple(geochemkeys)}([abs(randn()*geochemdata[i].e) 
            for i in geochemkeys]
        )
        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((bulkzero[i][randsample]), errs[i])) for i in geochemkeys]
        )

        # Get EarthChem data for all types
        # Metamorphic rock names don't give a lot of information about type, so just give
        # the matching algorithm everything
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