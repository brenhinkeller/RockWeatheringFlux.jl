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
    @info "Loading matched tuples $(Dates.format(now(), "HH:MM"))"
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

    # minorsed, minorign, minormet = get_minor_types()
    # for type in minorsed
    #     bulk_cats.sed .|= bulk_cats[type]
    # end
    # for type in minorign
    #     bulk_cats.ign .|= bulk_cats[type]
    # end
    # for type in minormet
    #     bulk_cats.met .|= bulk_cats[type]
    # end

    # EarthChem rock names
    rocknames = read(fid["vars"]["bulk_lookup_head"])
    data = read(fid["vars"]["bulk_lookup"])
    data = @. data > 0
    bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([data[:,i] for i in eachindex(rocknames)])

    close(fid)


## --- Load Macrostrat data
    @info "Loading Macrostrat data ($macrostrat_io) $(Dates.format(now(), "HH:MM"))"
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

    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
    #     macrostrat.rockdescrip, unmultimatch=false, inclusive=false, source=:macrostrat
    # )
    
## --- Load Earthchem bulk geochemical data
    @info "Loading EarthChem data $(Dates.format(now(), "HH:MM"))"
    bulkfid = h5open("output/bulk.h5", "r")

    # Bulk
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Bulk rock name, type, and material
    path = bulkfid["bulktext"]["sampledata"]
    header = read(path["header"])
    index = read(path["index"])

    target = ["Rock_Name", "Type", "Material"]
    targetind = [findall(==(i), header)[1] for i in target]

    bulktext = NamedTuple{Tuple(Symbol.(target))}(
        [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] 
            for i in eachindex(target)
        ]
    )
    close(bulkfid)

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


## --- Get the relative abundance of minor types in Macrostrat
    minorsed, minorign, minormet = get_minor_types()
    minortypes = (sed = minorsed, ign = minorign, met = minormet)
    pw = (
        sed = float.([count(macro_cats[i]) for i in minorsed]),
        ign = float.([count(macro_cats[i]) for i in minorign]),
        met = float.([count(macro_cats[i]) for i in minormet])
    )

    pw.sed ./= nansum(pw.sed)
    pw.ign ./= nansum(pw.ign)
    pw.met ./= nansum(pw.met)


## --- Find matching Earthchem sample for each Macrostrat sample
    # Preallocate
    matches = zeros(Int64, length(macro_cats.sed))

    # Define
    geochemkeys, = get_elements()               # Major elements
    bulk_idxs = collect(1:length(bulk.SiO2))    # Indices of bulk samples
    typelist = get_rock_class(false, false)     # Types, majors do not include minors

    # Zero-NaN version of the major elements in bulk
    bulkzero = deepcopy(bulk)
    bulkzero = NamedTuple{Tuple(geochemkeys)}(
        [zeronan!(bulkzero[i]) for i in geochemkeys]
    )

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

    @info "Starting sample matching $(Dates.format(now(), "HH:MM"))"

    p = Progress(length(matches), desc="Matching samples...")
    @timev for i in eachindex(matches)
        # Get the rock type and randomly select one sample rock name
        type = get_type(macro_cats, i, all_keys=true)
        if type==(:cover,) || type===nothing
            next!(p)
            continue
        end

        # Pick a random sample to act as the geochemistry for that sample
        samplenames = get_type(name_cats, i, all_keys=true)
        randname = rand(samplenames)
        randsample = rand(bulk_idxs[bulk_lookup[randname]])
        geochemdata = geochem_lookup[randname]

        errs = NamedTuple{Tuple(geochemkeys)}([abs(randn()*geochemdata[i].e) 
            for i in geochemkeys]
        )
        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((bulkzero[i][randsample]), errs[i])) for i in geochemkeys]
        )

        # Get EarthChem data for that type
        bulksamples = falses(length(bulk_cats[1]))
        type = replace_major(type, minortypes, pw)
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