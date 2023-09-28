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
    # class_names = string.(collect(keys(typelist)))
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
    
    # Minor rock types
    p_type = (
        sed = float.([count(macro_cats[i]) for i in minorsed]),
        ign = float.([count(macro_cats[i]) for i in minorign]),
        met = float.([count(macro_cats[i]) for i in minormet])
    )
    [p_type[i] ./= sum(p_type[i]) for i in keys(p_type)]

    # Descriptive rock names
    minortypes = (minorsed..., minorign..., minormet...)
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
    sampletypes = Array{Symbol}(undef, length(macrostrat.age), 1)
    samplenames = Array{Symbol}(undef, length(macrostrat.age), 1)

    # Metamorphic rock names without useful geochemical information
    ambig_met = nondescriptive()

    # Major types should not include minor types, otherwise class_up will give majors
    typelist = get_rock_class(inclusive=false)
    minortypes = (sed = minorsed, ign=minorign, met=minormet)

    p = Progress(length(sampletypes) ÷ 10, desc="Sanitizing types...")
    for i in eachindex(sampletypes)
        # Unweighted random selection of a rock name
        allnames = get_type(name_cats, i, all_keys=true)
        if allnames===nothing
            sampletypes[i] = samplenames[i] = :none
            i%10==0 && next!(p)
            continue
        end
        s_name = rand(allnames)

        # Unweighted random selection of a corresponding type. If the type is cover, pick again
        alltypes = class_up(typelist, string(s_name), all_types=true)
        notcover = .![t==:cover for t in alltypes]
        if count(notcover)==0
            sampletypes[i] = samplenames[i] = :none
            i%10==0 && next!(p)
            continue
        end
        s_type = rand(alltypes[notcover])

        # Ign and sed types get weighted-random assignment to a minor type and name.
        # Metamorphic rocks get re-assigned only if the name gives no useful information 
        # about its geochemistry
        if s_type==:ign || s_type==:sed
            s_type = minortypes[s_type][weighted_rand(p_type[s_type])]
            s_name = typelist[s_type][weighted_rand(p_name[s_type])]
        elseif s_type==:met && string(s_name) in ambig_met
            s_type = minortypes.met[weighted_rand(p_type.met)]
            s_name = typelist[s_type][weighted_rand(p_name[s_type])]
        end
        
        # Assign
        sampletypes[i] = s_type
        samplenames[i] = Symbol(s_name)

        i%10==0 && next!(p)
    end


## --- Initialize for EarthChem sample matching
    # Definitions
    geochemkeys, = get_elements()                   # Major elements
    bulk_idxs = collect(1:length(bulk.SiO2))        # Indices of bulk samples

    # Zero-NaN version of the major elements in bulk
    bulkzero = deepcopy(bulk)
    bulkzero = NamedTuple{Tuple(geochemkeys)}([zeronan!(bulkzero[i]) for i in geochemkeys])


## --- TEMPORARY: test matching algorithm for carbonates only
    using Plots
    
    filter = findall(x -> x==:carb, sampletypes)
    matches = zeros(Int64, length(filter))

    p = Progress(length(matches)÷10, desc="Matching samples...")
    for i in eachindex(matches)
        # Pick a random sample as the assumed geochemistry
        name = samplenames[filter[i]]
        geochemdata = geochem_lookup[name]
        errs = NamedTuple{Tuple(geochemkeys)}(
            zeronan!([abs(randn()*geochemdata[i].e) for i in geochemkeys])
        )

        randsample = rand(bulk_idxs[bulk_lookup[name]])
        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((bulkzero[i][randsample]), errs[i])) for i in geochemkeys]
        )

        # Get the EarthChem data
        bulksamples = bulk_cats.carb
        EC = (
            bulklat = bulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = bulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = bulk.Age[bulksamples],                 # EarthChem age
            sampleidx = bulk_idxs[bulksamples],              # Indices of EarthChem samples
        )
        bulkgeochem = NamedTuple{Tuple(geochemkeys)}([bulkzero[i][bulksamples] 
            for i in geochemkeys]
        )

        # Match
        matches[i] = likelihood(EC.bulkage, macrostrat.age[i], EC.bulklat, EC.bulklon, 
            macrostrat.rocklat[i], macrostrat.rocklon[i], bulkgeochem, geochemdata, 
            EC.sampleidx
        )

        i$10==0 && next!(p)
    end

    # Plot results
    bins = (0, 100, 100)

    c, n = bincounts(bulk.SiO2[matches], bins...)
    n = float(n) ./ nansum(float(n) .* step(c))
    SiO2ₕ = plot(c, n, seriestype=:bar, framestyle=:box, label="", ylabel="Weight", 
        xlabel="SiO2 [wt.%]", ylims=(0, round(maximum(n), digits=2)+0.01) 
    )

    c, n = bincounts(bulk.CaO[matches], bins...)
    n = float(n) ./ nansum(float(n) .* step(c))
    CaOₕ = plot(c, n, seriestype=:bar, framestyle=:box, label="", 
        xlabel="CaO [wt.%]", ylims=(0, round(maximum(n), digits=2)+0.01) 
    )

    c, n = bincounts(bulk.CaCO3[matches], bins...)
    n = float(n) ./ nansum(float(n) .* step(c))
    CaCO3ₕ = plot(c, n, seriestype=:bar, framestyle=:box, label="",
        xlabel="CaCO₃ [wt.%]", ylims=(0, round(maximum(n), digits=2)+0.01) 
    )

    h = plot(SiO2ₕ, CaOₕ, CaCO3ₕ, layout=(1, 3), size=(1800, 400),
        bottom_margin=(50, :px), left_margin=(50, :px)
    )
    display(h)

    # How many matches are the same sample?
    c, n = bincounts(bulk.SiO2[matches], bins...)
    i = findmax(n)[2]                                   # Index of largest bin
    s = step(c)/2                                       # Half-step to get bin edges
    tᵢ = @. c[i]-s <= bulk.SiO2[matches] <= c[i]+s      # Filter for samples in that bin
    ind = matches[tᵢ]                                   # Indices of the bin

    j, f = modal(ind)                                   # Value and frequency of modal index
    j₁, f₁ = modal(matches)                             # Value and frequency of modal sample

    @info """ Frequency of single sample in:
    Modal bin: $(round(f/length(ind)*100, sigdigits=3))%
    All matches: $(round(f₁/length(matches)*100, sigdigits=3))%
    """


## --- Find matching EarthChem sample for each Macrostrat sample
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

    @info "Starting sample matching $(Dates.format(now(), "HH:MM"))"
    p = Progress(length(matches), desc="Matching samples...")
    @timev for i in eachindex(matches)
        type = sampletypes[i]
        name = samplenames[i]
        if type==:none 
            next!(p)
            continue
        end

        # Pick a random sample to act as the geochemistry for that sample:
        randsample = rand(bulk_idxs[bulk_lookup[name]])
        
        geochemdata = geochem_lookup[name]
        errs = NamedTuple{Tuple(geochemkeys)}(
            zeronan!([abs(randn()*geochemdata[i].e) for i in geochemkeys])
        )
        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((bulkzero[i][randsample]), errs[i])) for i in geochemkeys]
        )

        # Get EarthChem data. EarthChem major types are not inclusive of minor types, but
        # all metamorphic rocks are sent for every metamorphic rock type
        # Would it maybe be easier just to not make metamorphic rocks exclusive wrt types?
        bulksamples = bulk_cats[type]

        if type==:metased || type==:metaign
            bulksamples .|= bulk_cats.met
        elseif type==:met 
            bulksamples .|= bulk_cats.metased
            bulksamples .|= bulk_cats.metaign
        end

        EC = (
            bulklat = bulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = bulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = bulk.Age[bulksamples],                 # EarthChem age
            sampleidx = bulk_idxs[bulksamples],              # Indices of EarthChem samples
        )

        # Get all EarthChem samples for that rock type
        bulkgeochem = NamedTuple{Tuple(geochemkeys)}([bulkzero[i][bulksamples] 
            for i in geochemkeys]
        )

        # Find match
        matches[i] = likelihood(EC.bulkage, macrostrat.age[i], EC.bulklat, EC.bulklon, 
            macrostrat.rocklat[i], macrostrat.rocklon[i], bulkgeochem, geochemdata, 
            EC.sampleidx
        )

        next!(p)
    end

    # Write data to a file
    writedlm("$matchedbulk_io", [matches string.(sampletypes)], '\t')

    # End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).

    Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """


## --- End of File