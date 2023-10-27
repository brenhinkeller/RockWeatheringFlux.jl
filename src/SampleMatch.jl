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
    using LogExpFunctions

    # Local utilities
    include("utilities/Utilities.jl")

    # Start timer
    start = now()
    @info """ Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))
    Input: $macrostrat_io
    Output: $matchedbulk_io
    """


## --- Load Macrostrat data
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
    # For each rock name in this Tuple, get a BitVector of all samples that could
    # potentially match to that rock name. The rock names which match to sample (i) are
    # not indicative of the rock names which match to that sample.
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

    # Metamorphic rocks are inclusive, but a rock that is both metased and metaign is
    # re-assigned to just metamorphic
    for type in minormet
        macro_cats.met .&= .!macro_cats[type]
        bulk_cats.met .&= .!bulk_cats[type]     # Avoid metased / metaign cross contamination
    end

    # Add multimatches back in to metamorphic EarthChem samples
    # Metamorphic rocks with unknown protoliths are allowed to be metaseds and metaigns
    bulk_cats.metased .|= bulk_cats.met
    bulk_cats.metaign .|= bulk_cats.met
    bulk_cats.met .|= (bulk_cats.metased .& bulk_cats.metaign)


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


## --- Match each Macrostrat sample to a single informative rock name and type
    # Each sample can technically be only one rock type: matches to more than one type
    # are from grouping rocks together on geologic maps. 
    # 
    # Major classifications of sedimentary and igneous, and in some cases, metamorphic,
    # cannot be used to infer geochemical composition. Samples matched only with major
    # types are technically minor types: e.g., an igneous rock is either volcanic or 
    # plutonic, but this information is unknown for major-only matches.
    
    # If there are no minor types matched with the rock, randomly select a minor type and
    # associated rock name from those mapped to the major type. Select such that the
    # probability is directly proportional to the abundance of that type exposed on the
    # crust, determined by its abundance in the Macrostrat data.

    # Preallocate
    sampletypes = Array{Symbol}(undef, length(macrostrat.age), 1)
    samplenames = Array{Symbol}(undef, length(macrostrat.age), 1)

    # Metamorphic rock names without useful geochemical information
    uninformative_met = nondescriptive()

    # Major types should not include minor types, otherwise class_up will give majors
    typelist = get_rock_class(inclusive=false)
    minortypes = (sed = minorsed, ign=minorign, met=minormet)

    p = Progress(length(sampletypes) ÷ 10, desc="Sanitizing types...")
    for i in eachindex(sampletypes)
        # # Get names and types matched with the sample 
        # allnames = get_type(name_cats, i, all_keys=true)
        # alltypes = get_type(macro_cats, i, all_keys=true)

        # Unweighted random selection of a rock name
        allnames = get_type(name_cats, i, all_keys=true)
        if allnames===nothing
            sampletypes[i] = samplenames[i] = :none
            i%10==0 && next!(p)
            continue
        end
        s_name = rand(allnames)

        # Unweighted random selection of a corresponding type. If cover, pick again
        alltypes = class_up(typelist, string(s_name), all_types=true)
        notcover = .![t==:cover for t in alltypes]
        if count(notcover)==0
            sampletypes[i] = samplenames[i] = :none
            i%10==0 && next!(p)
            continue
        end
        s_type = rand(alltypes[notcover])

        # Re-assign major types
        if s_type==:ign || s_type==:sed || s_type==:met
            # Ign and sed types get weighted-random assignment to a minor type and name.
            s_type = minortypes[s_type][weighted_rand(p_type[s_type])]
            s_name = typelist[s_type][weighted_rand(p_name[s_type])]

        # elseif s_type==:met && string(s_name) in uninformative_met
        #     # Metamorphic rocks get re-assigned only if the name gives no useful information 
        #     # about its geochemistry

        #     # Otherwise, pick randomly
        #     s_type = minortypes.met[weighted_rand(p_type.met)]
        #     s_name = typelist[s_type][weighted_rand(p_name[s_type])]
        end
        
        # Assign
        sampletypes[i] = s_type
        samplenames[i] = Symbol(s_name)

        i%10==0 && next!(p)
    end


## --- Initialize for EarthChem sample matching
    # Definitions
    geochemkeys = get_elements()[1][1:end-1]        # Major non-volatile elements
    bulk_idxs = collect(1:length(bulk.SiO2))        # Indices of bulk samples

    # Zero-NaN version of the major elements in bulk
    bulkzero = deepcopy(bulk)
    bulkzero = NamedTuple{Tuple(geochemkeys)}([zeronan!(bulkzero[i]) for i in geochemkeys])


## --- Find matching EarthChem sample for each Macrostrat sample
    # As part of this process, we'll need to assume the geochemistry of the Macrostrat 
    # sample.
    # 
    # To preserve any multi-modal distributions in the data, we'll randomly 
    # pick one rock name matched with the sample, and randomly select one EarthChem sample
    # that was also matched with that rock name.
    # 
    # The error for each major element will be randomly sampled from a normal distribution 
    # with a mean and standard deviation equal to the mean and standard deviation for that
    # major element within the selected rock name.

    # Preallocate
    matches = zeros(Int64, length(macro_cats.sed))

    ismet = macro_cats.met .| macro_cats.metaign .| macro_cats.metased;

    @info "Starting sample matching $(Dates.format(now(), "HH:MM"))"
    p = Progress(length(matches), desc="Matching samples...")
    @timev for i in eachindex(matches)
        # if !ismet[i]
        #     next!(p)
        #     continue
        # end

        type = sampletypes[i]
        if type == :none
            next!(p)
            continue
        end

        # Pick a random EarthChem sample as the assumed geochem of the Macrostrat sample
        name = samplenames[i]
        randsample = rand(bulk_idxs[bulk_lookup[name]])
        geochemdata = geochem_lookup[name]
        errs = NamedTuple{Tuple(geochemkeys)}(
            # nanunzero!([abs(randn()*geochemdata[j].e) for j in geochemkeys], 1.0)
            nanunzero!([geochemdata[j].e for j in geochemkeys], 1.0)
        )
        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((bulkzero[j][randsample]), errs[j])) for j in geochemkeys]
        )

        # Get EarthChem data
        bulksamples = bulk_cats[type]
        EC = (
            bulklat = bulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = bulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = bulk.Age[bulksamples],                 # EarthChem age
            sampleidx = bulk_idxs[bulksamples],              # Indices of EarthChem samples
        )
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


## --- Plot
    # # All non-zero samples are metamorphic, so we don't have to restrict
    # t = @. matches > 0;
    # h1 = histogram(bulk.SiO2[matches[t]], bins=100, norm=:pdf,
    #     color=colors.met, lcolor=:match, 
    #     label="", xlabel="SiO2 [wt.%]", ylabel="Abundance", title="Metamorphic",
    #     framestyle=:box
    # )
    # ymin, ymax = ylims(h1)
    # ylims!(0, ymax*1.05)

    # t .&= macro_cats.metaign;
    # h2 = histogram(bulk.SiO2[matches[t]], bins=100, norm=:pdf,
    #     color=colors.metaign, lcolor=:match, 
    #     label="", xlabel="SiO2 [wt.%]", ylabel="Abundance", title="Metaigneous",
    #     framestyle=:box
    # )
    # ymin, ymax = ylims(h2)
    # ylims!(0, ymax*1.05)

    # h = Plots.plot(h1, h2, layout=(2, 1), size=(600,800), left_margin=(25,:px))
    # display(h)


## --- End of File