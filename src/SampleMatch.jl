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


## --- Do terrible things to the lists of matched samples
    # If I get the types that are matched to a given sample, I want them to be as specific
    # as possible so I can characterize the sample as accurately as possible.
    
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    allrocks = collect(keys(typelist))

    # Don't match with a major type if you can match with a minor type. Also, don't match
    # to volcanic or plutonic if you can do better.
    for type in minorsed
        macro_cats.sed .&= .!macro_cats[type]
        bulk_cats.sed .&= .!bulk_cats[type]
    end
    for type in minorvolc
        macro_cats.volc .&= .!macro_cats[type]
        bulk_cats.volc .&= .!bulk_cats[type]
    end
    for type in minorplut
        macro_cats.plut .&= .!macro_cats[type]
        bulk_cats.plut .&= .!bulk_cats[type]
    end
    for type in minorign
        macro_cats.ign .&= .!macro_cats[type]
        bulk_cats.ign .&= .!bulk_cats[type]
    end

    # Don't match with metamorphic if you can help it: assume matches to other things
    # describes the protolith (at least, describes it better than guessing randomly)
    for type in allrocks
        type == :met && continue
        macro_cats.met .&= .!macro_cats[type]
        bulk_cats.met .&= .!bulk_cats[type]
    end

    # If it's cover and something else, it can just be the something else
    for type in allrocks
        type == :cover && continue
        macro_cats.cover .&= .!macro_cats[type]
        bulk_cats.cover .&= .!bulk_cats[type]
    end

    # If it IS just cover, it's not useful, so just take it out completely
    macro_cats.cover .= false
    bulk_cats.cover .= false

    # All granodiorites will also match with diorites, so take out those matches
    macro_cats.granodiorite .&= .!macro_cats.diorite
    bulk_cats.granodiorite .&= .!bulk_cats.diorite


## --- Calculate relative abundance of each type in the lithological map
    subminor_ign = (:volc, :plut, :carbonatite)     # Volc and plut MUST include subtypes

    # Counts
    psed = float.([count(macro_cats[i]) for i in minorsed])
    pvolc = float.([count(macro_cats[i]) for i in minorvolc])
    pplut = float.([count(macro_cats[i]) for i in minorplut])
    pign = float.([count(macro_cats[i]) for i in subminor_ign])
    
    # Relative abundance / fraction
    psed ./= nansum(psed)
    pvolc ./= nansum(pvolc)
    pplut ./= nansum(pplut)
    pign ./= nansum(pign)


## --- Match each Macrostrat sample to a single informative rock name and type
    # # Metamorphic rocks (with no known protolith) could be.... anything? Sure, I guess so
    # # Technically, not anything. It's probably not from a chert protolith. Metacarbonates
    # # (or metacarbonatites??) are also probably not defined as a gneiss.
    # # So exclude carbonates, evaporites, chert, phosphorite, coal, and carbonatites
    # protolith = (:siliciclast, :shale, :sed, minorvolc..., minorplut..., :ign);
    # for type in protolith
    #     macro_cats.met .|= macro_cats[type]
    #     bulk_cats.met .|= bulk_cats[type]
    # end

    # Preallocate
    bigtypes = Array{Symbol}(undef, length(macrostrat.age), 1)
    littletypes = Array{Symbol}(undef, length(macrostrat.age), 1)

    p = Progress(length(bigtypes) ÷ 10, desc="Sanitizing types...")
    for i in eachindex(bigtypes)
        alltypes = get_type(macro_cats, i, all_keys=true)

        if alltypes===nothing
            bigtypes[i] = littletypes[i] = :none
            i%10==0 && next!(p)
            continue
        end

        # Randomly select a rock name, unless there aren't anything
        samplename = rand(alltypes)

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


## --- Load resampled data
    # # Load data
    # fid = h5open("output/resampled/resampled_rocknames.h5", "r")
    # header = read(fid["vars"]["header"])
    # data = read(fid["vars"]["simout"])
    # simout = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # # Get rocknames from kz number
    # kz_names = read(fid["vars"]["rocknames_kz"])
    # kz_i = collect(1:length(kz_names))
    # simout_names = NamedTuple{Tuple(Symbol.(kz_names))}([falses(length(simout.Kz_1)) 
    #     for _ in eachindex(kz_names)]
    # );
    
    # kz_out = (:Kz_1, :Kz_2, :Kz_3, :Kz_4, :Kz_5, :Kz_6)
    # for i in kz_i
    #     for k in kz_out
    #         # Get all samples that match with the given rock name / kz number 
    #         t = @. simout[k] == i
    #         simout_names[Symbol(kz_names[i])][t] .|= true
    #     end
    # end

    # # For names with no matches, jump up a class and match with all mapped names
    # # Typelist is inclusive, because if something is an unknown sed, we want it to match
    # # with every sedimentary rock
    # typelist = get_rock_class(inclusive=true)
    # for k in keys(simout_names)
    #     if count(simout_names[k]) < 3
    #         upper = class_up(typelist, string(k))
    #         for r in typelist[upper]
    #             simout_names[k] .|= simout_names[Symbol(r)]
    #         end
    #         # println(k)
    #     end
    # end

    # # For major and minor subtypes, match with all names for that subtype 
    # for k in keys(bulk_cats)
    #     if haskey(simout_names, k)
    #         for r in typelist[k]
    #             if haskey(simout_names, Symbol(r))
    #                 simout_names[k] .|= simout_names[Symbol(r)]
    #             end
    #         end
    #     end
    # end


## --- Initialize for EarthChem sample matching
    # Definitions
    geochemkeys = get_elements()[1][1:end-1]        # Major non-volatile elements
    bulk_idxs = collect(1:length(bulk.SiO2))        # Indices of bulk samples
    simout_ind = collect(1:length(simout.SiO2))     # Indices of resampled data

    # # Zero-NaN version of the major elements in bulk
    bulkzero = NamedTuple{Tuple(geochemkeys)}([zeronan(bulkzero[i]) for i in geochemkeys])

    # Zero-NaN version of the major elements in the resampled dataset
    simout_zeronan = NamedTuple{Tuple(geochemkeys)}(zeronan(simout[i]) for i in geochemkeys)

    # # Get average geochemistry for each rock name
    # geochem_lookup = NamedTuple{keys(name_cats)}([major_elements(bulk, bulk_lookup[i]) 
    #     for i in eachindex(bulk_lookup)]
    # )

    # Average geochemistry for each rock name
    simout_geochem = NamedTuple{keys(simout_names)}(
        [major_elements(simout, simout_names[i]) for i in eachindex(simout_names)]
    );


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

    # ismet = macro_cats.met .| macro_cats.metaign .| macro_cats.metased;

    @info "Starting sample matching $(Dates.format(now(), "HH:MM"))"
    p = Progress(length(matches), desc="Matching samples...")
    @timev for i in eachindex(matches)
        # if !ismet[i]
        #     next!(p)
        #     continue
        # end

        type = sampletypes[i]
        # if type == :none
        #     next!(p)
        #     continue
        # end

        if type != :shale
            next!(p)
            continue
        end

        # # Pick a random EarthChem sample as the assumed geochem of the Macrostrat sample
        # name = samplenames[i]
        # randsample = rand(bulk_idxs[bulk_lookup[name]])
        # geochemdata = geochem_lookup[name]
        # errs = NamedTuple{Tuple(geochemkeys)}(
        #     nanunzero!([geochemdata[j].e for j in geochemkeys], 1.0)
        # )
        # geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
        #     tuple.((bulkzero[j][randsample]), errs[j])) for j in geochemkeys]
        # )

        # Pick a random resampled Earthchem sample as the assumed geochemistry of the 
        # Burwell sample
        name = samplenames[i]
        randsample = rand(simout_ind[simout_names[name]])
        uncertainty = NamedTuple{Tuple(geochemkeys)}(
            nanunzero!([simout_geochem[name][j].e for j in geochemkeys], 1.0)
        )
        geochemdata = NamedTuple{Tuple(geochemkeys)}([NamedTuple{(:m, :e)}(
            tuple.((simout_zeronan[j][randsample]), uncertainty[j])) for j in geochemkeys]
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

    # # Write data to a file
    # writedlm("$matchedbulk_io", [matches string.(sampletypes)], '\t')

    # # End timer
    # stop = now()
    # @info """
    # Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    # Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    # """

## --- Plot
    t = @. matches > 0;
    h = histogram(bulk.SiO2[matches[t]], bins=100, norm=:pdf,
        color=colors.shale, lcolor=:match, 
        label="Shale", xlabel="SiO₂ [wt.%]", ylabel="Weight", framestyle=:box
    )
    ymin, ymax = ylims(h)
    ylims!(0, ymax*1.05)
    display(h)

    # Would you still have unexpected modes if I was a snail :(
    snails = Symbol.(typelist.shale)
    snailfig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(snails))
    for i in eachindex(snails)
        # How many Burwell samples do we have, and what are they matched to?
        n = count(name_cats[snails[i]])
        samples = matches[t .& name_cats[snails[i]]]

        # Plot
        h = histogram(bulk.SiO2[samples], bins=100, norm=:pdf,
            color=colors.shale, lcolor=:match, framestyle=:box,
            title="$(snails[i]); n=$n", xlabel="SiO₂ [wt.%]", ylabel="Weight", label="",
            xlims=(0,100)
        )
        ymin, ymax = ylims(h)
        ylims!(0, ymax*1.05)
        # display(h)
        snailfig[i] = h

        # Get the ratio of number of samples to number of matched samples
        s = t .& name_cats[snails[i]]
        @info """ $(snails[i]):
        n burwell samples  = $n
        n possible samples = $(count(simout_names[snails[i]]))
        
        matches:
        n unique samples   = $(length(unique(samples)))
        n total samples    = $(count(s))
        """
    end

    nrows = round(Int, length(snails)/2)
    h = Plots.plot(snailfig..., layout=(nrows, 2), size=(1200,nrows*400),
        xlabel="", ylabel="",
        left_margin=(25,:px),
        legendfont=15
    )
    display(h)


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