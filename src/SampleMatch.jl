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


## --- Remove misleading positives from the matches
    # This is important for getting an accurate representation of the abundance of each
    # rock type on Earth's surfaces. We're also doing this to the EarthChem rocks for 
    # consistency

    # If it's cover and something else, it can just be the something else
    allrocks = collect(keys(macro_cats))
    for type in allrocks
        type == :cover && continue
        macro_cats.cover .&= .!macro_cats[type]
        bulk_cats.cover .&= .!bulk_cats[type]
    end

    # If it IS just cover, it's not useful, so just take it out completely
    macro_cats.cover .= false
    bulk_cats.cover .= false

    # Don't match with metamorphic if you can help it: assume matches to other things
    # describes the protolith (at least, describes it better than guessing randomly)
    for type in allrocks
        type == :met && continue
        macro_cats.met .&= .!macro_cats[type]
        bulk_cats.met .&= .!bulk_cats[type]
    end

    # All granodiorites will also match with diorites, so take out those matches
    macro_cats.granodiorite .&= .!macro_cats.diorite
    bulk_cats.granodiorite .&= .!bulk_cats.diorite


## --- Calculate relative abundance of each type in the lithological map
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    # For volcanic / plutonic abundance volcanic and plutonic rocks MUST include subtypes...
    subminor_ign = (:volc, :plut, :carbonatite) 
    for type in minorvolc
        macro_cats.volc .|= macro_cats[type]
    end
    for type in minorplut
        macro_cats.plut .|= macro_cats[type]
    end

    # Count number of rocks in each subtype
    ptype = (
        sed = float.([count(macro_cats[i]) for i in minorsed]),
        volc = float.([count(macro_cats[i]) for i in minorvolc]),
        plut = float.([count(macro_cats[i]) for i in minorplut]),
        ign = float.([count(macro_cats[i]) for i in subminor_ign]),
    )
    
    # Calculate relative abundance / fraction
    ptype.sed ./= nansum(ptype.sed)
    ptype.volc ./= nansum(ptype.volc)
    ptype.plut ./= nansum(ptype.plut)
    ptype.ign ./= nansum(ptype.ign)

    # I want to calculate the relative abundance of protoliths that could get turned into
    # metamorphic rocks. Metamorphic rocks (with no known protolith) could be.... 
    # anything? Or more technically, not anything. It's probably not from a chert 
    # protolith. Metacarbonates (or metacarbonatites??) are also probably not defined as 
    # a gneiss. So exclude carbonates, evaporites, chert, phosphorite, coal, and carbonatites
    protolith = (:siliciclast, :shale, minorvolc..., minorplut...,)
    pprotolith = float.([count(macro_cats[i]) for i in protolith])
    pprotolith ./= nansum(pprotolith)


## --- If samples are matched to a rock subtype and a rock major type, don't
    # This is mostly important for figuring out what minor type to assign each sample

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


## --- Match each Macrostrat sample to a single informative rock name and type
    # Doing this in several passes over the sample set means that I can optimize sections
    # that can be optimized, which will make this faster... by two orders of magnitude lol

    # Preallocate
    bigtypes = Array{Symbol}(undef, length(macrostrat.age), 1)      # Sed, volc, ign, etc.
    littletypes = Array{Symbol}(undef, length(macrostrat.age), 1)   # Shale, chert, etc.

    # Pass one: randomly pick a type for each sample
    for i in eachindex(littletypes)
        alltypes = get_type(macro_cats, i, all_keys=true)
        littletypes[i] = rand(alltypes)
    end

    # Pass two: reassign major types to a minor subtype
    for i in eachindex(littletypes)
        if littletypes[i] == :met
        # Metamorphic: assign a protolith
            littletypes[i] = protolith[weighted_rand(pprotolith)]
            continue
        elseif littletypes[i] == :sed
        # Sedimentary: assign a minor type
            littletypes[i] = minorsed[weighted_rand(ptype.sed)]
            continue
        elseif littletypes[i] == :ign
        # Igneous: assign volcanic / plutonic / carbonatite
            littletypes[i] = subminor_ign[weighted_rand(ptype.ign)]
        end

        # Volcanic / plutonic samples should be assigned an appropriate subtype
        if littletypes[i] == :volc 
            littletypes[i] = minorvolc[weighted_rand(ptype.volc)]
        elseif littletypes[i] == :plut 
            littletypes[i] = minorplut[weighted_rand(ptype.plut)]
        end
    end

    # Pass three: figure out if the type is a sedimentary or igneous type
    for i in eachindex(bigtypes)
        if littletypes[i]==:none 
            bigtypes[i] = :none 
        else
            bigtypes[i] = class_up(littletypes[i], minorsed, minorvolc, minorplut, minorign)
        end
    end


## --- Initialize for EarthChem sample matching
    # Definitions
    geochemkeys = get_elements()[1][1:end-1]        # Major non-volatile elements
    bulk_inds = collect(1:length(bulk.SiO2))        # Indices of bulk samples

    # # Zero-NaN version of the major elements in bulk
    bulkzero = NamedTuple{Tuple(geochemkeys)}([zeronan(bulk[i]) for i in geochemkeys])

    # Average geochemistry of each rock type 
    realrocks = deleteat!(allrocks, findall(x->x==:cover,allrocks))     # That isn't cover
    geochem_lookup = NamedTuple{Tuple(realrocks)}([major_elements(bulk, bulk_cats[i])
        for i in eachindex(realrocks)]
    );

    # Make major types inclusive of minor types?
    for type in minorsed
        macro_cats.sed .|= macro_cats[type]
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorvolc
        macro_cats.volc .|= macro_cats[type]
        bulk_cats.volc .|= bulk_cats[type]
    end
    for type in minorplut
        macro_cats.plut .|= macro_cats[type]
        bulk_cats.plut .|= bulk_cats[type]
    end
    for type in minorign
        macro_cats.ign .|= macro_cats[type]
        bulk_cats.ign .|= bulk_cats[type]
    end


## --- Find matching EarthChem sample for each Macrostrat sample
    # As part of this process, we'll need to assume the geochemistry of the Macrostrat 
    # sample.
    # 
    # To preserve any multi-modal distributions in the data, we'll randomly 
    # pick one rock name matched with the sample, and randomly select one EarthChem sample
    # that was also matched with that rock name.

    # Preallocate
    matches = zeros(Int64, length(macro_cats.sed))

    @info "Starting sample matching $(Dates.format(now(), "HH:MM"))"
    p = Progress(length(matches), desc="Matching samples...")
    @timev for i in eachindex(matches)
        ltype = littletypes[i] 
        btype = bigtypes[i]
        if ltype == :none
            next!(p)
            continue
        end

        # Assume the geochemical composition of the lithological sample: pick randomly
        randsample = rand(bulk_inds[bulk_cats[ltype]])
        uncert = nanunzero!([geochem_lookup[ltype][elem].e for elem in geochemkeys], 1.0)
        values = [bulkzero[elem][randsample] for elem in geochemkeys]

        geochemdata = NamedTuple{Tuple(geochemkeys)}(
            NamedTuple{(:m, :e)}((values[j], uncert[j])) for j in eachindex(values)
        )

        # Get EarthChem data
        bulksamples = bulk_cats[btype]
        EC = (
            bulklat = bulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = bulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = bulk.Age[bulksamples],                 # EarthChem age
            sampleinds = bulk_inds[bulksamples],             # Indices of EarthChem samples
        )
        bulkgeochem = NamedTuple{Tuple(geochemkeys)}([bulkzero[i][bulksamples] 
            for i in geochemkeys]
        )

        # Find match
        matches[i] = likelihood(EC.bulkage, macrostrat.age[i], EC.bulklat, EC.bulklon, 
            macrostrat.rocklat[i], macrostrat.rocklon[i], bulkgeochem, geochemdata, 
            EC.sampleinds
        )

        next!(p)
    end

    # Write data to a file
    writedlm("$matchedbulk_io", [matches string.(littletypes)], '\t')

    # End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """

## --- Quick visualization of distributions
    using Plots

    # Make matches nice and inclusive
    for type in minorsed
        macro_cats.sed .|= macro_cats[type]
    end
    for type in minorvolc
        macro_cats.volc .|= macro_cats[type]
    end
    for type in minorplut
        macro_cats.plut .|= macro_cats[type]
    end
    for type in minorign
        macro_cats.ign .|= macro_cats[type]
    end

    # Matched silica 
    t = matches .> 0;
    silica = bulk.SiO2[matches[t]]

    # Send this to three of your friends to totally visualize them
    get_visualized = (:volc, :plut, :ign, :siliciclast, :shale, :carb, :chert, :sed)
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(get_visualized)) 
    for i in eachindex(get_visualized)
        # Figure out what kind of bins we want
        r = get_visualized[i]
        if r==:sed || r in minorsed
            nb = (0,100,100)
        elseif r==:ign || r in minorign
            nb = (40,80,40)
        end

        # Lol get visualized!
        c, n = bincounts(silica[macro_cats[r][t]], nb...)
        n = float(n) ./ nansum(float(n) .* step(c))
        h = plot(c, n, seriestype=:bar, framestyle=:box, color=colors[r], linecolor=colors[r],
            title="$(string(r)); n = $(count(macro_cats[r]))", label="",     
            ylims=(0, round(maximum(n), digits=2)+0.01) 
        )
        fig[i] = h
    end

    h = plot(fig..., layout=(3,3), size=(1800, 1200))
    display(h)


## --- End of File