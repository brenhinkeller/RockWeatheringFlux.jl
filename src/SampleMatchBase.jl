## --- Load Macrostrat data
    fid = h5open(filemacrostrat, "r")
    
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
    # header = read(fid["type"]["macro_cats_head"])
    # data = read(fid["type"]["macro_cats"])
    # data = @. data > 0
    # macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
        macrostrat.rockdescrip, showprogress=show_progress
    )
    close(fid)


## --- Load Earthchem bulk geochemical data
    fid = h5open(filebulk, "r")

    # Bulk
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Rock type matches
    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    close(fid)


## --- Update matches in mapped samples and geochemical samples
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    # I hate cover
    macro_cats = delete_cover(macro_cats)
    bulk_cats = delete_cover(bulk_cats)

    # Metamorphic rocks are only metamorphic if we cannot infer a protolith
    for type in keys(macro_cats)
        type==:met && continue
        macro_cats.met .&= .!macro_cats[type]
        bulk_cats.met .&= .!bulk_cats[type]
    end


## --- Calculate relative abundance of each type in the lithological map
    # To calculate total volcanic / plutonic abundance, must include subtypes
    include_minor!(macro_cats)

    # Absolute abundance (count) of each rock type
    ptype = (;
        sed = float.([count(macro_cats[i]) for i in minorsed]),
        volc = float.([count(macro_cats[i]) for i in minorvolc]),
        plut = float.([count(macro_cats[i]) for i in minorplut]),
        ign = float.([count(macro_cats[i]) for i in minorign]),
    )

    # Calculate fractional abundance (fraction of total)
    ptype.sed ./= nansum(ptype.sed)
    ptype.volc ./= nansum(ptype.volc)
    ptype.plut ./= nansum(ptype.plut)
    ptype.ign ./= nansum(ptype.ign)

    # Calculate the relative abundance of protoliths that could become metamorphic rocks.
    # Probably not chert / evaporite / coal / carbonatite? Metacarbonates are unlikely to
    # end up as the type of rocks we have as uncategorized metamorphic (e.g., gneiss)
    protolith = (:siliciclast, :shale, :volc, :plut)
    p_protolith = float.([count(macro_cats[i]) for i in protolith])
    p_protolith ./= nansum(p_protolith)


    # If samples are matched to a rock subtype and a rock major type, don't!
    # This is mostly important for figuring out what minor type to assign each sample
    # If we know what kind of rock we have... we don't want to lose that information
    exclude_minor!(macro_cats)
    exclude_minor!(bulk_cats)


## --- Match each Macrostrat sample to a single informative rock name and type
    # Doing this in several passes over the sample set means that I can optimize sections
    # that can be optimized, which will make this faster... by two orders of magnitude
    # I love coding. Affirm: I AM optimization

    # Preallocate
    littletypes = Array{Symbol}(undef, length(macrostrat.age), 1)   # Shale, chert, etc.

    # Pass one: randomly pick a type for each sample
    for i in eachindex(littletypes)
        alltypes = get_type(macro_cats, i, all_keys=true)
        littletypes[i] = rand(alltypes)
    end

    # Pass two: assign uncategorized metamorphic rocks to a protolith 
    for i in eachindex(littletypes)
        if littletypes[i] == :met 
            littletypes[i] = protolith[weighted_rand(p_protolith)]
        end
    end

    # Pass three: reassign major types to a minor subtype
    for i in eachindex(littletypes)
        if littletypes[i] == :sed
            littletypes[i] = minorsed[weighted_rand(ptype.sed)]

        elseif littletypes[i] == :volc 
            littletypes[i] = minorvolc[weighted_rand(ptype.volc)]

        elseif littletypes[i] == :plut 
            littletypes[i] = minorplut[weighted_rand(ptype.plut)]

        elseif littletypes[i] == :ign 
            # Pick a sub-class (volcanic / plutonic / carbonatite) and re-assign volc / plut
            littletypes[i] = minorign[weighted_rand(ptype.ign)]

            if littletypes[i] == :volc
                littletypes[i] = minorvolc[weighted_rand(ptype.volc)]

            elseif littletypes[i] == :plut 
                littletypes[i] = minorplut[weighted_rand(ptype.plut)]
                
            end

        end
    end


## --- Initialize for EarthChem sample matching
    # Definitions
    geochemkeys = get_elements()[1][1:end-1]        # Major non-volatile elements
    bulk_inds = collect(1:length(bulk.SiO2))        # Indices of bulk samples

    # # Zero-NaN version of the major elements in bulk
    bulkzero = NamedTuple{Tuple(geochemkeys)}([zeronan(bulk[i]) for i in geochemkeys]);

    # Average geochemistry of each rock type
    geochem_lookup = NamedTuple{Tuple(keys(macro_cats))}([major_elements(bulk, bulk_cats[i])
        for i in eachindex(keys(macro_cats))]
    );

    # # Re-include minor types (this only matters if we don't reassign volc / plut)
    # include_minor!(bulk_cats)
    # include_minor!(macro_cats)


## --- Find matching EarthChem sample for each Macrostrat sample
    # Preallocate
    matches = zeros(Int64, length(macro_cats.sed))

    # @info "Starting sample matching $(Dates.format(now(), "HH:MM"))"
    p = Progress(length(matches), desc="Matching samples...", enabled=show_progress)
    @time for i in eachindex(matches)
        ltype = littletypes[i] 

        # Skip if no type assigned or no samples present
        if ltype == :none || isempty(bulk_inds[bulk_cats[ltype]]) || isempty(geochem_lookup[ltype])
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
        bulksamples = bulk_cats[ltype]
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


## --- End of File