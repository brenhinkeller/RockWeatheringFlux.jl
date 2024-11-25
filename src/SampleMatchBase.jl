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
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];

    # Volcaniclastic is a cursed category and I want it OUT
    i = findall(x->x==:volcaniclast, minorvolc)
    minorvolc = minorvolc[1:end .!= i]

    # I hate cover
    macro_cats = delete_cover(macro_cats)
    bulk_cats = delete_cover(bulk_cats)


## --- Calculate relative abundance of each descriptive rocktype
    # For random re-assignment of nondescriptive types to descriptive types!
    # To calculate accurate volcanic / plutonic abundances, include subtypes
    include_minor!(macro_cats)

    # Calculate relative abundance of each descriptive type. Each type sums to 1.0
    # Except volcaniclastic is a cursed category and I don't want it 
    
    type_abundance = (;
        sed = float.([count(macro_cats[i]) for i in minorsed]),
        volc = float.([count(macro_cats[i]) for i in minorvolc]),
        plut = float.([count(macro_cats[i]) for i in minorplut]),
        ign = float.([count(macro_cats[i]) for i in minorign]),
    )
    type_abundance.sed ./= nansum(type_abundance.sed)
    type_abundance.volc ./= nansum(type_abundance.volc)
    type_abundance.plut ./= nansum(type_abundance.plut)
    type_abundance.ign ./= nansum(type_abundance.ign)

    # Calculate relative abundance of possible undifferentiated protoliths.
    # Carbonates, evaporites, chert, phosphorites, and coal are unlikely to become 
    # undifferentiated metamorphic rocks (schist, gneiss, migmatite, etc.) 
    protolith = (:siliciclast, :shale, :volc, :plut)
    protolith_abundance = float.([count(macro_cats[i]) for i in protolith])
    protolith_abundance ./= nansum(protolith_abundance)

    # And metasedimentary protoliths 
    protolith_metased = (:siliciclast, :shale)
    protolith_metased_abundance = float.([count(macro_cats[i]) for i in protolith_metased])
    protolith_metased_abundance ./= nansum(protolith_metased_abundance)

    
## --- Deal with weird shit 
    # For undifferentiated metamorphic rocks: 
    # If we have information about descriptive rock types, then take that as the protolith 
    # and ignore the metamorphic part of it. This doesn't matter as much for metaigneous rocks, 
    # but we don't want to take "sed" from metasedimentary rocks and then assign our migmatite 
    # to an evaporite!
    # TL;DR For sedimentary rocks, if we know there's an associated descriptive type, stick with
    # that. Save undifferentiated sed + met as metasedimentary rocks. If it's igneous, stick with 
    # that
    for type in minorsed 
        macro_cats.met .&= .!macro_cats[type]
    end
    metased = macro_cats.met .& macro_cats.sed;

    for type in (:ign, minorign..., minorvolc..., minorplut...) 
        macro_cats.met .&= .!macro_cats[type]
    end
    
    # For volcaniclastic rocks: 
    # The only terms to match as volcaniclastic could be LITERALLY anything (it's like... 
    # lahar. ash that deposited in a swamp.) So if we have better information about what 
    # this is... use that instead. We also don't want to match with volcaniclastic rocks. 
    # As an aside -- we kinda... know what all the volcaniclastic rocks are. This wipes out 
    # the entire category
    for type in (:ign, minorign..., minorvolc..., minorplut...) 
        macro_cats.volcaniclast .&= .!macro_cats[type]
    end


## --- Match each Macrostrat sample to ONE informative (descriptive) rock name and type
    # Doing this in several passes over the sample set means that I can optimize sections
    # that can be optimized, which will make this faster... by two orders of magnitude
    # I love coding. Affirm: I AM optimization

    # Preallocate
    littletypes = Array{Symbol}(undef, length(macrostrat.age), 1)   # Descriptive types

    # Make sure we're only looking at descriptive types -- if we know something's a shale, 
    # we don't want to see "sed" pop up and thing it needs to be reassigned
    exclude_minor!(macro_cats)
    exclude_minor!(bulk_cats)

    # Pass one: randomly pick a type for each sample
    for i in eachindex(littletypes)
        alltypes = get_type(macro_cats, i, all_keys=true)
        littletypes[i] = rand(alltypes)
    end

    # Pass two: assign weird shit (undifferentiated metamorphic + undifferentiated 
    # volcaniclastic) to a protolith. Make sure metasedimentary rocks are assigned 
    # to a metasedimentary protolith. Keep in mind that metased rocks may be 
    # currently assigned only to sed!!
    for i in eachindex(littletypes)
        if (littletypes[i] == :met) && !metased[i]
            littletypes[i] = protolith[weighted_rand(protolith_abundance)]

        elseif (littletypes[i] == :met) && metased[i]
            littletypes[i] = protolith_metased[weighted_rand(protolith_metased_abundance)]

        elseif (littletypes[i] == :sed) && metased[i]
            littletypes[i] = protolith_metased[weighted_rand(protolith_metased_abundance)]

        elseif littletypes[i] == :volcaniclast 
            littletypes[i] = minorvolc[weighted_rand(type_abundance.volc)]

        end
    end

    # Pass three: reassign nondescriptive types
    for i in eachindex(littletypes)
        if littletypes[i] == :sed
            littletypes[i] = minorsed[weighted_rand(type_abundance.sed)]

        elseif littletypes[i] == :volc 
            littletypes[i] = minorvolc[weighted_rand(type_abundance.volc)]

        elseif littletypes[i] == :plut 
            littletypes[i] = minorplut[weighted_rand(type_abundance.plut)]

        elseif littletypes[i] == :ign 
            # Pick a sub-class (volcanic / plutonic / carbonatite) and re-assign volc / plut
            littletypes[i] = minorign[weighted_rand(type_abundance.ign)]

            if littletypes[i] == :volc
                littletypes[i] = minorvolc[weighted_rand(type_abundance.volc)]

            elseif littletypes[i] == :plut 
                littletypes[i] = minorplut[weighted_rand(type_abundance.plut)]
                
            end

        end
    end


## --- Initialize for geochemical sample matching
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


## --- Find matching geochemical sample for each Macrostrat sample
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