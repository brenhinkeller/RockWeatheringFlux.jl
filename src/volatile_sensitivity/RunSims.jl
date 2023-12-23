# Test the sensitivity of upper crustal composition estimates to the amount of assumed
# volatiles in sedimentary rocks.

# Each simulation will further restrict the wt.% of volatiltes we assume are not 
# reported in sedimentary rocks. For example, the first run will allow samples to pass 
# the filter if they have less than 90 wt.% assumed volatiles (i.e., they report at 
# minimum 10 wt.%). 

# We also test the sensitivity of the estimtate to normalization. Setting a maximum
# allowed wt.% assumed volatiles of 16% means that we only assume volatiles if the 
# sample would already be allowed through the filter. In this case, the effect is that
# the reported composition is not normalized to 100%. We compare this to a simulation
# where we do not add any assumed volatiles.

# nohup julia test/volatile_sensitivity/VolatileSensitivityTesting.jl &

## --- Set up 
    # Packages
    using MAT
    using StatGeochem
    using LoopVectorization
    using HDF5
    using Static
    using Measurements
    using ProgressMeter
    using StatsBase
    using DelimitedFiles
    using StaticArrays
    using LogExpFunctions
    using Dates

    # Local utilites
    include("../../src/utilities/Utilities.jl")

    # Start timer
    start = now()
    @info """ Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))
    Input: $macrostrat_io
    Output: $matchedbulk_io
    """
    

## --- Load the Macrostrat / burwell source file
    fid = h5open("output/250K_responses.h5", "r")
    
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

    
## --- Calculate the total reported weight percent and assumed volatiles for each sample
    include("sim_ScreenBulk.jl")

    # How many samples did we start with?
    t_init = @. 84 <= bulkweight <= 104;
    init = count(t_init)
    t_init = vec(t_init)

    # Normalization etc. is always the same, we only vary *which* samples we keep
    # This shouldn't matter for the bulkweight calculation since that isn't recalculated
    simvolatiles = volatiles .+ additional
    normbulk = merge(bulk, (Volatiles=simvolatiles,))
    normbulk = NamedTuple{Tuple(allkeys)}([normbulk[i] for i in allkeys])

    contributing = [allelements; :Volatiles]                # Re-include volatiles!
    for i in eachindex(normbulk.SiO2)
        sample = [normbulk[j][i] for j in contributing]     # Get it
        normalize!(sample)                                  # Normalize it
        for j in eachindex(contributing)                    # Put it back
            normbulk[contributing[j]][i] = sample[j]
        end
    end
    

## --- Open a file to save data
    fid = h5open("test/volatile_sensitivity/simout.h5", "w")
    g = create_group(fid, "vars")

    # majors, minors = get_elements()       # This was done in sim_ScreenBulk.jl
    # allelements = [majors; minors]

    # Save the initial conditions
    g_init = create_group(g, "init")
    g_init["elements"] = string.([majors; minors])  # Order of crustal composition
    g_init["bulkweight"] = bulkweight               # Computed wt.%
    g_init["volatiles_known"] = volatiles           # Reported volatiles
    g_init["additional"] = additional               # Assumed sedimentary volatiles
    g_init["nsamples"] = init                       # Number of samples after initial filter

    # Create a group for the simulations
    g = create_group(fid, "sims")


## --- Pre-parse some of the data that should be invariant between simulations
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    # I hate cover
    macro_cats = delete_cover(macro_cats)
    bulk_cats = delete_cover(bulk_cats)

    # All granodiorites will also match with diorite. Don't do that.
    macro_cats.diorite .&= .!macro_cats.granodiorite
    bulk_cats.diorite .&= .!bulk_cats.granodiorite

    # Metamorphic rocks are only metamorphic if we cannot infer a protolith
    rox = keys(macro_cats)
    for type in rox
        type==:met && continue
        macro_cats.met .&= .!macro_cats[type]
        bulk_cats.met .&= .!bulk_cats[type]
    end

    # Calculate relative abundance of each rock type
    for type in minorvolc
        macro_cats.volc .|= macro_cats[type]
    end
    for type in minorplut
        macro_cats.plut .|= macro_cats[type]
    end

    # Absolute abundance of each rock type (count)
    ptype = (
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

    # Un-include all minor types in major types
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


## --- Pick meaningful major / minor types for all samples
    # Law of alrge numbers and all that, but I want to be sure that I'm only changing one thing.
    # Preallocate
    bigtypes = Array{Symbol}(undef, length(macrostrat.age), 1)      # Sed, volc, ign, etc.
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
        t = littletypes[i]
        if t == :sed
            littletypes[i] = minorsed[weighted_rand(ptype.sed)]
            bigtypes[i] = :sed

        elseif t == :ign 
            bigtypes[i] = minorign[weighted_rand(ptype.ign)]
            if bigtypes[i] == :volc
                littletypes[i] = minorvolc[weighted_rand(ptype.volc)]

            elseif bigtypes[i] == :plut 
                littletypes[i] = minorplut[weighted_rand(ptype.plut)]

            elseif bigtypes[i] == :carbonatite
                littletypes[i] = :carbonatite

            end

        elseif littletypes[i] != :none
            bigtypes[i] = class_up(littletypes[i], minorsed, minorvolc, minorplut, minorign)

        else
            bigtypes[i] = :none
        
        end
    end


## --- Now re-include all the types
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


## --- Run simulations with assumed volatiles
    geochemkeys = get_elements()[1][1:end-1]        # Major non-volatile elements

    simvaluesin = [90, 80, 70, 60, 50, 40, 30, 20, 16]

    for j in eachindex(simvaluesin)
        # Initialize
        @info "Starting sim for volatiles at $(simvaluesin[j]) wt.% at $(Dates.format(now(), "HH:MM"))"

        sim_fid = "sim_" * string(round(Int, simvaluesin[j]))
        sim_g = create_group(g, "$sim_fid")

        # Restrict data to 84-104 wt.% and cap assumed volatiles at limit to be tested
        t = @. 84 <= bulkweight .+ additional <= 104;   # Keep if adding volatiles brings into range
        t .&= additional .<= simvaluesin[j]             # UNLESS the volatiltes are above the cutoff
        sim_g["added"] = count(t) - init                # How many samples did we add?

        # Restrict data to this set of samples
        t = vec(t)
        simbulk = NamedTuple{keys(normbulk)}(normbulk[k][t] for k in keys(normbulk))
        bulk_kittens = NamedTuple{keys(bulk_cats)}(bulk_cats[k][t] for k in keys(bulk_cats))

        # Prepare to match samples
        bulk_inds = collect(1:length(simbulk.SiO2))     # Indices of bulk samples
        bulkzero = NamedTuple{Tuple(geochemkeys)}([zeronan(simbulk[i]) for i in geochemkeys])
        geochem_lookup = NamedTuple{Tuple(rox)}([major_elements(simbulk, bulk_kittens[i])
            for i in eachindex(rox)]
        );

        # Match samples
        matches = zeros(Int64, length(macro_cats.sed))
        for i in eachindex(matches)
            ltype = littletypes[i] 
            btype = bigtypes[i]
            if ltype == :none
                continue
            end
    
            # Assume the geochemical composition of the lithological sample: pick randomly
            randsample = rand(bulk_inds[bulk_kittens[ltype]])
            uncert = nanunzero!([geochem_lookup[ltype][elem].e for elem in geochemkeys], 1.0)
            values = [bulkzero[elem][randsample] for elem in geochemkeys]
    
            geochemdata = NamedTuple{Tuple(geochemkeys)}(
                NamedTuple{(:m, :e)}((values[j], uncert[j])) for j in eachindex(values)
            )
    
            # Get EarthChem data
            bulksamples = bulk_kittens[ltype]
            EC = (
                bulklat = simbulk.Latitude[bulksamples],            # EarthChem latitudes
                bulklon = simbulk.Longitude[bulksamples],           # EarthChem longitudes
                bulkage = simbulk.Age[bulksamples],                 # EarthChem age
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
        end
        sim_g["matches"] = matches

        # Restrict the bulk dataset down to just the matched samples
        t = @. matches > 0
        matchbulk = NamedTuple{Tuple(allkeys)}(
            [zeronan(simbulk[k][matches[t]]) for k in eachindex(allkeys)])

        # Calculate the bulk composition of continental crust
        # majors, minors = get_elements()       # This was done in sim_ScreenBulk.jl
        # allelements = [majors; minors]
        simUCC = [nanmean(matchbulk[i]) for i in allelements]
        sim_g["UCC"] = simUCC

        # Standard error
        n = length(matchbulk.SiO2)
        simSTD = [nanstd(matchbulk[i]) for i in allelements]
        sim_g["SEM"] = simSTD ./ sqrt(n)
    end


## --- Run a simulation assuming no volatiles
    sim_g = create_group(g, "sim_0")
    sim_g["added"] = 0

    # Restrict the dataset to filtered samples
    simbulk = merge(bulk, (Volatiles=volatiles,))
    simbulk = NamedTuple{Tuple(allkeys)}([simbulk[i][t_init] for i in allkeys])

    # Normalize compositions to 100%
    contributing = [allelements; :Volatiles]                # Re-include volatiles!
    for i in eachindex(simbulk.SiO2)
        sample = [simbulk[j][i] for j in contributing]      # Get it
        normalize!(sample)                                  # Normalize it
        for j in eachindex(contributing)                    # Put it back
            simbulk[contributing[j]][i] = sample[j]
        end
    end

    # Get rock types for this set of samples
    bulk_kittens = NamedTuple{keys(bulk_cats)}(bulk_cats[k][t_init] for k in keys(bulk_cats))

    # Prepare to match samples
    bulk_inds = collect(1:length(simbulk.SiO2))     # Indices of bulk samples
    bulkzero = NamedTuple{Tuple(geochemkeys)}([zeronan(simbulk[i]) for i in geochemkeys])
    geochem_lookup = NamedTuple{Tuple(rox)}([major_elements(simbulk, bulk_kittens[i])
        for i in eachindex(rox)]
    );
    
    # Match samples
    matches = zeros(Int64, length(macro_cats.sed))
    for i in eachindex(matches)
        ltype = littletypes[i] 
        btype = bigtypes[i]
        if ltype == :none
            continue
        end

        # Assume the geochemical composition of the lithological sample: pick randomly
        randsample = rand(bulk_inds[bulk_kittens[ltype]])
        uncert = nanunzero!([geochem_lookup[ltype][elem].e for elem in geochemkeys], 1.0)
        values = [bulkzero[elem][randsample] for elem in geochemkeys]

        geochemdata = NamedTuple{Tuple(geochemkeys)}(
            NamedTuple{(:m, :e)}((values[j], uncert[j])) for j in eachindex(values)
        )

        # Get EarthChem data
        bulksamples = bulk_kittens[btype]
        EC = (
            bulklat = simbulk.Latitude[bulksamples],            # EarthChem latitudes
            bulklon = simbulk.Longitude[bulksamples],           # EarthChem longitudes
            bulkage = simbulk.Age[bulksamples],                 # EarthChem age
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
    end
    sim_g["matches"] = matches

    # Restrict the bulk dataset down to just the matched samples
    t = @. matches > 0
    matchbulk = NamedTuple{Tuple(allkeys)}(
        [zeronan!(simbulk[k][matches[t]]) for k in eachindex(allkeys)])

    # Calculate the bulk composition of continental crust
    # majors, minors = get_elements()       # This was done in sim_ScreenBulk.jl
    # allelements = [majors; minors]
    simUCC = [nanmean(matchbulk[i]) for i in allelements]
    sim_g["UCC"] = simUCC

    # Standard error
    n = length(matchbulk.SiO2)
    simSTD = [nanstd(matchbulk[i]) for i in allelements]
    sim_g["SEM"] = simSTD ./ sqrt(n)


## --- End simulation: close file, stop timer
    close(fid)

    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """


## --- End of file