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
    init_macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    close(fid)
    

## --- Calculate the total reported weight percent and assumed volatiles for each sample
    include("sim_ScreenBulk.jl")

    # How many samples did we start with?
    t_init = @. 84 <= bulkweight <= 104;
    init = count(t)

    # Open a file to save data
    fid = h5open("simout.h5", "w")
    g = create_group(fid, "vars")

    # Save the initial conditions
    g_init = create_group(g, "init")
    g_init["elements"] = string.([majors; minors])  # Order of crustal composition
    g_init["bulkweight"] = bulkweight               # Computed wt.%
    g_init["volatiles_known"] = volatiles           # Reported volatiles
    g_init["additional"] = additional               # Assumed sedimentary volatiles
    g_init["nsamples"] = init                       # Number of samples after initial filter


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
    bulk_cats = match_rocktype(
        lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name][t]),
        lowercase.(bulktext.elements.Type[bulktext.index.Type][t]), 
        lowercase.(bulktext.elements.Material[bulktext.index.Material][t]), 
        source=:earthchem
    )

    # Make a copy of the Macrostrat rock type matches to avoid weird things happening
    macro_cats = deepcopy(init_macro_cats)

    # Match samples
    include("sim_SampleMatch.jl")
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


## --- Run simulations with assumed volatiles
    simvaluesin = [90, 80, 70, 60, 50, 40, 30, 20, 16]

    for j in eachindex(simvaluesin)
        # Create a file group for this simulation
        sim_fid = "sim_" * string(round(Int, simvaluesin[j]))
        sim_g = create_group(g, "$sim_fid")

        # Restrict data to 84-104 wt.% and cap assumed volatiles at limit to be tested
        t = @. 84 <= bulkweight .+ additional <= 104;   # Keep if adding volatiles brings into range
        t .&= additional .<= simvaluesin[j]             # UNLESS the volatiltes are above the cutoff

        simadded = count(t) - init                # How many samples did we add?
        sim_g["added"] = simadded

        # Restrict the dataset to filtered samples
        simvolatiles = volatiles .+ additional
        simbulk = merge(bulk, (Volatiles=simvolatiles,))
        simbulk = NamedTuple{Tuple(allkeys)}([simbulk[i][t] for i in allkeys])

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
        bulk_cats = match_rocktype(
            lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name][t]),
            lowercase.(bulktext.elements.Type[bulktext.index.Type][t]), 
            lowercase.(bulktext.elements.Material[bulktext.index.Material][t]), 
            source=:earthchem
        )

        # Make a copy of the Macrostrat rock type matches to avoid weird things happening
        macro_cats = deepcopy(init_macro_cats)

        # Match samples
        include("sim_SampleMatch.jl")
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
    end

    # End simulation: close file, stop timer
    close(fid)

    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """


## --- End of file