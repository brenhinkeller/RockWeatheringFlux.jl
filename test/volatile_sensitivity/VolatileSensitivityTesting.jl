# Test the sensitivity of upper crustal composition estimates to the amount of assumed
# volatiles in sedimentary rocks.

# Calculate the reported weight percent of all EarthChem samples. Iteratively match
# samples and calculate the bulk composition of the continental crust, progressively
# increasing the reported analyzed weight percent required to assume the rest of the rock
# is unreported / unanalyzed volatiles. Save the compositions to a file for comparison.
# Ideally, we want there to be a relatively large "buffer area," such that there is a 
# region where the composition of continental crust is not a function of the fraction of 
# assumed volatiles.

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

    # Local utilites
    include("../../src/utilities/Utilities.jl")
    
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
    t = @. 84 <= bulkweight <= 104
    init = count(t)


## --- Run simulations
    # Minimum reported weight percent required to assume the rest are volatiles. That is,
    # if only 5 wt.% is required, we will assign an additional 95% volatiles
    simvaluesin = [5, 10, 20, 30, 40, 50, 60, 70, 80, 84]

    for i in eachindex(simvaluesin)
        # Restrict data to 84-104 wt.% and cap assumed volatiles at limit to be tested
        t = @. 84 <= bulkweight .+ additional <= 104
        t .&= additional .< simvaluesin[i]
        simadd = count(t) - init                # How many samples did we add?

        # Get rock types for this set of samples
        bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name][t])
        bulkrocktype = lowercase.(bulktext.elements.Type[bulktext.index.Type][t])
        bulkmaterial = lowercase.(bulktext.elements.Material[bulktext.index.Material][t])
        bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial, source=:earthchem)

        # Restrict the dataset and normalize compositions to 100%
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

        # Make a copy of the Macrostrat rock type matches to avoid weird things happening
        macro_cats = deepcopy(init_macro_cats)

        # Match samples
        include("sim_SampleMatch.jl")

        # Calculate the bulk composition of continental crust

        # Save data to a file
    end

## --- End of file