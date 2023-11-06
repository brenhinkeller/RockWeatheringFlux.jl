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

    # Local utilites
    include("../../src/utilities/Utilities.jl")
    

    # Calculate the reported weight and additional sedimentary volatiles needed to get 
    # to 100% for each sample
    include("ReportedWeight.jl")

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

        # How many samples did we add?
        simadd = count(t) - init

        # Get rock types for the samples we're keeping
        bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name][t])
        bulkrocktype = lowercase.(bulktext.elements.Type[bulktext.index.Type][t])
        bulkmaterial = lowercase.(bulktext.elements.Material[bulktext.index.Material][t])

        bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial, source=:earthchem)

        # Match samples

        # Calculate the bulk composition of continental crust

        # Save data to a file
    end

## --- End of file