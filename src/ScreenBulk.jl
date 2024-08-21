# This file is a wrapper for ScreenBulkBase.jl, which reads bulk.mat and bulktext.mat,
# converts all units to wt.%, restricts data to whole-rock analyses, and saves the data as
# a new file.

# Using a wrapper allows me to call ScreenBulkBase.jl in sensitivity tests, meaning I can 
# use the same code to run my tests as I do in the actual computations. This should reduce
# unintended errors from editing or not updating the simulation code!

## --- Code
    # Packages
    using RockWeatheringFlux
    using MAT, HDF5, Dates
    using LoopVectorization: @turbo
    
    # File name IO
    fileout = "output/geochemistry/bulk.h5"

    # Start timer 
    start = now()
    @info """ Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))
    Output: $fileout
    """

    # Calculate a reasonable assumption for wt.% volatiles
    # Differentiating dolomite from evaporites doesn't actually matter: see sensitivity
    # testing results
    dol = (12.01+2*16)/((24.869+40.08)/2+12.01+16*3)*100          # Dolomite
    gyp = (32.07+16*3+2*(18))/(40.08+32.07+16*4+2*(18))*100       # Gypsum
    bas = (32.07+16*3+0.5*(18))/(40.08+32.07+16*4+0.5*(18))*100   # Bassanite (2CaSO₄⋅H₂O)

    # Call file
    show_progress = true
    include("ScreenBulkBase.jl")

    # End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Program runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """

## --- End of File