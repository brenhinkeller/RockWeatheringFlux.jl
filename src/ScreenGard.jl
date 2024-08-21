# This file is a wrapper for ScreenGardBase.jl, which reads data from Gard et al., 2019
# (10.5194/essd-11-1553-2019) converts all units to wt.%, restricts data to whole-rock 
# analyses, and saves the data as a new file.

## --- Code
    # Packages
    using RockWeatheringFlux
    using HDF5, Dates
    using LoopVectorization: @turbo

    # File name IO
    fileout = "output/geochemistry/gard.h5"

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
    include("ScreenGardBase.jl")

    # End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Program runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """
    

## --- End of file 