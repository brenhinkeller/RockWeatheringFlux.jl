# Wrapper for SampleMatchBase.jl, whihc matches Macrostrat samples to the most likely 
# EarthChem sample

## --- Code
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, StaticArrays, Dates

    # File name IO
    filemacrostrat = macrostrat_io
    filebulk = "output/bulk.h5"

    # Start timer
    start = now()
    @info """ Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))
    Input: $macrostrat_io
    Output: $matchedbulk_io
    """

    # Call file and write output to a file
    include("SampleMatchBase.jl")
    writedlm("$matchedbulk_io", [matches string.(littletypes)], '\t')

    # End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Program runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """
    

## --- End of File