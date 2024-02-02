# Wrapper for SampleMatchBase.jl, which matches Macrostrat samples to the most likely 
# EarthChem sample

## --- Code
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, StaticArrays, Dates

    # File name IO
    filemacrostrat = macrostrat_io
    filebulk = geochem_fid

    # Start timer
    start = now()
    @info """ Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))
    Input: $macrostrat_io | $geochem_fid
    Output: $matchedbulk_io
    """

    # Call file and write output to a file
    show_progress = true
    include("SampleMatchBase.jl")
    writedlm("$matchedbulk_io", [matches string.(littletypes)], '\t')

    # End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Program runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """
    

## --- End of File