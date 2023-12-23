## --- Compile all code 
    using RockWeatheringFlux
    using Test, StatsBase
    # using Statistics


## --- Run tests    
    # Need unit tests for:
        # Analysis 
        # Definitions 
        # Macrostrat 
        # Slope
    # To do: remove the functions I'm not using anymore. Organize debugging functions 
    # into a better-named debugging file (this is currently Analysis.jl)
    @testset "NaN Measurements" begin include("testNaNMeasurements.jl") end
    @testset "Other Utiliities" begin include("testUtilities.jl") end

## --- End of file