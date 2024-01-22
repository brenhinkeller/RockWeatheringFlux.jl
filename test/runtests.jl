## --- Compile all code 
    using RockWeatheringFlux
    using Test, StatsBase
    # using Statistics


## --- Run tests    
    # To do: remove the functions I'm not using anymore. 
    # Organize debugging functions into a better-named debugging file 
    # (this is currently Analysis.jl)

    # @testset "Analysis                   "
    # @testset "Definitions                "
    @testset "GTS                        " begin include("testGTS.jl") end
    # @testset "Macrostrat                 "
    @testset "NaN Measurements           " begin include("testNaNMeasurements.jl") end
    @testset "Screen Geochemical Outliers" begin include("testScreenOutliers.jl") end
    # @testset "Slope                      "
    @testset "Other Utiliities           " begin include("testUtilities.jl") end


## --- End of file