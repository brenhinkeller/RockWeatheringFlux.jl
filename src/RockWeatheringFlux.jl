module RockWeatheringFlux

    using Reexport
    @reexport using StatGeochem
    @reexport using Measurements
    @reexport using ProgressMeter: @showprogress, Progress, next!
    
    # Vectorization tools
    using LoopVectorization: @turbo

    # File reading
    using HDF5, DelimitedFiles

    # Macrostrat API query 
    using JLD, HTTP, JSON

    # General requirements
    using Static 
    using Colors: RGB, Colorant
    using StatsBase: percentile, countmap
    using LogExpFunctions: logsumexp
    
    # Utilities
    include("../src/utilities/Definitions.jl")

    include("../src/utilities/Slope.jl")
    include("../src/utilities/NaNMeasurements.jl")
    include("../src/utilities/Macrostrat.jl")
    include("../src/utilities/Analysis.jl")
    include("../src/utilities/GTS.jl")
    include("../src/utilities/ScreenOutliers.jl")
    include("../src/utilities/Matching.jl")

    include("../src/utilities/Utilities.jl")

end # module