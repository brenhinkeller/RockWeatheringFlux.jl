module RockWeatheringFlux

    using Reexport
    @reexport using StatGeochem
    @reexport using Measurements
    @reexport using ProgressMeter: @showprogress, Progress, next!
    
    # Vectorization tools
    using LoopVectorization: @turbo

    # General requirements
    using Static 
    using Colors: RGB
    using StatsBase: percentile
    using LogExpFunctions: logsumexp
    
    # Utilities
    include("../src/utilities/Definitions.jl")

    include("../src/utilities/Slope.jl")
    include("../src/utilities/NaNMeasurements.jl")
    include("../src/utilities/Macrostrat.jl")
    include("../src/utilities/Analysis.jl")
    include("../src/utilities/GTS.jl")

    include("../src/utilities/Utilities.jl")

end # module