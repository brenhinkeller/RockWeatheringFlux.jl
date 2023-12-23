# To do: make this so I can just put using RWF in my files?
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
    
    # Utilities
    include("../src/utilities/Definitions.jl")

    include("../src/utilities/Slope.jl")
    include("../src/utilities/NaNMeasurements.jl")
    include("../src/utilities/Macrostrat.jl")
    include("../src/utilities/Analysis.jl")
    include("../src/utilities/Spidergram.jl")

    include("../src/utilities/Utilities.jl")

end # module