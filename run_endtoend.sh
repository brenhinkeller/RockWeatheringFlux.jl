#! /bin/bash

# Run all code, end to end.

# Fail if any command fails
set -e

## --- Ready, set... 
    # Get lithologic data from Macrostrat API
    julia --project="Project.toml" src/ParseMacrostrat.jl

    # Filter geochemical data
    # julia --project="Project.toml" src/ScreenCombined.jl

    # Calculate slope from SRTM15+
    # julia --project="Project.toml" src/SRTMSlope.jl

    # Construct slope / erosion model (note model is hardcoded in src/utilities/Slope.jl)
    # julia --project="Project.toml" src/ModelSlope.jl


## --- ... go!
    julia --project="Project.toml" src/SampleMatch.jl       

    # Eroded material and mapped surfical abundance 
    julia --project="Project.toml" src/CalculateFlux.jl

    # Composition of exposed continental crust
    julia --project="Project.toml" src/UpperCrustComposition.jl     


## --- End of File 