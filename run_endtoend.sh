#! /bin/bash

# Run all code, end to end.

# Fail if any command fails
set -e

## --- Generate data files.
    julia --project="Project.toml" src/ParseMacrostrat.jl   # Lithologic data
    julia --project="Project.toml" src/ScreenCombined.jl    # Filtered geochemical data
    julia --project="Project.toml" src/SRTMSlope.jl         # Slope


## --- Construct slope / erosion model (note model is hardcoded in src/utilities/Slope.jl)
    julia --project="Project.toml" src/ModelSlope.jl


## --- Match lithology and geochemistry
    julia --project="Project.toml" src/SampleMatch.jl       


## --- Results
    # Mass of eroded material: undifferentiated, by element, and by lithologic class
    # Fractional contribution of each lithologic class to total eroded material
    # Composition of eroded material by lithologic class
    # Surficial abundance vs. contribution to total eroded material by lithologic class
    # Mapped surficial abundance of each lithologic class
    julia --project="Project.toml" src/CalculateFlux.jl

    # Composition of exposed continental crust
    # Expected surficial distribution of major lithologic classes
    julia --project="Project.toml" src/UpperCrustComposition.jl     


## --- End of File 