#! /bin/bash

# Run all code, end to end.

# Fail if any command fails
set -e

## --- Generate data files. Optionally skip this step and use existing data files
    # Get lithologic map data for N points on the continental crust
    # Takes ~4 days for the 1M sample request. Use [output/N_*/*_responses.h5]
    julia --project="Project.toml" src/ParseMacrostrat.jl

    # Dataset used is controlled by a switch in src/utilities/Definitions.jl
    julia --project="Project.toml" src/ScreenBulk.jl        # [output/geochemistry/bulk.h5]
    julia --project="Project.toml" src/ScreenGard.jl        # [output/geochemistry/gard.h5]
    julia --project="Project.toml" src/ScreenCombined.jl    # [output/geochemistry/combined.h5]

    # Compute slope at each point on Earth. 
    # Takes a couple hours. Use [output/srtm15plus_maxslope.h5]
    julia --project="Project.toml" src/SRTMSlope.jl


## --- Construct slope / erosion model
    # This is hardcoded in src/utilities/Slope.jl
    julia --project="Project.toml" src/ModelSlope.jl

## --- Match lithology and geochemistry, estimate composition of crust and eroded material
    # [output/N_*/*_bulkidx_[gard/bulk].tsv]
    julia --project="Project.toml" src/SampleMatch.jl       

    # Absolute (Gt) composition [results/*_eroded_absolute_[gard/bulk].tsv]
    # Fractional contribution by class  [results/*_eroded_fraction_[gard/bulk].tsv]
    # Composition (wt.%) [results/*_eroded_composition_[gard/bulk].tsv]
    julia --project="Project.toml" src/CalculateFlux.jl

    # [results/*_exposedcrust_[gard/bulk].tsv]
    julia --project="Project.toml" src/UpperCrustComposition.jl     


## --- End of File 