#! /usr/bin/bash

# nohup bash src/volatile_sensitivity/RunSims.sh

# Fail if any command fails 
set -e 

# Start timer
printf "Starting simulations $(date +'%D %T')\n"


## Initialize variables
    # File names 
    simout="src/volatile_sensitivity/simout.h5"     # Simulation output file for results
    stem="src/volatile_sensitivity/simbulk_"        # Temporary simulation data file

    # Set array for volatile cutoffs
    volatiles=( 90 80 70 60 50 40 30 20 16 )

    # Current experimental setup must be declared separately
    gyp=47.5863523076   # Dolomite   (12.01+2*16)/((24.869+40.08)/2+12.01+16*3)*100 
    dol=67.4237583502   # Gypsum     (32.07+16*3+2*(18))/(40.08+32.07+16*4+2*(18))*100
    bas=61.3641060971   # Bassanite  (32.07+16*3+0.5*(18))/(40.08+32.07+16*4+0.5*(18))*100

## Set up output file
    julia --project="../RockWeatheringFlux.jl/Project.toml" src/volatile_sensitivity/DefineOutput.jl $simout

## Run simulations
    # Pass arguments in the order:
        # 1) Simulation output file name
        # 2) Temporary simulation data file name "stem"
        # 3) File name extension for simulation
        # 4) Volatile restrictions for gypsum, dolomite, and bassanite

    # Run simulations with volatiles 
    for i in "${volatiles[@]}"
    do
        printf "\nRunning simulation with ${i} wt.%% volatiles.\n"

        julia --project="../RockWeatheringFlux.jl/Project.toml" src/volatile_sensitivity/Simulation.jl\
        $simout $stem $i $i $i $i
    done

    # Run current experimental conditions
    printf "\nRunning simulation for current experimental conditions.\n"
    fname='initial'
    julia --project="../RockWeatheringFlux.jl/Project.toml" src/volatile_sensitivity/Simulation.jl\
    $simout $stem $fname $gyp $dol $bas


# Stop timer
printf "Ending simulations $(date +'%D %T')\n"

# TO DO:
# We also test the sensitivity of the estimtate to normalization. Setting a maximum
# allowed wt.% assumed volatiles of 16% means that we only assume volatiles if the 
# sample would already be allowed through the filter. In this case, the effect is that
# the reported composition is not normalized to 100%. We compare this to a simulation
# where we do not add any assumed volatiles.