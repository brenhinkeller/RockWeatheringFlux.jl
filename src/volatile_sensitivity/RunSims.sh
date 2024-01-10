#! /usr/bin/bash

# Takes roughly an hour and a half to run
# nohup bash src/volatile_sensitivity/RunSims.sh &

# Fail if any command fails 
set -e 

# Start timer
printf "Starting simulations $(date +'%D %T')\n"


## Initialize variables
    # File names 
    simout="src/volatile_sensitivity/simout.h5"             # Simulation output file for results
    simout_prop="src/volatile_sensitivity/simout_prop.h5"   # Proportional volatiles output
    simout_prop2="src/volatile_sensitivity/simout_prop.h5"  # Proportional volatiles output version 2
    stem="src/volatile_sensitivity/simbulk_"                # Temporary simulation data file

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


## Run a set of simulations with proportional general / evaporite 
    # Bassanite / Dolomite = 0.91

    # Set up output file
    julia --project="../RockWeatheringFlux.jl/Project.toml" src/volatile_sensitivity/DefineOutput.jl $simout_prop

    # Define volatile measurements
    evap=( 90 80 70 60 50 40 30 20 16 )
    seds=( 81.9 72.8 63.7 54.6 45.5 36.4 27.3 18.2 14.56)

    for i in "${!evap[@]}"
    do
        printf "\nRunning simulation with ${evap[i]} wt.%% evaporite volatiles.\n"

        julia --project="../RockWeatheringFlux.jl/Project.toml" src/volatile_sensitivity/Simulation.jl\
        $simout_prop $stem ${seds[i]} ${seds[i]} ${evap[i]} ${evap[i]}
    done


## Run a set of simulations with proportional general sed / general evaporite / gypsum
    # Gypsum / Bassanite = 0.77
    # Gypsum / Dolomite = 0.71

    gyp_sim=( 90 80 70 60 50 40 30 20 16 )
    bas_sim=( 69.3 61.6 53.9 46.2 38.5 30.8 23.1 15.4 12.32 )
    dol_sim=( 63.9 56.8 49.7 42.6 35.5 28.4 21.3 14.2 11.36 )

    # Initialize output file
    julia --project="../RockWeatheringFlux.jl/Project.toml" src/volatile_sensitivity/DefineOutput.jl $simout_prop2

    # Run simulations
    for i in "${!gyp_sim[@]}"
    do
        printf "\nRunning simulation with ${gyp_sim[i]} wt.%% evaporite (gypsum) volatiles.\n"

        julia --project="../RockWeatheringFlux.jl/Project.toml" src/volatile_sensitivity/Simulation.jl\
        $simout_prop2 $stem ${dol_sim[i]} ${gyp_sim[i]} ${dol_sim[i]} ${bas_sim[i]}
    done

# Stop timer
printf "Ending simulations $(date +'%D %T')\n"

# TO DO:
# We also test the sensitivity of the estimtate to normalization. Setting a maximum
# allowed wt.% assumed volatiles of 16% means that we only assume volatiles if the 
# sample would already be allowed through the filter. In this case, the effect is that
# the reported composition is not normalized to 100%. We compare this to a simulation
# where we do not add any assumed volatiles.