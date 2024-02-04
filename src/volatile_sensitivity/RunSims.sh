#! /usr/bin/bash

# Takes roughly two and a half hours to run
# nohup bash src/volatile_sensitivity/RunSims.sh &

# Fail if any command fails 
set -e 

## Start timer
    printf "Starting simulations $(date +'%D %T')\n"


## Initialize variables
    # File names (remember to change screening file in Simulation.jl too!!)
    simout='src/volatile_sensitivity/simout_gard.h5'             # Simulation output file for results
    simout1='src/volatile_sensitivity/simout_prop_gard.h5'       # Proportional volatiles output
    simout2='src/volatile_sensitivity/simout_prop2_gard.h5'      # Proportional volatiles output version 2
    stem='src/volatile_sensitivity/simbulk_'                # Temporary simulation data file

    # Project location
    proj='../RockWeatheringFlux.jl/Project.toml'

    # Set arrays for volatile cutoffs
    # Dolomite = 1; Gypsum / Dolomite = 1.42; Bassanite / Dolomite = 1.29
    dol_sim=( 90 80 70 60 50 40 30 20 16 )
    gyp_sim=( 127.8  113.6  99.4  85.2  71.0  56.8  42.6  28.4  22.72 )
    bas_sim=( 116.1  103.2  90.3  77.4  64.5  51.6  38.7  25.8  20.64 )

    # Current experimental setup must be declared separately
    dol=47.5863523076   # Dolomite   (12.01+2*16)/((24.869+40.08)/2+12.01+16*3)*100 
    gyp=67.4237583502   # Gypsum     (32.07+16*3+2*(18))/(40.08+32.07+16*4+2*(18))*100
    bas=61.3641060971   # Bassanite  (32.07+16*3+0.5*(18))/(40.08+32.07+16*4+0.5*(18))*100


## Set up output files
    julia --project=$proj src/volatile_sensitivity/DefineOutput.jl $simout
    julia --project=$proj src/volatile_sensitivity/DefineOutput.jl $simout1
    julia --project=$proj src/volatile_sensitivity/DefineOutput.jl $simout2


## Run simulations
    # Pass arguments in the order:
        # 1) Simulation output file name
        # 2) Temporary simulation data file name "stem"
        # 3) File name extension for simulation (dolomite value)
        # 4) Volatile restrictions for gypsum, dolomite, and bassanite

    for i in "${!dol_sim[@]}"
    do
        printf "\nRunning simulation with ${dol_sim[i]} wt.%% volatiles $(date +'%D %T').\n"

        # Equal volatiles
        julia --project=$proj src/volatile_sensitivity/Simulation.jl\
        $simout $stem ${dol_sim[i]} ${dol_sim[i]} ${dol_sim[i]} ${dol_sim[i]}

        # Proportional volatiles (undifferentiated evaporites)
        julia --project=$proj src/volatile_sensitivity/Simulation.jl\
        $simout1 $stem ${dol_sim[i]} ${bas_sim[i]} ${dol_sim[i]} ${bas_sim[i]}

        # Proportional volatiles (differentiated evaporites)
        julia --project=$proj src/volatile_sensitivity/Simulation.jl\
        $simout2 $stem ${dol_sim[i]} ${gyp_sim[i]} ${dol_sim[i]} ${bas_sim[i]}
    done


## Run current experimental conditions
    printf "\nRunning simulation for current experimental conditions $(date +'%D %T').\n"

    julia --project=$proj src/volatile_sensitivity/Simulation.jl\
    $simout $stem 'initial' $gyp $dol $bas


## Stop timer
    printf "\nEnding simulations $(date +'%D %T')\n"


## --- End of file