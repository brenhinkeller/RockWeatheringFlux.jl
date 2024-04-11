#! /usr/bin/bash

# Takes roughly 2.5 hours to run for EarthChem, and 7-ish for Gard
# nohup bash src/volatile_sensitivity/RunSims.sh &

# Fail if any command fails 
set -e 

## Start timer
    printf "Starting simulations $(date +'%D %T')\n"


## Output file names -- remember to change the ScreenBase script in Simulation.jl
    # Note that if we want to move faster, we only actually plot simout2, although we do 
    # use simout to store the current experimental conditions

    # # EarthChem
    # simout='src/volatile_sensitivity/simout.h5'             # Simulation output file for results
    # simout1='src/volatile_sensitivity/simout_prop.h5'       # Proportional volatiles output
    # simout2='src/volatile_sensitivity/simout_prop2.h5'      # Proportional volatiles output version 2
    
    # # Gard 
    # simout='src/volatile_sensitivity/simout_gard.h5'
    # simout1='src/volatile_sensitivity/simout_prop_gard.h5'
    # simout2='src/volatile_sensitivity/simout_prop2_gard.h5'

    # Combined
    simout='src/volatile_sensitivity/simout_combo.h5'
    simout1='src/volatile_sensitivity/simout_prop_combo.h5'
    simout2='src/volatile_sensitivity/simout_prop2_combo.h5'


## Initialize other variables
    # Temporary simulation data file 
    stem='src/volatile_sensitivity/simbulk_'  

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
    # julia --project=$proj src/volatile_sensitivity/DefineOutput.jl $simout
    # julia --project=$proj src/volatile_sensitivity/DefineOutput.jl $simout1
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

        # # Equal volatiles
        # julia --project=$proj src/volatile_sensitivity/Simulation.jl\
        # $simout $stem ${dol_sim[i]} ${dol_sim[i]} ${dol_sim[i]} ${dol_sim[i]}

        # # Proportional volatiles (undifferentiated evaporites)
        # julia --project=$proj src/volatile_sensitivity/Simulation.jl\
        # $simout1 $stem ${dol_sim[i]} ${bas_sim[i]} ${dol_sim[i]} ${bas_sim[i]}

        # Proportional volatiles (differentiated evaporites)
        julia --project=$proj src/volatile_sensitivity/Simulation.jl\
        $simout2 $stem ${dol_sim[i]} ${gyp_sim[i]} ${dol_sim[i]} ${bas_sim[i]}
    done


## Run current experimental conditions
    printf "\nRunning simulation for current experimental conditions $(date +'%D %T').\n"

    julia --project=$proj src/volatile_sensitivity/Simulation.jl\
    $simout2 $stem 'initial' $gyp $dol $bas


## Stop timer
    printf "\nEnding simulations $(date +'%D %T')\n"


## --- End of file