# Check that required data and directories are present, and create them if not

## --- Data files
    path = homedir()

    # Etopo 
    if !isfile("data/etopo/etopo1.h5")
        get_etopo()

        @info "Moving ETOPO1 data to $path/resources/ to $path/RockWeatheringFlux.jl/data/"
        run(`mv ../resources/etopo ../RockWeatheringFlux.jl/data`)
        run(`rm -rf mv ../resources/etopo`)
    end

    # SRTM15+ 
    if !isfile("data/srtm15plus.h5")
        get_srtm15plus()

        @info "Moving SRTM15+ data to $path/resources/ to $path/RockWeatheringFlux.jl/data/"
        run(`mv ../resources/srtm15plus ../RockWeatheringFlux.jl/data`)
        run(`rm -rf mv ../resources/srtm15plus`)
    end

    
## --- End of File