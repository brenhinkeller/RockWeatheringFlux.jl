# Look at the simulated crust compositions and see what happens

## --- Set up
    # Packages
    using RockWeatheringFlux
    using HDF5, Plots


## --- Load simulation with equal assumed volatiles and current experimental conditions
    # Load data and find SiO₂ index
    fid = h5open("src/volatile_sensitivity/simout.h5", "r")
    allelements = read(fid["vars"]["elements"])
    sᵢ = findfirst(x->x=="SiO2", allelements)

    # Preallocate 
    simid = keys(fid["sims"]["UCC"])
    sim = Array{Float64}(undef, length(simid), 1)
    silica = Array{Float64}(undef, length(simid), 1)
    stdev = Array{Float64}(undef, length(simid), 1)
    c = 0

    # Find bulk silica distribution for iterative experiments
    for i in eachindex(simid)
        silica[i] = read(fid["sims"]["UCC"][simid[i]])[sᵢ]
        stdev[i] = read(fid["sims"]["UCC_stdev"][simid[i]])[sᵢ]

        try 
            sim[i] = parse(Float64, split(simid[i],"_")[2]) 
        catch
            # Catch current experimental conditions
            c = i
            sim[i] = 47
        end
    end
    close(fid)


## --- Load simulation with proportional assumed volatiles (bassanite:dolomite)
    fid = h5open("src/volatile_sensitivity/simout_prop.h5", "r")
    allelements = read(fid["vars"]["elements"])
    sᵢ = findfirst(x->x=="SiO2", allelements)
    
    # Preallocate 
    simid = keys(fid["sims"]["UCC"])
    sim_prop = Array{Float64}(undef, length(simid), 1)
    silica_prop = Array{Float64}(undef, length(simid), 1)
    stdev_prop = Array{Float64}(undef, length(simid), 1)

    # Find bulk silica distribution for iterative experiments
    for i in eachindex(simid)
        sim_prop[i] = parse(Float64, split(simid[i],"_")[2]) 
        silica_prop[i] = read(fid["sims"]["UCC"][simid[i]])[sᵢ]
        stdev_prop[i] = read(fid["sims"]["UCC_stdev"][simid[i]])[sᵢ]
    end
    close(fid)


## --- Load simulation with proportional assumed volatiles (gypsum:dolomite:bassanite)
    fid = h5open("src/volatile_sensitivity/simout_prop.h5", "r")
    allelements = read(fid["vars"]["elements"])
    sᵢ = findfirst(x->x=="SiO2", allelements)

    # Preallocate 
    simid = keys(fid["sims"]["UCC"])
    sim_prop2 = Array{Float64}(undef, length(simid), 1)
    silica_prop2 = Array{Float64}(undef, length(simid), 1)
    stdev_prop2 = Array{Float64}(undef, length(simid), 1)

    # Find bulk silica distribution for iterative experiments
    for i in eachindex(simid)
        sim_prop2[i] = parse(Float64, split(simid[i],"_")[2]) 
        silica_prop2[i] = read(fid["sims"]["UCC"][simid[i]])[sᵢ]
        stdev_prop2[i] = read(fid["sims"]["UCC_stdev"][simid[i]])[sᵢ]
    end
    close(fid)


## --- Plot simulations
    # Build plot 
    h = Plots.plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica,
        ylabel="Average Crustal Silica [wt.%]",
        xlabel="Maximum Assumed Sedimentary Volatiles [wt.%]",
        legend=:bottomleft,
        fg_color_legend=:white,
    )

    # Equal volatiles
    Plots.plot!(sim[1:end .!= c], silica[1:end .!= c], ribbon=stdev[1:end .!= c], 
        # seriestype=:scatter,
        markershape=:circle,
        msc=:auto,
        fillalpha=0.5,
        label="Carbonate = Evaporite = Gypsum"
    )

    # Proportional volatiles (carbonate / evaporite)
    Plots.plot!(sim_prop, silica_prop, ribbon=stdev_prop, 
        markershape=:circle,
        msc=:auto,
        fillalpha=0.5,
        label="Carbonate ≠ Evaporite = Gypsum"
    )

    # Proportional volatiles (carbonate / evaporite / gypsum)
    Plots.plot!(sim_prop2, silica_prop2, ribbon=stdev_prop2, 
        markershape=:circle,
        msc=:auto,
        fillalpha=0.5,
        label="Carbonate ≠ Evaporite ≠ Gypsum"
    )

    # Current experimental conditions 
    Plots.plot!([sim[c]], [silica[c]], yerror=stdev[c], 
        markershape=:diamond,
        markersize=7,
        linewidth=2,
        color=:black, linecolor=:black, msc=:auto,
        label="Current Experimental Setup"
    )
    display(h)

## --- End of file 