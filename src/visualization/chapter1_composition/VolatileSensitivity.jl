## --- Set up 
    # Visualize simulated crust compositions as a function of assumed volatiles

    # Load data and base packages
    if !@isdefined(filepath)
        include("Definitions.jl")
    end

    # Define simulation filepaths
    sim1 = "src/volatile_sensitivity/simout.h5"
    sim2 = "src/volatile_sensitivity/simout_prop.h5"
    sim3 = "src/volatile_sensitivity/simout_prop2.h5"


## --- Load simulation with equal assumed volatiles and current experimental conditions
    # Load data and find SiO₂ index
    fid = h5open(sim1, "r")
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
        stdev[i] = read(fid["sims"]["UCC_sem"][simid[i]])[sᵢ]

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
    fid = h5open(sim2, "r")
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
        stdev_prop[i] = read(fid["sims"]["UCC_sem"][simid[i]])[sᵢ]
    end
    close(fid)


## --- Load simulation with proportional assumed volatiles (gypsum:dolomite:bassanite)
    fid = h5open(sim3, "r")
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
        stdev_prop2[i] = read(fid["sims"]["UCC_sem"][simid[i]])[sᵢ]
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
        xlims=(10,100),
        legend=:bottomleft,
        fg_color_legend=:white,
    )

    # Plot other estimates for comparison
    Plots.hline!([66.62], label="Rudnick and Gao, 2014", 
        linewidth=3, linestyle=:dash, color=colors.siliciclast
    )
    Plots.hline!([59.5], label="Pease et al., 2023 [EarthChem prior]", 
        linewidth=3, linestyle=:dash, color=colors.plut
    )

    # # Equal volatiles
    # Plots.plot!(sim[1:end .!= c], silica[1:end .!= c], ribbon=stdev[1:end .!= c], 
    #     markershape=:circle,
    #     msc=:auto,
    #     fillalpha=0.15,
    #     label="Carbonate = Evaporite = Gypsum",
    #     # label=""
    # )

    # # Proportional volatiles (carbonate / evaporite)
    # Plots.plot!(sim_prop, silica_prop, ribbon=stdev_prop, 
    #     markershape=:circle,
    #     msc=:auto,
    #     fillalpha=0.15,
    #     label="Carbonate ≠ Evaporite = Gypsum"
    # )

    # Proportional volatiles (carbonate / evaporite / gypsum)
    Plots.plot!(sim_prop2, silica_prop2, ribbon=stdev_prop2, 
        markershape=:circle,
        color=colors.evap, msc=:auto,
        fillalpha=0.15,
        # label="Carbonate ≠ Evaporite ≠ Gypsum"
        label="Simulations"
    )

    # Current experimental conditions 
    Plots.scatter!([sim[c]], [silica[c]], # yerror=stdev[c], 
        markershape=:star5,
        markersize=10,
        linewidth=1,
        color=colors.rhyolite, linecolor=:black, msc=:black,
        label="Current Experimental Conditions"
    )

    display(h)
    savefig(h, "$filepath/volatile_sensitivity.pdf")


## --- End of file