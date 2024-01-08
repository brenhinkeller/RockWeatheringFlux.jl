# Look at the simulated crust compositions and see what happens

## --- Set up
    # Packages
    using RockWeatheringFlux
    using HDF5, Plots

    # Load data and find SiO₂ index
    fid = h5open("src/volatile_sensitivity/simout.h5", "r")
    allelements = read(fid["vars"]["elements"])
    sᵢ = findfirst(x->x=="SiO2", allelements)


## --- Plot bulk silica vs. max. assumed volatiles 
    # Preallocate 
    simid = keys(fid["sims"]["UCC"])
    sim = Array{Float64}(undef, length(simid), 1)
    silica = Array{Float64}(undef, length(simid), 1)
    stdev = Array{Float64}(undef, length(simid), 1)

    # Find bulk silica distribution
    for i in eachindex(simid)
        silica[i] = read(fid["sims"]["UCC"][simid[i]])[sᵢ]
        stdev[i] = read(fid["sims"]["UCC_stdev"][simid[i]])[sᵢ]

        try 
            sim[i] = parse(Float64, simid[i]) 
        catch
            sim[i] = 47
        end
    end
    close(fid)

    # Build plot 
    h = Plots.plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica, 
        # xlims=(40,80),
        # xticks=(40:10:80, string.(40:10:80)),
        yticks=false,
        ylabel="Average Crustal Silica [wt.%]",
        xlabel="Maximum Assumed Sedimentary Volatiles [wt.%]",
    )
    Plots.plot!(simvalue, silica, yerror=stdev*2, 
        seriestype=:scatter,
        label=""
    )
    display(h)


## --- End of file 