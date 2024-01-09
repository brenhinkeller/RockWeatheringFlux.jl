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
    c = 0

    # Find bulk silica distribution for iterative experiments
    for i in eachindex(simid)
        silica[i] = read(fid["sims"]["UCC"][simid[i]])[sᵢ]
        stdev[i] = read(fid["sims"]["UCC_stdev"][simid[i]])[sᵢ]

        try 
            sim[i] = parse(Float64, split(simid[i],"_")[2]) 
        catch
            c = i
            sim[i] = 47
        end
    end
    close(fid)

    # Build plot 
    h = Plots.plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica,
        ylabel="Average Crustal Silica [wt.%]",
        xlabel="Maximum Assumed Sedimentary Volatiles [wt.%]",
        legend=:bottomleft
    )
    Plots.plot!(sim[1:end .!= c], silica[1:end .!= c], ribbon=stdev[1:end .!= c], 
        # seriestype=:scatter,
        markershape=:circle,
        msc=:auto,
        label=""
    )
    Plots.plot!([sim[c]], [silica[c]], yerror=stdev[c], 
        markershape=:diamond,
        markersize=7,
        msc=:auto,
        label="Current Experimental Setup"
    )
    display(h)


## --- End of file 