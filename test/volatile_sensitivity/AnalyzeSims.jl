# Look at the simulated crust compositions and see what happens

## --- Set up
    using StatGeochem
    using HDF5
    using Plots


## --- Load data 
    fid = h5open("test/volatile_sensitivity/simout.h5", "r")

    # Find index of SiO2 measurement
    allelements = read(fid["vars"]["init"]["elements"])
    sᵢ = findfirst(x->x=="SiO2", allelements)

    # How many simulations did we run?
    folders = keys(fid["vars"]["sims"])
    nsims = length(folders)

    # Preallcoate
    silica = Array{Float64}(undef, nsims)
    sem = Array{Float64}(undef, nsims)
    simvalue = Array{Float64}(undef, nsims)
    
    # Get data
    for i in folders
        k = folders[i]

        # Get the silica composition of the crust
        comp = read(fid["vars"]["sims"][k]["UCC"])
        err = read(fid["vars"]["sims"][k]["SEM"])
        silica[i] = comp[sᵢ]
        sem[i] = err[sᵢ]

        # Get the maximum wt.% assumed volatiles 
        val = split(k, "_")[2]
        simvalue[i] = parse(Float64, val)
    end

    close(fid)


## --- Plot dependence on assumed volatiles
    # We are plotting silica as a function of assumed volatiles
    h = plot(simvalue, silica, seriestype=:scatter, framestyle=:box,
        ylabel="Silica [wt.%]", xlabel="Maximum Assumed Sedimentary Volatiles [wt.%]",
        label=""
    )
    display(h)


## --- I want to add SEM but I don't want to run everything again. Load matches
    # fid = h5open("test/volatile_sensitivity/simout.h5", "r")
    # folders = keys(fid["vars"])
    # nsims = 0
    # for k in folders
    #     occursin("sim", k) && (nsims += 1)
    # end

    # # Find index of SiO2 measurement
    # allelements = read(fid["vars"]["init"]["elements"])
    # sᵢ = findfirst(x->x=="SiO2", allelements)

    # # Preallcoate
    # silica = Array{Float64}(undef, nsims)
    # sem = Array{Float64}(undef, nsims)
    # simvalue = Array{Float64}(undef, nsims)

    # # Re-compute crustal composition and get SEM
    # i = 1   # Generic counter
    # for k in folders
    #     !occursin("sim", k) && continue

    #     # Load matches
    #     matches = read(fid["vars"]["k"]["matches"])
    #     t = @. matches > 0
    #     matchbulk = NamedTuple{Tuple(allelements)}(
    #         [zeronan!(simbulk[e][matches[t]]) for e in eachindex(allelements)])

    #     n = length(matchbulk.SiO2)
    #     simUCC = [nanmean(matchbulk[i]) for i in allelements]
    #     simSEM = [nanstd(matchbulk[i]) for i in allelements] ./ sqrt(n)

    #     # Get the silica composition of the crust
    #     silica[i] = simUCC[sᵢ]
    #     sem[i] = simSEM[sᵢ]

    #     # Get the maximum wt.% assumed volatiles 
    #     val = split(k, "_")[2]
    #     simvalue[i] = parse(Float64, val)
    #     i += 1
    # end

    # # And visualize
    # h = plot(simvalue, silica, ribbbon=:sem, seriestype=:scatter, framestyle=:box,
    #     ylabel="Silica [wt.%]", xlabel="Maximum Assumed Sedimentary Volatiles [wt.%]",
    #     label=""
    # )
    # display(h)


## --- End of file 