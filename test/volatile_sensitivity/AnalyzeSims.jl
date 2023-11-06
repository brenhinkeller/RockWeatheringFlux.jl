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
    folders = keys(fid["vars"])
    nsims = 0
    for k in folders
        occursin("sim", k) && (nsims += 1)
    end

    # Preallcoate
    silica = Array{Float64}(undef, nsims)
    simvalue = Array{Float64}(undef, nsims)
    
    # Get data
    i = 1   # Counter
    for k in folders
        !occursin("sim", k) && continue

        # Get the silica composition of the crust
        comp = read(fid["vars"][k]["UCC"])
        silica[i] = comp[sᵢ]

        # Get the maximum wt.% assumed volatiles 
        val = split(k, "_")[2]
        simvalue[i] = parse(Float64, val)
        i += 1
    end

    close(fid)


## --- Plot dependence on assumed volatiles
    # We are plotting silica as a function of assumed volatiles
    h = plot(simvalue, silica, seriestype=:scatter, framestyle=:box,
        ylabel="Silica [wt.%]", xlabel="Maximum Assumed Sedimentary Volatiles [wt.%]",
        label=""
    )
    display(h)

    
## --- End of file 