## --- Run simulation for a single volatile restriction set passed to this script

## --- Set up
    # Packages
    using RockWeatheringFlux
    using MAT, HDF5, DelimitedFiles, StaticArrays, Dates
    using LoopVectorization: @turbo

    # List of elements 
    majors, minors = get_elements()
    allelements = [majors; minors]

    # Parse arguments
    simout = ARGS[1]                        # Simulation output file
    fname = ARGS[3]                         # File name extension
    fileout = ARGS[2] * fname * ".h5"       # Temporary simulation output file 
    gyp = parse(Float64, ARGS[4])           # Volatile restriction 
    dol = parse(Float64, ARGS[5])
    bas = parse(Float64, ARGS[6])

    # Suppress progress bars
    show_progress = false


## --- Run simulation
    # Open simulation output file
    fid_out = h5open(simout, "r+")

    # Run screening restrictions
    # include("../ScreenBulkBase.jl")
    # include("../ScreenGardBase.jl")
    include("../ScreenCombinedBase.jl")

    # Save BitVector and sample increase to simout file
    a = zeros(Int64, length(t))
    a[t] .= 1
    fid_out["sims"]["BitVectors"]["sim_$(fname)"] = a
    fid_out["sims"]["addition"]["sim_$(fname)"] = up

    # Run SampleMatchBase with restricted dataset
    filemacrostrat = macrostrat_io          # Macrostrat dataset
    filebulk = fileout                      # Temporary simulation output file 
    include("../SampleMatchBase.jl")

    # Write output to the simout file
    fid_out["sims"]["indices"]["sim_$(fname)"] = matches
    fid_out["sims"]["class_assigned"]["sim_$(fname)"] = string.(littletypes)

    # Compute composition of upper continental crust and save to simout file
    s = matches .!= 0;
    mbulk = NamedTuple{keys(bulk)}([bulk[k][matches[s]] for k in keys(bulk)])
    for k in allelements
        zeronan!(mbulk[k])
    end
    fid_out["sims"]["UCC"]["sim_$(fname)"] = [nanmean(mbulk[i]) for i in allelements]
    fid_out["sims"]["UCC_std"]["sim_$(fname)"] = [
        nanstd(mbulk[i]) for i in allelements
    ]

    # Close simulation output file and delete temporary simulation file
    close(fid_out)
    run(`rm $filebulk`)


## --- End of File 