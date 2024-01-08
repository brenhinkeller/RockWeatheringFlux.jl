# Test the sensitivity of upper crustal composition estimates to the amount of assumed
# volatiles in sedimentary rocks.

# nohup julia src/volatile_sensitivity/RunSims.jl &
# TO DO: I need to activate an environment to run this, but I can't seem to do that via 
# terminal to run a nohup process :(

## --- Set up 
    # using Pkg
    # Pkg.activate()

    # Packages
    using RockWeatheringFlux
    using MAT, HDF5, DelimitedFiles, StaticArrays, Dates
    using LoopVectorization: @turbo

    # Definitions
    simout = "src/volatile_sensitivity/simout.h5"
    stem = "src/volatile_sensitivity/simbulk_"

    # Start timer
    start = now()
    @info """ Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))
    Output: $simout
    """


## --- Initialize test conditions 
    # Current experimental setup
    dol = (12.01+2*16)/((24.869+40.08)/2+12.01+16*3)*100          # Dolomite
    gyp = (32.07+16*3+2*(18))/(40.08+32.07+16*4+2*(18))*100       # Gypsum
    bas = (32.07+16*3+0.5*(18))/(40.08+32.07+16*4+0.5*(18))*100   # Bassanite (2CaSO₄⋅H₂O)

    # Test conditions: maximum allowable wt.% total volatiles
    init = (
        ("90", (90,90,90)),         # Any sample with > 10 wt.% reported
        ("80", (80,80,80)),
        ("70", (70,70,70)),
        ("60", (60,60,60)),
        ("50", (50,50,50)),
        ("init", (dol,gyp,bas)),    # Current experimental setup
        ("40", (40,40,40)),
        ("30", (30,30,30)),
        ("20", (20,20,20)),
        ("16", (16,16,16)),         # Samples that would have gotten through anyway
    )

    # List of elements 
    majors, minors = get_elements()
    allelements = [majors; minors]

    # Open file 
    fid = h5open(simout, "w")
    vars = create_group(fid, "vars")
    vars["elements"] = string.(allelements)
    
    sims = create_group(fid, "sims")
    g_bitvec = create_group(sims, "BitVectors")
    g_index = create_group(sims, "indices")
    g_types = create_group(sims, "type_assigned")
    g_crust = create_group(sims, "UCC")
    g_err = create_group(sims, "UCC_stdev")
    g_increase = create_group(sims, "addition")


## --- Run ScreenBulkBase restrictions
    for i in eachindex(init)
        # @info "Starting $(init[i][1])% volatiles simulation."

        # Initialize and call file
        global fileout = stem * init[i][1] * ".h5"     # Output file name 
        global dol, gyp, bas = init[i][2]              # Volatile restriction 
        include("../ScreenBulkBase.jl")

        # Save BitVector and sample increase to simout file
        a = zeros(Int64, length(t))
        a[t] .= 1
        g_bitvec["sim_$(init[i][1])"] = a
        g_increase["sim_$(init[i][1])"] = up
    end


## --- Run SampleMatchBase with simulation files 
    for i in eachindex(init)
        # Initialize and call file 
        global filemacrostrat = macrostrat_io
        global filebulk = stem * init[i][1] * ".h5"
        include("../SampleMatchBase.jl")

        # Write output to the simout file
        g_index["sim_$(init[i][1])"] = matches
        g_types["sim_$(init[i][1])"] = string.(littletypes)

        # Remove the simbulk file to avoid generating 5GB of test files
        run(`rm $filebulk`)

        # Compute composition of upper continental crust
        s = matches .!= 0;
        mbulk = NamedTuple{keys(bulk)}([bulk[k][matches[s]] for k in keys(bulk)])
        for k in allelements
            zeronan!(mbulk[k])
        end

        g_crust["sim_$(init[i][1])"] = [nanmean(mbulk[i]) for i in allelements]
        g_err["sim_$(init[i][1])"] = [nanstd(mbulk[i]) for i in allelements]
    end

    close(fid)

    # Stop timer 
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Program runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """


## --- End of file

# TO DO:
# We also test the sensitivity of the estimtate to normalization. Setting a maximum
# allowed wt.% assumed volatiles of 16% means that we only assume volatiles if the 
# sample would already be allowed through the filter. In this case, the effect is that
# the reported composition is not normalized to 100%. We compare this to a simulation
# where we do not add any assumed volatiles.