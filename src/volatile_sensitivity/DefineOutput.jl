## Define the simulation output file 

## --- Set up
    using RockWeatheringFlux
    using HDF5
    
    # Parse file name
    simout = ARGS[1]


## --- Create file 
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
        g_types = create_group(sims, "class_assigned")
        g_crust = create_group(sims, "UCC")
        g_err = create_group(sims, "UCC_std")
        g_increase = create_group(sims, "addition")

    close(fid)


## --- End of file