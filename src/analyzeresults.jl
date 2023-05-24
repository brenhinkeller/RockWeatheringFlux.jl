## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using LoopVectorization

    # Local Utilities
    include("Utilities.jl")

    # Load results
    res = h5open("output/rwf_output2.h5", "r")


## --- Organize rock types and element types for easy indexing
    rocktypes = read(res["bulkrockflux"]["rocktypes"])
    idx = 1:length(rocktypes)
    rocktypes = NamedTuple{Tuple(Symbol.(rocktypes))}(idx)

    biglist = read(res["element_names"])
    idx = 1:length(biglist)
    biglist = NamedTuple{Tuple(Symbol.(biglist))}(idx)


## --- Calculate total global denundation (?) rate!
    netflux = read(res["bulkrockflux"]["val"])
    globalflux = netflux[rocktypes.sed] + netflux[rocktypes.ign] + netflux[rocktypes.met]


## --- End of File