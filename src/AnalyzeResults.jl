## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using DelimitedFiles
    using Measurements
    using Plots

    # Local Utilities
    include("Utilities.jl")
    include("NaNMeasurements.jl")

    # Conversion for kg to Gt. Data file is in units of kg/yr
    const kg_gt = 1e12


## --- Load and parse data
    # Macrostrat
    # macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)
    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    # Bulk denundation
    res = h5open("output/erodedmaterial.h5", "r")
    path = res["bulk_denundation"]
    bulk_denundation = read(path["values"]) .± read(path["errors"])

    # Element flux
    path = res["element_flux"]
    vals = read(path["values"])
    errs = read(path["errors"])
    header = read(path["header"])
    close(res)

    elementflux = NamedTuple{Tuple(Symbol.(header))}(
        [zeros(Measurement{Float64}, size(vals)[1]) for _ in 1:(size(vals)[2])]
    )
    for i in eachindex(header)
        elementflux[Symbol(header[i])] .= vals[:,i] .± errs[:,i]
    end

    
## --- Compute total...
    # Global denundation
    global_denun = nansum(bulk_denundation) / kg_gt
    @info "Total global denundation: $(round(global_denun, digits=2)) Gt/yr"

    # Global flux for each element
    global_elemflux = Dict{Symbol, Measurement{Float64}}()
    globalsum = 0.0
    for i in keys(elementflux)
        global_elemflux[i] = nansum(elementflux[i]) / kg_gt
        global globalsum += global_elemflux[i]
    end
    global_elemflux = NamedTuple{Tuple(keys(global_elemflux))}(values(global_elemflux))

    @info "Total element flux: $(round(globalsum, digits=2)) Gt/yr"

## --- End of File