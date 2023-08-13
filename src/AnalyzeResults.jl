## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using DelimitedFiles
    using Measurements
    using Plots
    using LoopVectorization
    using Static

    # Local Utilities
    include("utilities/Utilities.jl")

    # Conversion for kg to Gt. Data file is in units of kg/yr
    const kg_gt = 1e12


## --- Load and parse data
    # Indices of matched samples from SampleMatch.jl
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0     # Exclude samples with missing data

    # Macrostrat, if there's matched EarthChem data
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        type = read(macrofid["typecategory"])[t],
    )
    close(macrofid)
    macro_cats = match_rocktype(macrostrat.type)
    subcats = collect(keys(macro_cats))
    deleteat!(subcats, findall(x->x==:cover,subcats))

    # Bulk denundation at each point in kg/yr
    res = h5open("$eroded_out", "r")
    path = res["bulk_denundation"]
    bulk_denundation = read(path["values"]) .± read(path["errors"])

    # Global flux of each element at each point in kg/yr
    path = res["element_flux"]
    vals = read(path["values"])
    errs = read(path["errors"])
    header = Symbol.(read(path["header"]))
    close(res)

    elementflux = NamedTuple{Tuple(header)}(
        [vals[:,i] .± errs[:,i] for i in eachindex(header)]
    )

    
## --- Compute total...
    # Global denundation in Gt/yr
    global_denun = nansum(bulk_denundation) / kg_gt
    
    # Total global flux for each element in Gt/yr
    totalelementflux = NamedTuple{Tuple(keys(elementflux))}([nansum(i) / kg_gt for i in elementflux])

    # Sum of all element fluxes in Gt/yr
    globalflux = nansum(totalelementflux)

    # Print to terminal
    @info """
    Total global denundation: $global_denun Gt/yr
    Sum of all element fluxes: $globalflux Gt/yr
    """


## --- Export results
    # Preallocate
    result = zeros(length(header), length(macro_cats))
    rows = string.(header)
    cols = hcat("", string.(reshape(subcats, 1, length(subcats))), "global")

    # Absolute contribution of each rock type to element flux
    for i in eachindex(header)          # Filter by element (row)
        for j in eachindex(subcats)     # Filter by rock type (column)
            result[i,j] = nansum(elementflux[i][macro_cats[subcats[j]]]).val / kg_gt
        end
    end

    TEFval = unmeasurementify(totalelementflux)[1]
    result[:,end] = TEFval

    # Absolute contribution
    writedlm("$erodedabs_out", vcat(cols, hcat(rows, result)))

    # Relative contribution
    writedlm("$erodedrel_out", vcat(cols, hcat(rows, result ./ TEFval)))


## --- End of File