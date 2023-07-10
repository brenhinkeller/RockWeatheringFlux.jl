## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using DelimitedFiles
    using Measurements
    using Plots

    # Local Utilities
    include("Utilities.jl")

    # Conversion for kg to Gt. Data file is in units of kg/yr
    const kg_gt = 1e12


## --- Load and parse data
    # Indices of matched samples from SampleMatch.jl
    bulkidx = Int.(vec(readdlm("output/bulkidx.tsv")))
    t = @. bulkidx != 0     # Exclude samples with missing data

    # Macrostrat, if there's matched EarthChem data
    macrostrat = importdataset("output/pregenerated_responses.tsv", '\t', importas=:Tuple)
    macro_cats = match_rocktype(macrostrat.rocktype[t], macrostrat.rockname[t], 
        macrostrat.rockdescrip[t]
    )

    # Bulk denundation at each point
    res = h5open("output/erodedmaterial.h5", "r")
    path = res["bulk_denundation"]
    bulk_denundation = read(path["values"]) .± read(path["errors"])

    # Global flux of each element at each point
    path = res["element_flux"]
    vals = read(path["values"])
    errs = read(path["errors"])
    header = read(path["header"])
    close(res)

    elementflux = NamedTuple{Tuple(Symbol.(header))}(
        [vals[:,i] .± errs[:,i] for i in eachindex(header)]
    )

    
## --- Compute total...
    # Global denundation
    global_denun = nansum(bulk_denundation) / kg_gt
    
    # Total global flux for each element
    totalelementflux = NamedTuple{Tuple(keys(elementflux))}([nansum(i) / kg_gt for i in elementflux])

    # Sum of all element fluxes
    globalflux = nansum(totalelementflux)

    # Print to terminal
    @info """
    Total global denundation: $global_denun, digits=3)) Gt/yr
    Sum of all element fluxes: $(globalsum[1]), digits=3)) Gt/yr
    """


## --- Export results


## --- End of File