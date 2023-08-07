## --- Set up
    # Computational Packages
    using StatGeochem
    using HDF5
    using LoopVectorization
    using Measurements
    using Static

    # Plotting Packages
    using CairoMakie
    using GeoMakie
    using ImageMagick

    # Local utilities
    include("utilities/Utilities.jl")

    # Location data
    data = importdataset("dist.csv", ',', importas=:Tuple)
    macrostrat = (rocklat = data.lat_m, rocklon = data.lon_m)
    dist = data.dist

    # silica data
    silica = importdataset("silica.csv", ',', importas=:Tuple)


## --- Plot distance between sample and matched point
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, color=dist,
        markersize = 3, alpha=0.5
    )
    Colorbar(f[1,2], h, label = "Arc Degrees from Matched Point", height = Relative(0.9))
    display(f)


## --- Plot silica content of matched samples
    t = silica.SiO2 .> 50

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, silica.lon[t], silica.lat[t], color=silica.SiO2[t],
        markersize = 3, alpha=0.5
    )
    Colorbar(f[1,2], h, label = "SiO₂ [wt/%]", height = Relative(0.9))
    display(f)



## --- Plot histograms of silica content by igneous rock type. Load data:
    # Run UCC calculations
    using Plots
    using StatGeochem
    using HDF5
    using LoopVectorization
    using Measurements
    using Static
    using DelimitedFiles
    using ProgressMeter
    include("utilities/Utilities.jl")
    
    # Matched samples
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Macrostrat
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocktype = read(macrofid["rocktype"])[t],
        rockname = read(macrofid["rockname"])[t],
        rockdescrip = read(macrofid["rockdescrip"])[t]
    )
    close(macrofid)
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
        macrostrat.rockdescrip
    )

    # Earthchem
    bulkfid = h5open("output/bulk.h5", "r")
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(bulkfid)

    # Resampled EarthChem
    rsign = importdataset("output/bulk_ignsilica_rs.tsv", '\t', importas=:Tuple)
    rsvolc = importdataset("output/bulk_volcsilica_rs.tsv", '\t', importas=:Tuple)
    rsplut = importdataset("output/bulk_plutsilica_rs.tsv", '\t', importas=:Tuple)


## --- Sample data plots; normalizing distributions
    # All igneous
    c, n = bincounts(bulk.SiO2[macro_cats.ign], 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="All Igneous; n = $(count(macro_cats.ign))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, 
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_ign.png")

    # Volcanic
    c, n = bincounts(bulk.SiO2[macro_cats.volc], 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Volcanic; n = $(count(macro_cats.volc))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box,
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_volc.png")

    # Plutonic
    c, n = bincounts(bulk.SiO2[macro_cats.plut], 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Plutonic; n = $(count(macro_cats.plut))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box,
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_plut.png")


## --- Resampled data plots; normalizing distributions
    # All igneous
    c, n, = bincounts(rsign.SiO2, 40, 80, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="All Igneous (resample); n = $(length(rsign.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, color=:purple, 
        linecolor=:purple, ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rsign.png")

    # Volcanic
    c, n, = bincounts(rsvolc.SiO2, 40, 80, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Volcanic (resample); n = $(length(rsvolc.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, color=:purple, 
        linecolor=:purple, ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rsvolc.png")

    # Plutonic
    c, n, = bincounts(rsplut.SiO2, 40, 80, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Plutonic (resample); n = $(length(rsplut.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, color=:purple, 
        linecolor=:purple, ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rsplut.png")


## --- Check how many indices in a given bin are from one sample
    # Given an index i and bin centers c...
    filter = macro_cats.plut
    i = findmax(n)[2]

    s = step(c)/2
    s = @. c[i]-s <= bulk.SiO2[filter] <= c[i]+s

    ind = bulkidx[t][macro_cats.plut][s]
    counts = [count(==(i), ind) for i in unique(ind)]

    f = findmax(counts)[1] / length(ind)
    @info "$(round(f*100, digits=2))% indices are from one EarthChem sample"

## --- End of File