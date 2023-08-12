# Visualize global data; maps.
#
# CTRL+K CTR+[ to collapse all sections

## --- Set up
    # Base computational packages
    using StatGeochem
    using HDF5
    using LoopVectorization
    using Measurements
    using Static
    using DelimitedFiles
    using ProgressMeter

    using CairoMakie
    using GeoMakie
    using ImageMagick

    # Local utilities
    include("utilities/Utilities.jl")

    
## --- Location and age of EarthChem (bulk) points
    # Additional packages
    using MAT

    # From wt.% restricted to 84-104 total wt.% 
    bulkfid = h5open("output/bulk.h5", "r")
        header = read(bulkfid["bulk"]["header"])
        data = read(bulkfid["bulk"]["data"])
    close(bulkfid)
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, bulk.Longitude, bulk.Latitude, color=bulk.Age,
        colormap=clr_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("bulk_restricted.png", f)

    # All data
    bulk = matread("data/bulk.mat")["bulk"]
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, vec(bulk.Longitude), vec(bulk.Latitude), color=vec(bulk.Age),
        colormap=clr_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("bulk_all.png", f)


## --- Location and age of bulk points matched to Macrostrat samples
    # Get indices of matched samples
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Get matched Earthchem (bulk) samples
    bulkfid = h5open("output/bulk.h5", "r")
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(bulkfid)

    # Plot
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, vec(bulk.Longitude), vec(bulk.Latitude), color=vec(bulk.Age),
        colormap=clr_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("bulk_matched.png", f)


## --- Location and age of Macrostrat samples
    # Get indices of matched samples
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Get Macrostrat samples with a known match
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(macrofid["rocklat"])[t],
        rocklon = read(macrofid["rocklon"])[t],
        age = read(macrofid["age"])[t]
    )
    close(macrofid)

    # Plot
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, 
        color=macrostrat.age, colormap=clr_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("macrostrat.png", f)


## --- [DATA MISSING] Difference in SiO₂ and location between bulk and Macrostrat
    # # Location data
    # data = importdataset("dist.csv", ',', importas=:Tuple)
    # macrostrat = (rocklat = data.lat_m, rocklon = data.lon_m)
    # dist = data.dist

    # # silica data
    # silica = importdataset("silica.csv", ',', importas=:Tuple)

    # # Distance between sample and matched point
    # f = Figure(resolution = (1200, 600))
    # ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    # h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, color=dist,
    #     markersize = 3, alpha=0.5
    # )
    # Colorbar(f[1,2], h, label = "Arc Degrees from Matched Point", height = Relative(0.9))
    # display(f)

    # # Silica content of matched samples
    # t = silica.SiO2 .> 50

    # f = Figure(resolution = (1200, 600))
    # ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    # h = CairoMakie.scatter!(ax, silica.lon[t], silica.lat[t], color=silica.SiO2[t],
    #     markersize = 3, alpha=0.5
    # )
    # Colorbar(f[1,2], h, label = "SiO₂ [wt/%]", height = Relative(0.9))
    # display(f)


## --- End of file