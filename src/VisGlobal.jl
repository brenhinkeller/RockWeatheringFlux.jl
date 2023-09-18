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

    
## --- [DATA, PLOT] Location and age of EarthChem (bulk) points
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
        colormap=c_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("results/figures/bulk_restricted.png", f)

    # All data
    bulk = matread("data/bulk.mat")["bulk"]
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, vec(bulk.Longitude), vec(bulk.Latitude), color=vec(bulk.Age),
        colormap=c_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("results/figures/bulk_all.png", f)


## --- [DATA] Macrostrat and matched EarthChem samples
    # Indices
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Macrostrat
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(macrofid["vars"]["rocklat"])[t],
        rocklon = read(macrofid["vars"]["rocklon"])[t],
        age = read(macrofid["vars"]["age"])[t],
    )
    close(macrofid)

    # Matched EarthChem
    bulkfid = h5open("output/bulk.h5", "r")
    header = Tuple(Symbol.(read(bulkfid["bulk"]["header"])))
    data = read(bulkfid["bulk"]["data"])
    mbulk = NamedTuple{header}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(bulkfid)


## --- Location and age of bulk points matched to Macrostrat samples
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, vec(mbulk.Longitude), vec(mbulk.Latitude), 
        color=vec(mbulk.Age), colormap=c_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("results/figures/bulk_matched.png", f)


## --- Location and age of Macrostrat samples
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, 
        color=macrostrat.age, colormap=c_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("results/figures/macrostrat.png", f)


## --- Visualize the distribution of matched EarthChem samples
    # SiO₂ content
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, 
        color=mbulk.SiO2, colormap=c_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "SiO₂ [wt.%]", height = Relative(0.9))
    display(f)
    save("results/figures/globalsilica.png", f)

    # Distance between Macrostrat and matched EarthChem sample
    distance = Array{Float64}(undef, count(t), 1)
    distance .= haversine.(macrostrat.rocklat, macrostrat.rocklon, mbulk.Latitude, mbulk.Longitude)

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, 
        color=distance, colormap=c_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Distance to EarthChem sample [decimal degrees]", height = Relative(0.9))
    display(f)
    save("results/figures/distance.png", f)


## --- [DATA] SRTM15+ data in OCTOPUS basin polygons
#     # SRTM
#     fid = h5open("output/srtm15plus_maxslope.h5", "r")
#     srtm = read(fid["vars"])
#     srtm = NamedTuple{Tuple(Symbol.(keys(srtm)))}(values(srtm))
#     close(fid)

#     # Basin polygons
#     fid = h5open("output/basin_coordinates.h5", "r+")
#     g = create_group(fid, "basins")

#     p = Progress(length(keys(fid["vars"])))
#     i = 0
#     for obj in fid["vars"]
#         data = read(obj)

#         col, row = find_grid_inpolygon(srtm.x_lon_cntr, srtm.y_lat_cntr, data[:,2], data[:,1])
#         if length(row) > 0
#             val = [srtm.slope[row[i],col[i]] for i in eachindex(row,col)]
#             g[lpad(i, 4, "0")] = hcat(srtm.x_lon_cntr[col], srtm.y_lat_cntr[row], val)
#         end
#         i += 1
#         next!(p)
#     end
#     close(fid)


# ## --- OCTOPUS basins, or the most computationally intensive plot you've ever seen
#     fid = h5open("output/basin_coordinates.h5", "r")
#     npoints = 0
#     for obj in fid["basins"]
#         npoints += length(read(obj))
#     end

#     lat = Array{Float64}(undef, npoints)
#     lon = Array{Float64}(undef, npoints)
#     slope = Array{Float64}(undef, npoints)

#     start = stop = 1
#     for obj in fid["basins"]
#         data = read(obj)
#         stop = start + length(data[:,1]) - 1

#         lon[start:stop] .= data[:,1]
#         lat[start:stop] .= data[:,2]
#         slope[start:stop] .= data[:,3]

#         start = stop + 1
#     end
#     close(fid)

#     f = Figure(resolution = (1200, 600))
#     ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
#     CairoMakie.plot!(ax, lon, lat, color=slope, colormap=c_gradient, markersize = 3)
#     display(f)
#     save("results/figures/octopus_basins.pdf", f)

    

## --- End of file