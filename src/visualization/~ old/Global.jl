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
    include("../utilities/Utilities.jl")

    
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


## --- EarthChem points by rock type
    fid = h5open("output/bulk.h5", "r")
        header = read(fid["bulk"]["header"])
        data = read(fid["bulk"]["data"])
        bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

        header = read(fid["bulktypes"]["bulk_cats_head"])
        data = read(fid["bulktypes"]["bulk_cats"])
        data = @. data > 0
        bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    minorsed, minorign, minormet = get_minor_types()
    for type in minorsed
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end
    for type in minormet
        bulk_cats.met .|= bulk_cats[type]
    end
    
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")

    # catkeys = reverse(collect(keys(bulk_cats)))
    catkeys = [:ign, :sed, :met]
    for i in eachindex(catkeys)
        CairoMakie.scatter!(ax, bulk.Longitude[bulk_cats[catkeys[i]]], 
            bulk.Latitude[bulk_cats[catkeys[i]]], color=colors[catkeys[i]], 
            markersize = 3, label="$(catkeys[i])"
        )
    end

    Legend(f[1, 2], ax)
    display(f)
    

## --- [DATA] Macrostrat and matched EarthChem samples
    # Indices
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")[:,1]))
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


## --- Locations of points matched with a specific sample
    # TEMPORARY: Allows this block to be run from SampleMatch.jl
    # sampletypes_temp = string.(vec(readdlm("$matchedbulk_io")[:,2]))
    bulkidx = matches;

    # Semi-temporary?
    target = :shale             # Target rock type
    i = 142971                   # Target index

    # sample_cats = match_rocktype(sampletypes_temp)
    sample_cats = match_rocktype(string.(sampletypes))
    typefilter = sample_cats[target]

    # Filters
    t = @. bulkidx != 0;
    s = @. bulkidx[t] == i;

    dist = haversine.(macrostrat.rocklat[typefilter], macrostrat.rocklon[typefilter],
        bulk.Latitude[bulkidx[t]], bulk.Longitude[bulkidx[t]]
    )

    # Make plot
    f = Figure(resolution = (1400, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h1 = CairoMakie.scatter!(ax, macrostrat.rocklon[typefilter], macrostrat.rocklat[typefilter], 
        color=dist, markersize = 3, alpha=0.1, label="All Macrostrat $target"
    )
    h2 = CairoMakie.scatter!(ax, macrostrat.rocklon[typefilter][s], macrostrat.rocklat[typefilter][s], 
        color=dist[s], markersize = 3, label="Matched with sample $i"
    )
    h3 = CairoMakie.scatter!(ax, bulk.Longitude[bulk_cats[target]], bulk.Latitude[bulk_cats[target]], 
        color=:red, markersize = 5, label="All EarthChem $target"
    )
    h4 = CairoMakie.scatter!(ax, bulk.Longitude[i], bulk.Latitude[i], 
        color=:red, markersize = 20, marker=:star5, label="Sample $i"
    )
    Colorbar(f[1, 2], h2, label = "Distance to matched point [arc-degrees]", vertical = true)
    Legend(f[1, 3], ax)
    display(f)
    

## --- TEMPORARY: distance between points and their matched sample
    t = @. matches != 0;
    dist = haversine.(macrostrat.rocklat[typefilter], macrostrat.rocklon[typefilter],
        bulk.Latitude[matches[t]], bulk.Longitude[matches[t]]
    )
    f = Figure(resolution = (1400, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h1 = CairoMakie.scatter!(ax, macrostrat.rocklon[typefilter], macrostrat.rocklat[typefilter], 
        color=dist, markersize = 3, label="Macrostrat $(target)s"
    )
    h2 = CairoMakie.scatter!(ax, bulk.Longitude[bulk_cats[target]], bulk.Latitude[bulk_cats[target]], 
        color=:red, markersize = 3, label="EarthChem $(target)s"
    )
    Colorbar(f[1, 2], h1, label = "Distance to matched point [arc-degrees]", vertical = true)
    Legend(f[1, 3], ax)
    display(f)


## --- [DATA] Slope at each point
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )
    rock_ersn = emmkyr.(rockslope)
    
    slope = unmeasurementify(rockslope)[1]
    ersn = unmeasurementify(rock_ersn)[1]


## --- Slope of rocks older than 3000 Ma
    t = @. macrostrat.age > 3000;

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon[t], macrostrat.rocklat[t], 
        color=slope[t], colormap=c_gradient, markersize = 5
    )
    Colorbar(f[1,2], h, label = "Hillslope [m/km]", height = Relative(0.9))
    display(f)


## --- Slope of every point on Earth
    # TO DO: I should not have to restrict slope like this... something is wrong I think
    # t = @. slope < 1000
    t = @. slope < 600
    # t = trues(length(slope))

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon[t], macrostrat.rocklat[t], 
        color=ersn[t], colormap=c_gradient, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Erosion [m/Myr]", height = Relative(0.9))
    display(f)
    # save("results/figures/slope.png", f)


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