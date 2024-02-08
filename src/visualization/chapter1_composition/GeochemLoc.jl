## --- Set up
    # Plot the locations of the geochemical sample dataset. This makes sure one region
    # isn't missing any samples.
    # To this end, also separate by rock type.

    # Unique packages
    using CairoMakie
    using GeoMakie      # Must be @0.5.0 until they fix v0.6.0
    using ImageMagick

    # Data and base packages 
    include("Definitions.jl")


## --- Overwrite a bunch of stuff so I can compare EarthChem and Gard 
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    fid = h5open("output/gard.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    gard = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    gard_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    # Major classes include minors
    include_minor!(bulk_cats)
    include_minor!(gard_cats)
    include_minor!(macro_cats)


## --- Main plot: global distribution and age of all geochemical samples
    f = Figure(resolution=(1200,1200), fontsize=24)
    ax1 = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true, title="EarthChem")
    h = CairoMakie.scatter!(ax1, 
        bulk.Longitude, bulk.Latitude, 
        color=bulk.Age, colormap=c_gradient, 
        markersize = 3
    )
    ax2 = GeoAxis(f[2,1]; dest = "+proj=wintri", coastlines=true, title="Gard et al., 2019")
    h = CairoMakie.scatter!(ax2, 
        gard.Longitude, gard.Latitude, 
        color=gard.Age, colormap=c_gradient, 
        markersize = 3,
        
    )
    Colorbar(f[1:2,2], h, label = "Age [Ma]", height = Relative(0.75))
    display(f)
    # save("$filepath/bulk_locations.pdf", f)


## --- Sedimentary rocks 
    f = Figure(resolution=(1200,1200), fontsize=24, title="Sedimentary Rocks")
    ax1 = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true, title="EarthChem")
    h = CairoMakie.scatter!(ax1, 
        bulk.Longitude[bulk_cats.sed], bulk.Latitude[bulk_cats.sed],
        color=bulk.Age[bulk_cats.sed], colormap=c_gradient, 
        markersize = 4,
    )
    ax2 = GeoAxis(f[2,1]; dest = "+proj=wintri", coastlines=true, title="Gard et al., 2019")
    h = CairoMakie.scatter!(ax2, 
        gard.Longitude[gard_cats.sed], gard.Latitude[gard_cats.sed], 
        color=gard.Age[gard_cats.sed], colormap=c_gradient, 
        markersize = 4,
    )
    Colorbar(f[1:2,2], h, label = "Age [Ma]", height = Relative(0.75))
    display(f)


## --- Volcanic / Plutonic 
    f = Figure(resolution=(1200,1200), fontsize=24, title="Sedimentary Rocks")
    ax1 = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true, title="EarthChem")
    h = CairoMakie.scatter!(ax1, 
        bulk.Longitude[bulk_cats.volc], bulk.Latitude[bulk_cats.volc],
        color=:red,
        markersize = 4,
    )
    h = CairoMakie.scatter!(ax1, 
        bulk.Longitude[bulk_cats.plut], bulk.Latitude[bulk_cats.plut],
        color=:darkblue,
        markersize = 4,
    )

    ax2 = GeoAxis(f[2,1]; dest = "+proj=wintri", coastlines=true, title="Gard et al., 2019")
    h = CairoMakie.scatter!(ax2, 
        gard.Longitude[gard_cats.volc], gard.Latitude[gard_cats.volc], 
        color=:red,
        markersize = 4,
    )
    h = CairoMakie.scatter!(ax2, 
        gard.Longitude[gard_cats.plut], gard.Latitude[gard_cats.plut], 
        color=:darkblue,
        markersize = 4,
    )
    display(f)


## --- Lithology distributions
    # Metamorphic rocks are only metamorphic if we cannot infer a protolith
    for type in keys(macro_cats)
        type==:met && continue
        macro_cats.met .&= .!macro_cats[type]
    end

    f = Figure(resolution=(1500,800), fontsize=24,)
    ax = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true)
    h = CairoMakie.scatter!(ax, 
        macrostrat.rocklon[macro_cats.sed], macrostrat.rocklat[macro_cats.sed],
        color=colors.sed, label="Sedimentary",
        markersize = 2,
    )
    h = CairoMakie.scatter!(ax, 
        macrostrat.rocklon[macro_cats.ign], macrostrat.rocklat[macro_cats.ign],
        color=colors.ign, label="Igneous",
        markersize = 2,
    )
    h = CairoMakie.scatter!(ax,
        macrostrat.rocklon[macro_cats.met], macrostrat.rocklat[macro_cats.met],
        color=colors.met, 
        markersize = 2, label="Metamorphic"
    )

    # Legend
    elem1 = MarkerElement(color=colors.sed, marker=:circle, markersize=35,
        points = Point2f[(0.5, 0.5)]
    )
    elem2 = MarkerElement(color=colors.ign, marker=:circle, markersize=35,
        points = Point2f[(0.5, 0.5)]
    )
    elem3 = MarkerElement(color=colors.met, marker=:circle, markersize=35,
        points = Point2f[(0.5, 0.5)]
    )
    Legend(f[1, 2], 
        [elem1, elem2, elem3], 
        ["Sedimentary", "Igneous", "Metamorphic"],
        patchsize = (35, 35), rowgap = 10, framevisible=false
    )
    display(f)
    save("$filepath_png/global_lithology.png", f)


## --- Lithology distributions in the bulk geochemical dataset 
    ## --- Lithology distributions
    # Metamorphic rocks are only metamorphic if we cannot infer a protolith
    for type in keys(gard_cats)
        type==:met && continue
        gard_cats.met .&= .!gard_cats[type]
    end

    f = Figure(resolution=(1500,800), fontsize=24,)
    ax = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true)
    h = CairoMakie.scatter!(ax, 
        gard.Longitude[gard_cats.ign], gard.Latitude[gard_cats.ign],
        color=colors.ign, label="Igneous",
        markersize = 2,
    )
    h = CairoMakie.scatter!(ax,
        gard.Longitude[gard_cats.met], gard.Latitude[gard_cats.met],
        color=colors.met, 
        markersize = 2, label="Metamorphic"
    )
    h = CairoMakie.scatter!(ax, 
        gard.Longitude[gard_cats.sed], gard.Latitude[gard_cats.sed],
        color=colors.sed, label="Sedimentary",
        markersize = 2,
    )

    # Legend
    elem1 = MarkerElement(color=colors.sed, marker=:circle, markersize=35,
        points = Point2f[(0.5, 0.5)]
    )
    elem2 = MarkerElement(color=colors.ign, marker=:circle, markersize=35,
        points = Point2f[(0.5, 0.5)]
    )
    elem3 = MarkerElement(color=colors.met, marker=:circle, markersize=35,
        points = Point2f[(0.5, 0.5)]
    )
    Legend(f[1, 2], 
        [elem1, elem2, elem3], 
        ["Sedimentary", "Igneous", "Metamorphic"],
        patchsize = (35, 35), rowgap = 10, framevisible=false
    )
    display(f)
    save("$filepath_png/global_lithology_gard.png", f)


## --- Crustal silica?
    f = Figure(resolution=(1500,800), fontsize=24,)
    ax = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true)
    h = CairoMakie.scatter!(ax, 
        macrostrat.rocklon, macrostrat.rocklat,
        color=mbulk.SiO2, colormap=:bluesreds, 
        markersize = 5, alpha=0.75
    )
    Colorbar(f[1,2], h, label = "SiO₂ [wt.%]", height = Relative(0.75))
    display(f)
    save("$filepath_png/global_silica.png", f)


## --- Erosion?
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point (do not propagate variance in slope)
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )
    rockslope = Measurements.value.(rockslope)

    # Calculate all erosion rates (mm/kyr) (propagate uncertainty)
    rock_ersn = emmkyr.(rockslope)
    rock_ersn = Measurements.value.(rock_ersn)

## --- Make figure...
    t = rock_ersn .< percentile(rock_ersn, 95);

    f = Figure(resolution=(1500,800), fontsize=24,)
    ax = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true)
    h = CairoMakie.scatter!(ax, 
        macrostrat.rocklon[t], macrostrat.rocklat[t],
        color=rock_ersn[t], colormap=c_gradient, 
        markersize = 5, alpha=0.75
    )
    Colorbar(f[1,2], h, label = "Erosion rate [m/Myr]", height = Relative(0.75))
    display(f)
    save("$filepath_png/global_slope.png", f)


## --- End of file 