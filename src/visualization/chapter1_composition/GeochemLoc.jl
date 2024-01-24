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

    # Major classes include minors
    include_minor!(bulk_cats)
    include_minor!(gard_cats)


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


## --- End of file 