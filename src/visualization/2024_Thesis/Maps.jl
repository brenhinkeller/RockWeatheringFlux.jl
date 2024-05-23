## --- Set up
    # Demonstrate spatial sampling bias

    # Packages
    using CairoMakie
    using GeoMakie      # Must be @0.5.0 until they fix v0.6.0
    using ImageMagick
    using KernelDensity

    # Load data and base packages
    include("Definitions.jl")

    # Load unmatched geochemical dataset
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    
## --- Resample ages of matched geochemical dataset 
    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
        uncert_rel=0.05, uncert_abs=50
    )
    k = invweight_age(sampleage)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    sim_age = bsresample(sampleage, ageuncert, Int(1e6), p)


## --- Global distribution and age of geochemical compilation 
    t = rand(1:length(bulk.Age), 40_000)

    f = Figure(resolution=(1500,800), fontsize=24,)
    ax = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true)
    h = CairoMakie.scatter!(ax, 
        bulk.Longitude[t], bulk.Latitude[t], 
        color=bulk.Age[t], colormap=cscheme, 
        markersize = 5
    )
    Colorbar(f[1,2], label = "Age [Ma.]", 
        height = Relative(0.75), 
        limits=(0,3800), 
        colormap=cscheme
    )
    display(f)
    save("$filepath/map_geochemical.pdf", f)


## --- Global distribution and age of mapped samples 
    t = rand(1:length(macrostrat.age), 40_000)    

    f = Figure(resolution=(1500,800), fontsize=24,)
    ax = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true)
    h = CairoMakie.scatter!(ax, 
        macrostrat.rocklon[t], macrostrat.rocklat[t], 
        color=macrostrat.age[t], colormap=cscheme, 
        markersize = 5
    )
    Colorbar(f[1,2], label = "Age [Ma.]", 
        height = Relative(0.75), 
        limits=(0,3800), 
        colormap=cscheme
    )
    display(f)
    save("$filepath/map_mapped.pdf", f)


## --- Set up histograms 
    # Define a color scheme!
    pal2 = (
        # bulk = pal[2],
        mbulk = pal[2],
        macrostrat = pal[6],
        sim = pal[9],
    )

    # Bins and bars 
    agemin, agemax, agebins = 0,3800,38
    latmin, latmax, latbins = -90, 90, 18
    lonmin, lonmax, lonbins = -180, 180, 36

    agebar = ((agemax-agemin)/agebins)
    latbar = ((latmax-latmin)/latbins)
    lonbar = ((lonmax-lonmin)/lonbins)
    

## --- Age 
    h = Plots.plot(
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="Relative Abundance",
        grid=false,
        fg_color_legend=:white,
    )

    # Mapped samples
    c,n = bincounts(macrostrat.age, agemin, agemax, agebins)
    n = float(n) ./ nansum(float(n) .* step(c))
    u = kde(macrostrat.age[.!isnan.(macrostrat.age)])
    Plots.plot!(u.x, u.density, label="Mapped", 
        color=pal2.macrostrat,
        linewidth=2,
    )

    # Matched samples
    c,n = bincounts(mbulk.Age, agemin, agemax, agebins)
    n = float(n) ./ nansum(float(n) .* step(c))
    u = kde(mbulk.Age[.!isnan.(mbulk.Age)])
    Plots.plot!(u.x, u.density, label="Matched", 
        color=pal2.mbulk,
        linewidth=2,
    )
    ymin, ymax = ylims(h)
    Plots.xlims!(agemin, agemax)

    # Temporal resample
    c,n = bincounts(sim_age, agemin, agemax, agebins)
    n = float(n) ./ nansum(float(n) .* step(c))
    u = kde(sim_age[.!isnan.(sim_age) .& (sim_age .> 0) .& (sim_age .< 3800)])
    Plots.plot!([-1], [-1], label="Resampled", 
        color=pal2.sim,
        linewidth=2,
    )
    Plots.ylims!(ymin, ymax)

    Plots.plot!(twinx(), u.x, u.density, label="", 
        color=pal2.sim,
        linewidth=2,
        ylabel="Relative Abundance"
    )
    Plots.xlims!(agemin, agemax)

    display(h)
    savefig(h, "$filepath/map_hist_age.pdf")
    

## --- Latitude 
    h = Plots.plot(
        framestyle=:box,
        xlabel="Latitude", ylabel="Abundance",
        grid=false,
        fg_color_legend=:white,
    )

    c,n = bincounts(macrostrat.age, agemin, agemax, agebins)
    n = float(n) ./ nansum(float(n) .* step(c))
    u = kde(macrostrat.age[.!isnan.(macrostrat.age)])
    Plots.plot!(u.x, u.density, label="Mapped", 
        color=pal2.macrostrat,
        linewidth=2,
    )

    c,n = bincounts(mbulk.Age, agemin, agemax, agebins)
    n = float(n) ./ nansum(float(n) .* step(c))
    u = kde(mbulk.Age[.!isnan.(mbulk.Age)])
    Plots.plot!(u.x, u.density, label="Matched", 
        color=pal2.mbulk,
        linewidth=2,
    )

    display(h)
    savefig(h, "$filepath/map_hist_lat.pdf")


## --- Longitude

## --- End of file