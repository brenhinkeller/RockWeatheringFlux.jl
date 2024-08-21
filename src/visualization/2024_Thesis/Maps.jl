## --- Set up
    # Demonstrate spatial sampling bias

    # Packages
    using CairoMakie
    using GeoMakie      # Must be @0.5.0 until they fix v0.6.0
    using ImageMagick
    using KernelDensity

    # Load data and base packages
    include("Definitions.jl")

    # Overwrite filepath 
    filepath = "results/figures/burial"

    # Load unmatched geochemical dataset
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    
## --- Spatiotemporal resample age and location of unmatched geochemical dataset 
    # Resampling weights
    k = invweight(bulk.Latitude, bulk.Longitude, bulk.Age);
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Age uncertainty (5% or 50 Ma.)
    ageuncert = nanadd.(bulk.Age_Max, .- bulk.Age_Min) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = nanmaximum([ageuncert[i], bulk.Age[i] .* 0.05, 50])
    end

    # Resample 
    data = hcat(bulk.Latitude, bulk.Longitude, bulk.Age)
    uncert = hcat(bulk.Loc_Prec, bulk.Loc_Prec, ageuncert)
    simout = bsresample(data, uncert, Int(1e5), p)
    simout = (;
        Latitude = simout[:,1],
        Longitude = simout[:,2],
        Age = simout[:,3],
    )

    simout.Age[(0 .> simout.Age) .| (simout.Age .> 3800)] .= NaN;
    simout.Latitude[(-90 .> simout.Latitude) .| (simout.Latitude .> 90)] .= NaN;
    simout.Longitude[(-180 .> simout.Longitude) .| (simout.Longitude .> 180)] .= NaN;


## --- Global distribution and age of resampled geochemical compilation
    t = rand(1:length(simout.Age), 40_000)

    f = Figure(resolution=(1500,800), fontsize=24,)
    ax = GeoAxis(f[1,1]; dest = "+proj=wintri", coastlines=true)
    h = CairoMakie.scatter!(ax, 
        simout.Longitude[t], simout.Latitude[t], 
        color=simout.Age[t], colormap=cscheme, 
        markersize = 5
    )
    Colorbar(f[1,2], label = "Age [Ma.]", 
        height = Relative(0.75), 
        limits=(0,3800), 
        colormap=cscheme
    )
    display(f)
    save("$filepath/map_resampled.pdf", f)


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


## --- End of file