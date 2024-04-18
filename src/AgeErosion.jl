## --- Model erosion (slope) as a function of rock age
    # TO DO: correlation plot
    # TO DO: PCA?

## -- Set up
    # Packages
    using RockWeatheringFlux
    # using StatGeochem
    using HDF5
    using DelimitedFiles
    using StatsBase
    # using CurveFit; using Isoplot

    using Plots
    # using StatsPlots
    # using CairoMakie
    # using GeoMakie
    # using ImageMagick

    # using LoopVectorization
    # using Static
    # using Measurements

    # Local utilities
    # include("utilities/Utilities.jl")


## --- Load Macrostrat data
    # Indices of matched EarthChem samples from SampleMatch.jl
    fid = readdlm(matchedbulk_io)
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    # Macrostrat data
    fid = h5open(macrostrat_io, "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
    )
    close(fid)

    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();
    

## --- Calculate erosion rate at each point of interest
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point with a known EarthChem sample
    # Modify this function to return an error as well
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )

    rockslope, rockslope_uncert = unmeasurementify(rockslope)

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)
    rock_ersn, rock_ersn_uncert = unmeasurementify(rock_ersn)


## --- Resample matched data by major rock type
    # Calculate rock age uncertainty, samples without uncert are assigned 5% rock age
    sampleage = copy(macrostrat.age)
    ageuncert = nanadd.(macrostrat.agemax, .- macrostrat.agemin) ./ 2;

    t = isnan.(ageuncert)
    ageuncert[t] .= sampleage[t] .* 0.05


    # Set up
    # simitemsout = [macrostrat.rocklat, macrostrat.rocklon, sampleage, rockslope]
    # simitemsuncert = (
    #     zeros(length(macrostrat.rocklat)),
    #     zeros(length(macrostrat.rocklon)),
    #     ageuncert,
    #     rockslope_uncert,
    # )

    nsims = Int(1e6)

    # Samples without min / max bounds are assigned an uncertainty of 5% of the rock age
    # for i in eachindex(simitemsuncert[3])
    #     simitemsuncert[3][i] = ifelse(isnan(simitemsuncert[3][i]), 0.05 * macrostrat.age[i], 
    #         simitemsuncert[3][i]
    #     )
    # end

    # Preallocate
    n_out = 2
    simout = (
        sed = Array{Float64}(undef, nsims, n_out),
        ign = Array{Float64}(undef, nsims, n_out),
        volc = Array{Float64}(undef, nsims, n_out),
        plut = Array{Float64}(undef, nsims, n_out),
    )

    # Run simulation for each rock type
    for key in keys(simout)
        s = match_cats[key]
        # k = invweight_age(sampleage[s])
        # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        p = ones(count(s))

        data = hcat(sampleage[s], rockslope[s])
        uncert = hcat(ageuncert[s],rockslope_uncert[s])
        simout[key] .= bsresample(data, uncert, nsims, p)
    end

    # # Map columns to data
    # c_lat = findfirst(x -> x==macrostrat.rocklat, simitemsout)
    # c_lon = findfirst(x -> x==macrostrat.rocklon, simitemsout)
    # c_slp = findfirst(x -> x==rockslope, simitemsout)
    # c_age = 0
    # for i in eachindex(simitemsout)
    #     c_age = i
    #     filter(!isnan, simitemsout[i]) == filter(!isnan, macrostrat.age) && return c_age
    # end

    # c_lat = 1; c_lon = 2; 
    c_slp = 2; c_age = 1;

    # # Remove any physically impossible data
    # # slope only: age and lat/lon should self-select when relevant
    # t = @. (0 < simout.sed[:,c_slp] < 1000)
    # simsed = simout.sed[t[:],:]
    
    # t = @. (0 < simout.ign[:,c_slp] < 1000)
    # simign = simout.ign[t[:],:]
    
    # t = @. (0 < simout.met[:,c_slp] < 1000)
    # simmet = simout.met[t[:],:]


## --- Everything everywhere over 3800 million years?
    # (Plot everything on the same axis)

    h = plot(xlabel="Bedrock Age [Ma]", ylabel="Hillslope [m/km]", framestyle=:box,
        legend=:topright, # yaxis=:log10,
        ylims=(20,175)
    )

    for k in keys(simout)
        c,m,e = binmeans(simout[k][:,c_age], simout[k][:,c_slp], 0, 3800, 38)
        plot!(h, c, m, ribbon=e, 
            color=colors[k], msc=:auto, markershape=:circle, label="$k"
        )
    end

    # c,m,e = binmeans(simout.sed[:,c_age], simout.sed[:,c_slp], 0, 3800, 38)
    # Plots.plot!(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
    #     label="Sed", markershape=:diamond,
    # )

    # c,m,e = binmeans(simout.ign[:,c_age], simout.ign[:,c_slp], 0, 3800, 38)
    # Plots.plot!(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
    #     label="Ign", markershape=:circle,
    # )

    # c,m,e = binmeans(simout.met[:,c_age], simout.met[:,c_slp], 0, 3800, 38)
    # Plots.plot!(c, m, yerror=e, color=:purple, lcolor=:purple, msc=:purple,
    #     label="Met", markershape=:star5,
    # )

    display(h)
    # savefig(h, "results/figures/ageslope_stack.png")


## --- Yorkfit in log space
    include("utilities/yorkfit.jl")

    # Transform to log-space
    ersn_sed, ersn_sed_e = unmeasurementify(emmkyr.(simsed[:,c_slp]))
    ersn_ign, ersn_ign_e = unmeasurementify(emmkyr.(simign[:,c_slp]))
    ersn_met, ersn_met_e = unmeasurementify(emmkyr.(simmet[:,c_slp]))

    c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    xσ = ones(length(c))
    fobj = yorkfit(collect(c), xσ, m, e)
    esed(x) = x * (fobj.slope) + (fobj.intercept)

    c,m,e = binmeans(simign[:,c_age], ersn_ign, 0, 3800, 38)
    fobj = yorkfit(collect(c), xσ, m, e)
    eign(x) = x * (fobj.slope) + (fobj.intercept)

    c,m,e = binmeans(simmet[:,c_age], ersn_met, 0, 3800, 38)
    fobj = yorkfit(collect(c), xσ, m, e)
    emet(x) = x * (fobj.slope) + (fobj.intercept)


## --- Erosion rate over time in linear space
    x = collect(0:10:3800)
    c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    hₛ = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
        label="Sed", markershape=:circle,)
    Plots.plot!(x, unmeasurementify(esed.(x))[1], label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simign[:,c_age], simign[:,c_slp], 0, 3800, 38)
    hᵢ = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
        label="Ign", ylabel="Hillslope [m/km]", markershape=:circle,)
    Plots.plot!(x, unmeasurementify(eign.(x))[1], label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simmet[:,c_age], simmet[:,c_slp], 0, 3800, 38)
    hₘ = Plots.plot(c, m, yerror=e, color=:orange, lcolor=:orange, msc=:orange,
        label="Met", xlabel="Bedrock Age [Ma]", markershape=:circle,)
    Plots.plot!(x, unmeasurementify(emet.(x))[1], label="Model", color=:black, linewidth=2)

    h = Plots.plot(hₛ, hᵢ, hₘ, layout=(3,1), size=(600, 1200), left_margin=(30, :px), 
        framestyle=:box, legend=:topright
    )
    display(h)
    

## --- Exponential fit to erosion
    ersn_sed, ersn_sed_e = unmeasurementify(emmkyr.(simsed[:,c_slp]))
    ersn_ign, ersn_ign_e = unmeasurementify(emmkyr.(simign[:,c_slp]))
    ersn_met, ersn_met_e = unmeasurementify(emmkyr.(simmet[:,c_slp]))

    x = 0:10:3800

    c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    sedexp = exp_fit(c, m)
    yₛ = @. sedexp[1] * exp(sedexp[2] * x);
    hₛ = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
        label="Sed", markershape=:circle,)
    Plots.plot!(x, yₛ, label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simign[:,c_age], ersn_ign, 0, 3800, 38)
    ignexp = exp_fit(c, m)
    yᵢ = @. ignexp[1] * exp(ignexp[2] * x);
    hᵢ = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
        label="Ign", ylabel="Erosion [m/Myr]", markershape=:circle,)
    Plots.plot!(x, yᵢ, label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simmet[:,c_age], ersn_met,  0, 3800, 38)
    metexp = exp_fit(c, m)
    yₘ = @. metexp[1] * exp(metexp[2] * x);
    hₘ = Plots.plot(c, m, yerror=e, color=:orange, lcolor=:orange, msc=:orange,
        label="Met", xlabel="Bedrock Age [Ma]", markershape=:circle,)
    Plots.plot!(x, yₘ, label="Model", color=:black, linewidth=2)

    h = Plots.plot(hₛ, hᵢ, hₘ, layout=(3,1), size=(600, 1200), left_margin=(30, :px), 
        framestyle=:box, legend=:topright,
    )
    display(h)


## --- Try to figure out why Archean samples are eroding so fast
## --- Show that this exists in the real data and is not an artifact of resampling
    c,m,e = binmeans(macrostrat.age[macro_cats.sed], rockslope[macro_cats.sed], 0, 3800, 38)
    h1 = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue, framestyle=:box,
        label="Sed",
        markershape=:circle, yaxis=:log10, legend=:topright,
    )

    c,m,e = binmeans(macrostrat.age[macro_cats.ign], rockslope[macro_cats.ign], 0, 3800, 38)
    h2 = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red, framestyle=:box,
        label="Ign", ylabel="Hillslope [m/km]", 
        markershape=:circle, yaxis=:log10, legend=:topright,
    )

    c,m,e = binmeans(macrostrat.age[macro_cats.met], rockslope[macro_cats.met], 0, 3800, 38)
    h3 = Plots.plot(c, m, yerror=e, color=:purple, lcolor=:purple, msc=:purple, framestyle=:box,
        label="Met", xlabel="Bedrock Age [Ma]",
        markershape=:circle, yaxis=:log10, legend=:topright,
    )

    h = Plots.plot(h1, h2, h3, layout=(3,1), size=(600, 1200), left_margin=(30, :px))
    display(h)
    

## --- Characterize old and young Archean rocks
    # Archean rocks younger than 3000 Ma tend to have lower slopes (< 30 m/km) while
    # rocks older than 3000 Ma tend to have higher slopes (>80 m/km).
    #
    # This is mostly true for sed and ign rocks, and less true for mets.

    # I hate vowels
    old_archn = @. macrostrat.age >= 3000;
    yng_archn = @. 2500 <= macrostrat.age < 3000;


## --- Geospatial (Where are the Archean rocks?)
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h1 = CairoMakie.scatter!(ax, macrostrat.rocklon[old_archn], macrostrat.rocklat[old_archn], 
        color=:crimson, markersize = 3,)
    elem1 = MarkerElement(color=:crimson, marker=:circle, markersize=15,
        points = Point2f[(0.5, 0.5)]
    )
    
    h2 = CairoMakie.scatter!(ax, macrostrat.rocklon[yng_archn], macrostrat.rocklat[yng_archn], 
        color=:blueviolet, markersize = 3,)
    elem2 = MarkerElement(color=:blueviolet, marker=:circle, markersize=15, 
        points = Point2f[(0.5, 0.5)]
    )

    Legend(f[1, 2], [elem1, elem2], ["> 3000 Ma", "< 3000 Ma"], patchsize = (35, 35), rowgap = 10)
    display(f)

## --- Slope of Archean rocks
    archean = old_archn .| yng_archn;

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h1 = CairoMakie.scatter!(ax, macrostrat.rocklon[archean], macrostrat.rocklat[archean], 
        color=rockslope[archean], colormap=c_gradient, markersize = 3,)
    Colorbar(f[1,2], h1, label = "Hillslope [m/km]", height = Relative(0.9))
    display(f)


## --- Temporal (Where are the Archean rocks?)
    ageuncert = (macrostrat.agemax .- macrostrat.agemin) ./ 2;


## --- Geologic province (Who are the Archean rocks?)
    old_provs = decode_find_geolprov(find_geolprov(macrostrat.rocklat[old_archn], 
        macrostrat.rocklon[old_archn]))
    yng_provs = decode_find_geolprov(find_geolprov(macrostrat.rocklat[yng_archn], 
        macrostrat.rocklon[yng_archn]))
    archeanprovs = unique([old_provs; yng_provs])

    old_provs = float.([count(x -> x==name, old_provs) for name in archeanprovs])
    yng_provs = float.([count(x -> x==name, yng_provs) for name in archeanprovs])


    # We already know most rocks are shields, so it's more useful to look at, for each
    # province, the proportion of rocks greater / less than 3000 Ma.
    totalprovs = old_provs .+ yng_provs
    old_provs ./= totalprovs
    yng_provs ./= totalprovs

    x = 1:length(totalprovs)
    h = StatsPlots.groupedbar([yng_provs old_provs], bar_position=:stack,
        framestyle=:box, label=["< 3000 Ma" "> 3000 Ma"],
        ylabel="Abundance", xlabel="Geologic Province", xticks=(x, archeanprovs), 
        xrotation = 45, ylims = (0, 1.1),
        legend=:outertopright, bottom_margin=(30, :px)
    )


## --- Rock type (Who are the Archean rocks?)


## --- Look at just igneous rocks
    c,m,e = binmeans(macrostrat.age[macro_cats.ign], rockslope[macro_cats.ign], 2500, 3800, 13)
    h = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red, framestyle=:box,
        label="Ign", ylabel="Hillslope [m/km]", xlabel="Age [Ma]", seriestype=:scatter,
        markershape=:circle, yaxis=:log10, legend=:topright,
    )


## --- Average slope of each province
    allprov = decode_find_geolprov(find_geolprov(macrostrat.rocklat, macrostrat.rocklon))
    uniqueprovs = unique(allprov)
    inprov = NamedTuple{Tuple(Symbol.(uniqueprovs))}([allprov .== name for name in uniqueprovs])
    
    avg_slope = [nanmean(rockslope[t]) for t in inprov]
    # TO DO: errors aren't behaving as expected--e.g. maximum error isn't plotted at 317 even though
    # that's the maximum upper error :(
    # lower = avg_slope .- [percentile(rockslope[t], 5) for t in inprov]
    # upper = [percentile(rockslope[t], 95) for t in inprov] .- avg_slope

    # Average slope by geologic province
    x = 1:length(avg_slope)
    h = Plots.plot(x, avg_slope, seriestype=:bar, framestyle=:box, label="", 
        ylabel="Average Slope [m/km]", xlabel="Geologic Province", xticks=(x, uniqueprovs), 
        xrotation = 45, ylims = (0, maximum(upper) + 0.1*maximum(upper)),
        bottom_margin=(30, :px)
    )
    display(h)


## --- Abundance of each geologic eon in each province
    allprov = decode_find_geolprov(find_geolprov(macrostrat.rocklat, macrostrat.rocklon))
    uniqueprovs = unique(allprov)

    # TO DO: normalize to relative abundance relative to abundance of all rock ages
    # Or something to show how Archean rocks tend to be in shields...
    archean = @. macrostrat.age >= 2500;
    proterozoic = @. 541 <= macrostrat.age < 2500;
    phanerozoic = @. macrostrat.age < 541;

    # Total number of rocks by eon
    # count_archean = 
    # count_proterozoic = 
    # count_phanerozoic = 

    # Abundance by province
    arc_provs = [count(x -> x==name, allprov[archean]) for name in uniqueprovs]
    pro_provs = [count(x -> x==name, allprov[proterozoic]) for name in uniqueprovs]
    pha_provs = [count(x -> x==name, allprov[phanerozoic]) for name in uniqueprovs]
    total_provs = @. arc_provs + pro_provs + pha_provs

    x = 1:length(arc_avg_slope)
    h = StatsPlots.groupedbar([arc_provs pro_provs pha_provs], bar_position=:stack,
        framestyle=:box, label=["Archean" "Proterozoic" "Phanerozoic"],
        ylabel="Abundance", xlabel="Geologic Province", xticks=(x, uniqueprovs), 
        xrotation = 45, ylims = (0, maximum(total_provs) + 0.1*maximum(total_provs)),
        legend=:topright, bottom_margin=(30, :px)
    )

    
## --- End of file 