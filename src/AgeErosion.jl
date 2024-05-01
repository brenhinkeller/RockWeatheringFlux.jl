## -- Set up
    # Packages
    using RockWeatheringFlux
    using HDF5, DelimitedFiles
    using Plots
    using Isoplot: yorkfit
    # using StatsBase
    # using CurveFit; using Isoplot


## --- Load Macrostrat data
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

    # Indices of matched EarthChem samples from SampleMatch.jl
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Macrostrat data
    fid = h5open(macrostrat_io, "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    close(fid)

    # Geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)


## --- Slope / erosion rate at each coordinate point
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point (exclude slope > 1000 m/km)
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )
    rockslope = (;
        vals = Measurements.value.(rockslope),
        errs = Measurements.uncertainty.(rockslope),
    )
    t = rockslope.vals .>= 1000;
    rockslope.vals[t] = rockslope.errs[t] .= NaN

    # Calculate erosion rate [mm/kyr] (exclude erosion > 10_000 mm/kyr)
    rock_ersn = emmkyr.(rockslope.vals)
    rock_ersn = (
        vals = Measurements.value.(rock_ersn),
        errs = Measurements.uncertainty.(rock_ersn),
    )
    t = rock_ersn.vals .> 10_000
    rock_ersn.vals[t] = rock_ersn.errs[t] .= NaN


## --- Resample matched data by major rock type
    # Calculate rock age and uncertainty, using mapped age because we're basically 
    # interested in the relationship between location and age.
    # Rock uncertainty is the maximum of mapped uncertainty, 5% age, or 50 Ma.
    sampleage = copy(macrostrat.age)
    ageuncert = nanadd.(macrostrat.agemax, .- macrostrat.agemin) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = nanmaximum([ageuncert[i], sampleage[i] .* 0.05, 50])
    end

    # Preallocates
    nsims = Int(1e6)
    target = (:sed, :ign, :volc, :plut)
    simout = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

    # Run simulation for each rock type
    for key in keys(simout)
        s = match_cats[key]
        # k = invweight_age(sampleage[s])
        # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        p = ones(count(s))

        data = hcat(sampleage[s], rock_ersn.vals[s])
        uncert = hcat(ageuncert[s],rock_ersn.errs[s])
        simout[key] .= bsresample(data, uncert, nsims, p)
    end


## --- Plot erosion binned by age
    h = plot(xlabel="Bedrock Age [Ma]", ylabel="Erosion [m/Myr.]", framestyle=:box,
        legend=:topright,
        ylims=(0, 240),
        fg_color_legend=:white,
    )
    for k in keys(simout)
        c,m,e = binmeans(simout[k][:,1], simout[k][:,2], 0, 3800, 38)
        plot!(h, c, m, yerror=2e, 
            label="$k",
            color=colors[k], linecolor=colors[k], msc=:auto, 
            markershape=:circle, 
            
        )
    end
    display(h)
    savefig(h, "results/figures/burial/age_erosion.pdf")


## --- Sanity check volcanic results
    # Separate felsic / mafic volcanics
    comp = (;
        felsic = (62 .<= mbulk.SiO2 .< 73) .& match_cats.volc,
        mafic = (43 .<= mbulk.SiO2 .< 62) .& match_cats.volc,
    );

    # Resample 
    simvolc = NamedTuple{keys(comp)}(Array{Float64}(undef, nsims, 2) for _ in keys(comp))
    for key in keys(simvolc)
        s = comp[key]
        p = ones(count(s))

        data = hcat(sampleage[s], rock_ersn.vals[s])
        uncert = hcat(ageuncert[s],rock_ersn.errs[s])
        simvolc[key] .= bsresample(data, uncert, nsims, p)
    end

    # Plot erosion by age 
    color_volc = (felsic=:royalblue, mafic=:hotpink)
    h = plot(xlabel="Bedrock Age [Ma]", ylabel="Erosion [m/Myr.]", framestyle=:box,
        legend=:topright,
    )
    for k in keys(simvolc)
        c,m,e = binmeans(simvolc[k][:,1], simvolc[k][:,2], 0, 3800, 38)
        plot!(h, c, m, # ribbon=e, 
            color=color_volc[k], msc=:auto, markershape=:circle, label="$k",
            linewidth=2,
        )
    end
    display(h)
    savefig(h, "results/figures/burial/age_erosion_volcanic.pdf")


## --- Where are the fast-eroding igneous samples? 
    # Sanity check to make sure this spike appears for realsies 
    c,m,e = binmeans(sampleage[match_cats.ign], rock_ersn.vals[match_cats.ign], 0, 3800, 38)
    t = .!isnan.(m)
    h = plot(c[t], m[t], yerror=2e[t],
        label="",
        color=colors.ign, lcolor=colors.ign, msc=:auto,
        markershape=:circle,
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="Observed Erosion Rate [m/Myr]",
    )

    # Location of the spike for volcanic and plutonic rocks
    t = @. 550 < sampleage < 850;
    s = rock_ersn.vals .> nanmean(rock_ersn.vals[match_cats.ign]);
    r = s .& t .& match_cats.volc;
    h = mapplot(mbulk.Longitude[r], mbulk.Latitude[r], 
        label="",
        color=colors.volc, msc=:auto,
        markersize=1,
    )
    display(h)

    r = s .& t .& match_cats.plut;
    h = mapplot(mbulk.Longitude[r], mbulk.Latitude[r], 
        label="",
        color=colors.plut, msc=:auto,
        markersize=1,
    )
    display(h)

    # Histograms?
    r = s .& t .& match_cats.volc;
    h = histogram(rock_ersn.vals[r])
    display(h)

    r = s .& t .& match_cats.plut;
    h = histogram(rock_ersn.vals[r])
    display(h)

## --- Model erosion as a piecewise function of age
    # # Plot base 
    # h = plot(
    #     xlabel="Bedrock Age [Ma]", ylabel="[LOG] Erosion [m/Myr.]",
    #     framestyle=:box,
    #     fontfamily=:Helvetica, 
    #     titleloc=:left,
    # )

    
## --- Plutonic
    # # Yorkfit in logspace 0 - 1600 Ma; Constant 1600 - 3800 Ma.
    # h1 = deepcopy(h)
    # c,m,e = binmeans(simout.plut[:,1], simout.plut[:,2], 0, 3800, 38)
    # plot!(h1, c, log10.(m), label="",
    #     color=colors.plut, msc=:auto, markershape=:circle,
    #     title="Plutonic",
    # )

    # # Yorkfit 0 - 1600 Ma.
    # c,m,e = binmeans(simout.plut[:,1], simout.plut[:,2], 0, 1600, 16)
    # c_err = fill(100, length(c))
    # fobj = yorkfit(collect(c), c_err, log10.(m), log10.(e))
    
    # x = 0:50:1550
    # plot!(x, fobj.slope.val .* x .+ fobj.intercept.val, label="",
    #     color=:black,
    # )

    # # Constant 1600 - 3800 Ma.
    # c,m,e = binmeans(simout.plut[:,1], simout.plut[:,2], 1600, 3800, Int((3800-1600)/100))
    # m = nanmean(log10.(m))

    # x = 1550:50:3750
    # plot!(x, fill(m, length(x)), label="",
    #     color=:black,
    # )
    

## --- Volcanic
    # # Yorkfit 0 - 900 Ma.; Yorkfit 900 - 1700 Ma.; Constant 1700 - 3800 Ma.
    # h1 = deepcopy(h)
    # c,m,e = binmeans(simout.volc[:,1], simout.volc[:,2], 0, 3800, 38)
    # plot!(h1, c, log10.(m), label="",
    #     color=colors.volc, msc=:auto, markershape=:circle,
    #     title="Volcanic",
    # )

    # # Yorkfit 0 - 900 Ma.
    

    # # Yorkfit 900 - 1700 Ma.
    
    
    # # Constant 1700 - 3800 Ma.


## --- Yorkfit in log space
    # include("utilities/yorkfit.jl")

    # # Transform to log-space
    # ersn_sed, ersn_sed_e = unmeasurementify(emmkyr.(simsed[:,c_slp]))
    # ersn_ign, ersn_ign_e = unmeasurementify(emmkyr.(simign[:,c_slp]))
    # ersn_met, ersn_met_e = unmeasurementify(emmkyr.(simmet[:,c_slp]))

    # c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    # xσ = ones(length(c))
    # fobj = yorkfit(collect(c), xσ, m, e)
    # esed(x) = x * (fobj.slope) + (fobj.intercept)

    # c,m,e = binmeans(simign[:,c_age], ersn_ign, 0, 3800, 38)
    # fobj = yorkfit(collect(c), xσ, m, e)
    # eign(x) = x * (fobj.slope) + (fobj.intercept)

    # c,m,e = binmeans(simmet[:,c_age], ersn_met, 0, 3800, 38)
    # fobj = yorkfit(collect(c), xσ, m, e)
    # emet(x) = x * (fobj.slope) + (fobj.intercept)


## --- Erosion rate over time in linear space
    # x = collect(0:10:3800)
    # c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    # hₛ = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
    #     label="Sed", markershape=:circle,)
    # Plots.plot!(x, unmeasurementify(esed.(x))[1], label="Model", color=:black, linewidth=2)

    # c,m,e = binmeans(simign[:,c_age], simign[:,c_slp], 0, 3800, 38)
    # hᵢ = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
    #     label="Ign", ylabel="Hillslope [m/km]", markershape=:circle,)
    # Plots.plot!(x, unmeasurementify(eign.(x))[1], label="Model", color=:black, linewidth=2)

    # c,m,e = binmeans(simmet[:,c_age], simmet[:,c_slp], 0, 3800, 38)
    # hₘ = Plots.plot(c, m, yerror=e, color=:orange, lcolor=:orange, msc=:orange,
    #     label="Met", xlabel="Bedrock Age [Ma]", markershape=:circle,)
    # Plots.plot!(x, unmeasurementify(emet.(x))[1], label="Model", color=:black, linewidth=2)

    # h = Plots.plot(hₛ, hᵢ, hₘ, layout=(3,1), size=(600, 1200), left_margin=(30, :px), 
    #     framestyle=:box, legend=:topright
    # )
    # display(h)
    

## --- Exponential fit to erosion
    # ersn_sed, ersn_sed_e = unmeasurementify(emmkyr.(simsed[:,c_slp]))
    # ersn_ign, ersn_ign_e = unmeasurementify(emmkyr.(simign[:,c_slp]))
    # ersn_met, ersn_met_e = unmeasurementify(emmkyr.(simmet[:,c_slp]))

    # x = 0:10:3800

    # c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    # sedexp = exp_fit(c, m)
    # yₛ = @. sedexp[1] * exp(sedexp[2] * x);
    # hₛ = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
    #     label="Sed", markershape=:circle,)
    # Plots.plot!(x, yₛ, label="Model", color=:black, linewidth=2)

    # c,m,e = binmeans(simign[:,c_age], ersn_ign, 0, 3800, 38)
    # ignexp = exp_fit(c, m)
    # yᵢ = @. ignexp[1] * exp(ignexp[2] * x);
    # hᵢ = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
    #     label="Ign", ylabel="Erosion [m/Myr]", markershape=:circle,)
    # Plots.plot!(x, yᵢ, label="Model", color=:black, linewidth=2)

    # c,m,e = binmeans(simmet[:,c_age], ersn_met,  0, 3800, 38)
    # metexp = exp_fit(c, m)
    # yₘ = @. metexp[1] * exp(metexp[2] * x);
    # hₘ = Plots.plot(c, m, yerror=e, color=:orange, lcolor=:orange, msc=:orange,
    #     label="Met", xlabel="Bedrock Age [Ma]", markershape=:circle,)
    # Plots.plot!(x, yₘ, label="Model", color=:black, linewidth=2)

    # h = Plots.plot(hₛ, hᵢ, hₘ, layout=(3,1), size=(600, 1200), left_margin=(30, :px), 
    #     framestyle=:box, legend=:topright,
    # )
    # display(h)

    
## --- End of file 