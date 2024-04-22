## -- Set up
    # Packages
    using RockWeatheringFlux
    using HDF5, DelimitedFiles
    using Plots
    # using StatsBase
    # using CurveFit; using Isoplot


## --- Load Macrostrat data
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

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
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    close(fid)
    
    # Sample age and uncertainty (5% or 50 Myr. uncertainty)
    # sampleage, ageuncert = resampling_age(
    #     mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
    #     macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
    #     uncert_rel=5, uncert_abs=50
    # )


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
    h = plot(xlabel="Bedrock Age [Ma]", ylabel="Hillslope [m/km]", framestyle=:box,
        legend=:topright,
        # ylims=(20,175)
    )
    for k in keys(simout)
        c,m,e = binmeans(simout[k][:,1], simout[k][:,2], 0, 3800, 38)
        plot!(h, c, m, ribbon=e, 
            color=colors[k], msc=:auto, markershape=:circle, label="$k"
        )
    end
    display(h)
    # savefig(h, "results/figures/ageslope_stack.png")


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

    
## --- End of file 