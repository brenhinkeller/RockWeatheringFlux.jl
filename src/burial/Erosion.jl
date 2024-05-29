## --- Set up 
    # Model erosion (slope) as a function of age 

    # Packages 
    using RockWeatheringFlux
    using HDF5, DelimitedFiles
    using Plots
    using LsqFit: curve_fit
    using Distributions

    # Save figures to: 
    filepath = "results/figures/burial"

    # Definitions
    nsims = Int(1e6)
    xmin, xmax, nbins = 0, 3800, 38
    age_error = 0.05
    age_error_abs = 50


## --- Load data 
    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();
    classes = keys(match_cats)

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)

    # Macrostrat
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    close(fid)
    
    # SRTM15+ DEM 
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")


## --- Calculate mean slope at each spatial point 
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )
    rockslope = (;
        vals = Measurements.value.(rockslope),
        errs = Measurements.uncertainty.(rockslope),
    )
    t = rockslope.vals .>= 1000;
    rockslope.vals[t] = rockslope.errs[t] .= NaN

    # Calculate erosion rate [mm/kyr] and exclude erosion > 10_000 mm/kyr
    rock_ersn = emmkyr.(rockslope.vals)
    rock_ersn = (
        vals = Measurements.value.(rock_ersn),
        errs = Measurements.uncertainty.(rock_ersn),
    )
    t = rock_ersn.vals .> 10_000
    rock_ersn.vals[t] = rock_ersn.errs[t] .= NaN    


## --- Resample age / erosion relationship 
    # Ages and uncertainties (prefer mapped age since we're considering spatial relationships)
    # sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
    #     macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
    #     uncert_rel=age_error, uncert_abs=age_error_abs
    # )
    sampleage = copy(macrostrat.age)
    ageuncert = nanadd.(macrostrat.agemax, .- macrostrat.agemin) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = nanmaximum([ageuncert[i], sampleage[i] .* age_error, age_error_abs])
    end

    # Rock class assigned during matching 
    matched_in = Array{Int64}(undef, length(match_cats.sed), length(classes))
    for i in eachindex(classes)
        for j in eachindex(match_cats[classes[i]])
            matched_in[j,i] = ifelse(match_cats[classes[i]][j], 1, 0)
        end
    end

    # Resampling weights are equal since we're only concerned with what's exposed on the 
    # surface right now 
    p = ones(length(sampleage))

    # Resample! 
    data = hcat(sampleage, rock_ersn.vals, matched_in)
    uncert = hcat(ageuncert, rock_ersn.errs, zeros(size(matched_in)))
    resampled = bsresample(data, uncert, nsims, p)

    # Parse resampled data into usable Tuples 
    sim_ersn = (;
        age = resampled[:,1],
        ersn = resampled[:,2]
    )
    sim_cats = delete_cover(get_cats(false, nsims)[2]);
    i = 3
    for k in classes
        sim_cats[k] .= resampled[:,i] .> 0
        global i += 1
    end


## --- [PLOT] Age / erosion relationships 
    h = plot(
        xlabel="Age [Ma.]", ylabel="Erosion [m/Myr]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );

    target = (:sed, :ign, :volc, :plut)
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(figs)
        hᵢ = deepcopy(h)

        f = sim_cats[target[i]]
        c,m,e = binmeans(sim_ersn.age[f], sim_ersn.ersn[f], xmin, xmax, nbins)
        plot!(hᵢ, c, m, 
            yerror=2*e,
            label="",
            title="$(target[i])",
            color=colors[target[i]], lcolor=colors[target[i]], msc=:auto,
            markershape=:circle,
            linewidth=2,
        )
        figs[i] = hᵢ
    end

    h = plot(figs..., layout=(2,2), size=(1200,800))
    display(h)


## --- Model erosion as a function of age 
    # Data 
    f = sim_cats.sed;
    c,m,e = binmeans(sim_ersn.age[f], sim_ersn.ersn[f], xmin, xmax, nbins);
    plot(c, m, yerror=2e, markershape=:circle, color=:black, lcolor=:black, msc=:auto, label="")

    # σ = nanstd(log.(m))
    # μ = log.(nanmean(c))
    b = nanmean(m[22:38])

    d = LogNormal(nanmean(m), nanstd(m))
    μ=meanlogx(d)
    σ=stdlogx(d)

    d = fit(LogNormal, c)
    μ=meanlogx(d)
    σ=stdlogx(d)
    # LogNormal(μ, σ)

    # y = lognorm(c, [μ, σ, b])
    # plot!(c, y)

    # hline!([b])
    # vline!([σ, μ])

    # Start with a log-normal model, where 
        # p[1] = σ 
        # p[2] = μ 
        # p[3] = intercept / asymptote

    czoom = 0.01:0.01:4000
    @. lognorm(x, p) = p[1] / (x * p[3]) * exp(- (log(x) - p[2])^2 / (2*p[3]^2)) + p[4]
    p₀ = [100., 100., 10., 20.]
    plot(czoom, lognorm(czoom, p₀))

    fobj = curve_fit(lognorm, c, m, p₀)

    plot(c, lognorm(c, fobj.param))
    display(unique(lognorm(c, p₀)))


## ---
    p₀ = [1,0,nanmean(m[22:end])]
    x = 0.01:0.01:10
    y = lognorm(x, p₀)
    plot(x, y)
    vline!([1], label="mean")
    # vline!([0], label="stdev")


## --- End of file