## --- Set up 
    # Resample carbon isotope data and correct for post-depositional alteration 

    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots
    using Isoplot: yorkfit
    using LsqFit: curve_fit

    # Save figures to: 
    filepath = "results/figures/burial"

    # Definitions 
    nsims = Int(1e5)
    xmin, xmax, nbins = 0, 3800, 38
    isocolors = (;
        org_dark = :peru,
        org_light = :wheat,
        ct_dark = :black,
        carb_dark = :royalblue,
        carb_light = :lightblue,
    )

    # Load carbon isotope data 
    carbon = importdataset("data/carbonisotope/compilation.csv", ',', importas=:Tuple)


## --- Specific functions 
    # Function which describes Rayleigh curve
    @. r₀(HC, p) = p[1]/(HC/p[2] + p[3])^(p[4]-1) + p[5]

    # Find parameters to fit data to curve
    function fit_rayleigh(d13c, hc)
        p₀ = Float64[1, 2, 1e-3, 1.3, -30]
        lb = Float64[0, 1e-3, 1e-3, 1, -50]
        ub = Float64[Inf, 10, 1, 10, 0]

        fitted = curve_fit(r₀, hc, d13c, p₀, lower=lb, upper=ub)
        return fitted.param
    end

    # Calculate Rayleigh curve
    function rayleigh_curve(p, hc)
        x = 0:0.01:maximum(hc)
        y = r₀(x, p)
        return x, y
    end

    
## --- Model H/C as a function of age and assign H/C values to all organic carbon samples
    # Ages without uncertainty are assigned an uncertainty of 5%
    ageuncert = Array{Float64}(undef, length(carbon.age), 1)
    for i in eachindex(ageuncert)
        ageuncert[i] = ifelse(isnan(carbon.age_uncert[i]), carbon.age[i]*0.05, carbon.age_uncert[i])
    end

    # Fit model to real data; force a y intercept of 1.5
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc) & (carbon.hc > 0);
    t .&= carbon.std_fm_name .!= "Onverwacht Gp";

    x = [fill(0, 500); carbon.age[t]]
    x_sigma = [fill(1e-8, 500); ageuncert[t]]
    y = [log.(fill(1.5, 500)); log.(carbon.hc[t])]
    y_sigma = [fill(1e-8, 500); log.(fill(0.1, count(t)))]

    fobj = yorkfit(x, x_sigma, y, y_sigma)
    hc_age(age) = exp(age * (fobj.slope) + (fobj.intercept))

    # For any organic carbon value without a H/C ratio: 
        # Pick age randomly from a Gaussian distribution: mean = age, std = age uncert
        # Pick H/C randomly from a Gaussian distribution: mean = modeled HC, std = uncert 
    hc_assigned = Array{Float64}(undef, length(carbon.hc), 1)
    for i in eachindex(carbon.hc)
        hc = hc_age(randn() * ageuncert[i] + carbon.age[i])
        hc_assigned[i] = ifelse(!isnan(carbon.hc[i]), carbon.hc[i], randn() * hc.err + hc.val)
    end


## --- [PLOT] Model and plot H/C data 
    # Model 
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc) & (carbon.hc > 0);
    t .&= carbon.std_fm_name .!= "Onverwacht Gp";
    x = 50:50:3800
    y = Measurements.value.(hc_age.(x))
    e = Measurements.uncertainty.(hc_age.(x))
    h1 = plot(
        carbon.age[t], carbon.hc[t],
        label="Observed", 
        seriestype=:scatter, markersize=2,
        color=:black, msc=:auto, # alpha=0.5,
        xlabel="Age [Ma.]", ylabel="[LOG] H/C Ratio",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:bottomleft,
        yaxis=:log10
    )
    plot!(x, y, ribbon=2*e, 
        label="Modeled ± 2 s.d.", 
        linewidth=2, color=:goldenrod,
        fillalpha=0.25,
    )
    display(h1)
    savefig("$filepath/carbon_hc_model.pdf")

    # Modeled data 
    t = @. !isnan(carbon.d13c_org) 
    h2 = plot(hc_assigned[t], carbon.d13c_org[t], 
        seriestype=:scatter, label="Assigned", 
        color=:lightgrey, msc=:auto,
        markersize=1,
        framestyle=:box, 
        ylabel="d13c organic", xlabel="H/C ratio",
        legend=:bottomright,
        fg_color_legend=:white,
        fontfamily=:Helvetica,
    )
    plot!(carbon.hc[t], carbon.d13c_org[t], 
        seriestype=:scatter, label="Observed", 
        color=:black, msc=:auto,
        markersize=2
    )
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc)
    params = fit_rayleigh(carbon.d13c_org[t], carbon.hc[t])
    x,y = rayleigh_curve(params, carbon.hc[t])
    plot!(x,y, label="Rayleigh Model", color=:goldenrod, linewidth=3)
    display(h2)
    savefig("$filepath/carbon_hc_model_carbon.pdf")

    # Together -- temporary as a quick fix for thesis drafts
    title!(h1, "A. H/C Age Model", titleloc=:left)
    title!(h2, "B. H/C and δ13C", titleloc=:left)
    h = plot(h1, h2, layout=(2,1), size=(600,800), left_margin=(30,:px))
    display(h)
    savefig(h, "$filepath/carbon_hc_combined.pdf")


## --- Resample all data with spatiotemporal weights 
    # Isotope data without uncertainty is assigned an uncertainty of 0.02 ‰
    carbon.d13c_carb_uncert[isnan.(carbon.d13c_carb_uncert)] .= 0.02
    carbon.d13c_org_uncert[isnan.(carbon.d13c_org_uncert)] .= 0.02

    # Age data is assigned a minimum uncertainty of 10%
    for i in eachindex(carbon.age_uncert)
        carbon.age_uncert[i] = nanmax(carbon.age_uncert[i], carbon.age[i]*0.05)
    end

    # Inorganic carbon isotope data 
    t = .!isnan.(carbon.d13c_carb)
    k = invweight(carbon.lat[t], carbon.lon[t], carbon.age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    c,m,el,eu = bin_bsr_means(carbon.age[t], carbon.d13c_carb[t], xmin, xmax, nbins,
        x_sigma = carbon.age_uncert[t],
        y_sigma = carbon.d13c_carb_uncert[t],
        nresamplings = nsims,
        p = p
    )
    sim_carb = (c = c, m = m, el = el, eu = eu,)

    # Organic carbon isotope data 
    t = .!isnan.(carbon.d13c_org)
    k = invweight(carbon.lat[t], carbon.lon[t], carbon.age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    c,m,el,eu = bin_bsr_means(carbon.age[t], carbon.d13c_org[t], xmin, xmax, nbins,
        x_sigma = carbon.age_uncert[t],
        y_sigma = carbon.d13c_org_uncert[t],
        nresamplings = nsims,
        p = p
    )
    sim_org = (c = c, m = m, el = el, eu = eu,)

    
## --- Correct organic carbon for post-depositional alteration 
    # Model with existing data 
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc)
    params = fit_rayleigh(carbon.d13c_org[t], carbon.hc[t])

    # Correct observed data
    corrected_obs = carbon.d13c_org[t] .- (vec(r₀(carbon.hc[t], params)) .-  r₀(1.5, params))
    corrected_mod = carbon.d13c_org .- (vec(r₀(hc_assigned, params)) .-  r₀(1.5, params))

    # Correcte and resample observed data 
    t = .!isnan.(carbon.d13c_org)
    k = invweight(carbon.lat[t], carbon.lon[t], carbon.age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    c,m,el,eu = bin_bsr_means(carbon.age[t], corrected_mod[t], xmin, xmax, nbins,
        x_sigma = carbon.age_uncert[t],
        y_sigma = carbon.d13c_org_uncert[t],
        nresamplings = nsims,
        p = p
    )
    sim_corr = (c = c, m = m, el = el, eu = eu,)


## --- Save resampled and corrected data to a file 
    fid = h5open("src/burial/resampled_carbon.h5", "w")
    g = create_group(fid, "vars")
    g_carb = create_group(g, "carb")
        for k in keys(sim_carb)
            g_carb["$k"] = collect(sim_carb[k])
        end 
    g_org = create_group(g, "org")
        for k in keys(sim_org)
            g_org["$k"] = collect(sim_org[k])
        end 
    g_corr = create_group(g, "corrected")
        for k in keys(sim_org)
            g_corr["$k"] = collect(sim_corr[k])
        end 
    close(fid)


## --- [PLOT] Carbon isotope record over geologic time 
    h = plot(
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="δ13C [‰]",
        fg_color_legend=:white,
        legend=:bottomleft,
        size=(600,1000),
        ylims=(-50, 20),
        left_margin=(30,:px),
        fontfamily=:Helvetica,
        labelfontsize=18, tickfontsize=14,legendfontsize=14,
    )

    # Plot a random selection of observed data 
    indices = 1:length(carbon.age)
    t = rand(indices[.!isnan.(carbon.d13c_carb)], 5_000)
    plot!(carbon.age[t], carbon.d13c_carb[t], 
        label="",
        color=isocolors.carb_light, msc=:auto, 
        seriestype=:scatter,
        markersize=1
    )
    # t = rand(indices[.!isnan.(carbon.d13c_org)], 5_000)
    t = trues(length(carbon.age))
    plot!(carbon.age[t], carbon.d13c_org[t], 
        label="",
        color=isocolors.org_light, msc=:auto, 
        seriestype=:scatter,
        markersize=1
    )

    # Resampled means
    hline!([0], label="",
        color=isocolors.carb_dark,
        linestyle=:dash,
    )
    hline!([-25], label="",
        color=isocolors.ct_dark,
        linestyle=:dash,
    )
    plot!(sim_carb.c, sim_carb.m, 
        yerror=(2*sim_carb.el, 2*sim_carb.eu),  
        label="Carbonate", 
        color=isocolors.carb_dark, lcolor=isocolors.carb_dark, msc=:auto, 
        markershape=:circle,
        seriestype=:scatter,
    )
    plot!(sim_org.c, sim_org.m, 
        yerror=(2*sim_org.el, 2*sim_org.eu), 
        label="Organic [Observed]", 
        color=isocolors.org_dark, lcolor=isocolors.org_dark, msc=:auto, 
        markershape=:circle,
        seriestype=:scatter,
    )
    plot!(sim_corr.c, sim_corr.m, 
        yerror=(2*sim_corr.el, 2*sim_corr.eu), 
        label="Organic [Corrected]", 
        color=isocolors.ct_dark, lcolor=isocolors.ct_dark, msc=:auto, 
        markershape=:circle,
        seriestype=:scatter,
    )
    display(h)
    savefig(h, "$filepath/carbon_isotope_record.pdf")
    

## --- [PLOT] To correct or not to correct, and the consequences thereof 
    h = plot(
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="δ13C [‰]",
        fg_color_legend=:white,
        legend=:bottomleft,
        size=(600,600),
        ylims=(-60,-12),
        left_margin=(30,:px),
        fontfamily=:Helvetica,
    )
    c,m,e = binmeans(carbon.age, carbon.d13c_org, xmin, xmax, nbins)
    plot!(c[.!isnan.(m)], m[.!isnan.(m)], yerror=2e, 
        label="Observed", 
        color=isocolors.org_dark, lcolor=isocolors.org_dark, msc=isocolors.org_dark, 
        markershape=:circle,
        # linestyle=:dash,
    )

    # Modeled data 
    h1 = deepcopy(h)
    t = @. !isnan(carbon.d13c_org)
    c,m,e = binmeans(carbon.age[t], corrected_mod[t], xmin, xmax, nbins)
    plot!(h1, c[.!isnan.(m)], m[.!isnan.(m)], yerror=2e, 
        label="Corrected [Modeled H/C]", 
        # label="",
        color=isocolors.ct_dark, lcolor=isocolors.ct_dark, msc=:auto, 
        markershape=:circle,
        # linestyle=:dash,
    )
    display(h1)
    savefig(h1, "$filepath/correction_modeled.pdf")
    
    # Observed Data
    h2 = deepcopy(h)
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc)
    c,m,e = binmeans(carbon.age[t], corrected_obs, xmin, xmax, nbins)
    plot!(h2, c[.!isnan.(m)], m[.!isnan.(m)], yerror=2e, 
        label="Corrected [Observed H/C]", 
        # label="",
        color=isocolors.ct_dark, lcolor=isocolors.ct_dark, msc=:auto, 
        markershape=:circle,
        # linestyle=:dash,
    )
    display(h2)
    savefig(h2, "$filepath/correction_observed.pdf")


## --- [PLOT] Fraction of carbon buried as organic 
    # Definitions  
    mantle = -5
    carb = sim_carb.m .± nanmean([sim_carb.el sim_carb.eu], dims=2)
    org = sim_org.m .± nanmean([sim_org.el sim_org.eu], dims=2)
    org_corrected = sim_corr.m .± nanmean([sim_corr.el sim_corr.eu], dims=2)

    h = plot(
        xlabel="Age [Ma.]", ylabel="Fraction Buried as Organic",
        framestyle=:box,
        fontfamily=:Helvetica,
        ylims=(0.05,0.3),
        size=(400,500),
        legend=:topright,
        fg_color_legend=:white
    );

    # Significant ages 
    vline!([2300, 717, 541, 252], label="", 
        linestyle=:dot,
        color=:black,
        alpha=0.5
    )
    annotate!([(2300,0.055, text("GOE", 8, :left, :bottom, :black, rotation=90))])
    annotate!([(717,0.055, text("Snowball Earth Initiation", 8, :left, :top, :black, rotation=90))])
    annotate!([(541,0.055, text("Cambrian Explosion", 8, :left, :bottom, :black, rotation=90))])
    annotate!([(252,0.055, text("end-Permian Extinction", 8, :left, :bottom, :black, rotation=90))])

    # Des Marais curve 
    des_marais = (;
        age=[2650.0, 2495.5516180173463, 2047.2862072181294, 1949.6316329385436, 
            1853.1940688240234, 1747.3141844633037, 1646.8618856663252, 1553.2220460691974, 
            1451.8744754266536, 1350.582859274457, 1251.9162470079896, 1051.5664666811736, 
            957.5075512657118, 850.7471608277119, 756.0325287284863, 656.6550224336061],
        forg=[0.08958593677142582, 0.10020889676396522, 0.1494628368926606, 0.19049706238925668, 
            0.17004010071808257, 0.11977338431409115, 0.12975286766763028, 0.14035064813951315, 
            0.14039261400727404, 0.14105567471789607, 0.16009238967121547, 0.1497514947102282, 
            0.19076697804304346, 0.2210044219590288, 0.21133867232428727, 0.14991501911778413]
    )
    plot!(des_marais.age, des_marais.forg, label="Des Marais et al.",
        markershape=:circle, linestyle=:dash,
        color=:white, lcolor=:black, msc=:black,
    )

    # Measured carbonate record
    h1 = deepcopy(h)
    frog = (mantle .- carb) ./ (org_corrected .- carb)
    frog = (;
        val = Measurements.value.(frog),
        err = Measurements.uncertainty.(frog),
    )
    plot!(h1, c, frog.val,
        ribbon = 2*frog.err,
        label="This Study",
        color=:seagreen, lcolor=:seagreen, msc=:auto, 
        markershape=:circle,
        markersize=3,
        # seriestype=:scatter,
    )
    display(h1)
    savefig(h1, "$filepath/f_org.pdf")

    # 0‰ carbon record 
    h2 = deepcopy(h)
    carb = 0.0 ± 0.5
    frog = (mantle .- carb) ./ (org_corrected .- carb)
    frog = (;
        val = Measurements.value.(frog),
        err = Measurements.uncertainty.(frog),
    )
    plot!(h2, c, frog.val,
        ribbon = 2*frog.err,
        label="This Study",
        color=:seagreen, lcolor=:seagreen, msc=:auto, 
        markershape=:circle,
        markersize=3,
        # seriestype=:scatter,
    )
    display(h2)
    savefig(h2, "$filepath/f_org_carbvariant.pdf")


## --- End of file 