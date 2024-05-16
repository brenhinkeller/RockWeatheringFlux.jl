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
    nsims = Int(1e6)
    xmin, xmax, nbins = 0, 3800, 38
    isocolors = (;
        org_light = :navajowhite,
        org_dark = :darkorange,
        ct_dark = :seagreen,
        ct_light = :palegreen,
        carb_light = :lightblue,
        carb_dark = :royalblue,
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


## --- [PLOT] Model and modeled H/C data 
    # Model 
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc) & (carbon.hc > 0);
    t .&= carbon.std_fm_name .!= "Onverwacht Gp";
    x = 50:50:3800
    y = Measurements.value.(hc_age.(x))
    e = Measurements.uncertainty.(hc_age.(x))
    h = plot(
        carbon.age[t], carbon.hc[t],
        label="Observed", 
        seriestype=:scatter, markersize=3,
        color=:darkturquoise, msc=:auto,
        xlabel="Age [Ma.]", ylabel="[LOG] H/C Ratio",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:bottomleft,
        yaxis=:log10
    )
    plot!(x, y, ribbon=2*e, 
        label="Modeled ± 2 s.d.", 
        linewidth=2, color=:teal
    )
    display(h)
    savefig("$filepath/carbon_hc_model.pdf")

    # Modeled data 
    t = @. !isnan(carbon.d13c_org) 
    h = plot(hc_assigned[t], carbon.d13c_org[t], 
        seriestype=:scatter, label="Assigned", 
        color=:lightblue, msc=:auto,
        markersize=2,
        framestyle=:box, 
        ylabel="d13c organic", xlabel="H/C ratio",
        legend=:topright,
        fg_color_legend=:white,
    )
    plot!(carbon.hc[t], carbon.d13c_org[t], 
        seriestype=:scatter, label="Observed", 
        color=:red, msc=:auto,
        markersize=2
    )
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc)
    params = fit_rayleigh(carbon.d13c_org[t], carbon.hc[t])
    x,y = rayleigh_curve(params, carbon.hc[t])
    plot!(x,y, label="Rayleigh Model", color=:black, linewidth=2)
    display(h)
    savefig("$filepath/carbon_hc_model_carbon.pdf")


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
    resampled = bsresample([carbon.d13c_carb[t] carbon.age[t]],
        [carbon.d13c_carb_uncert[t] carbon.age_uncert[t]],
        nsims,p
    )
    sim_carb = (;
        d13c_carb = resampled[:,1],
        age = resampled[:,2]
    )

    # Organic carbon isotope data 
    t = .!isnan.(carbon.d13c_org)
    k = invweight(carbon.lat[t], carbon.lon[t], carbon.age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    resampled = bsresample([carbon.d13c_org[t] hc_assigned[t] carbon.age[t]],
        [carbon.d13c_org_uncert[t] fill(0.01, count(t)) carbon.age_uncert[t]],
        nsims,p
    )
    sim_org = (;
        d13c_org = resampled[:,1],
        hc = resampled[:,2],
        age = resampled[:,3],
    )


## --- Correct organic carbon for post-depositional alteration 
    # Model with existing data 
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc)
    params = fit_rayleigh(carbon.d13c_org[t], carbon.hc[t])

    # Correct observed and resampled data
    corrected_min = carbon.d13c_org[t] .- (vec(r₀(carbon.hc[t], params)) .-  r₀(1.5, params))
    corrected_obs = carbon.d13c_org .- (vec(r₀(hc_assigned, params)) .-  r₀(1.5, params))
    corrected_sim = sim_org.d13c_org .- (vec(r₀(sim_org.hc, params)) .-  r₀(1.5, params))


## --- Save resampled and corrected data to a file 
    fid = h5open("src/visualization/burial/resampled_carbon.h5", "w")
    g = create_group(fid, "vars")
    g_carb = create_group(g, "carb")
        for k in keys(sim_carb)
            g_carb["$k"] = sim_carb[k]
        end 
    g_org = create_group(g, "org")
        for k in keys(sim_org)
            g_org["$k"] = sim_org[k]
        end 
        g_org["d13c_org_corrected"] = corrected_sim
    close(fid)


## --- [PLOT] Impact of correcting post-depositional alteration on isotope record
    # Resampled data
    h_sim = plot(
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="d13c",
        fg_color_legend=:white,
        legend=:bottomleft,
        size=(600,1200),
        ylims=(-50, 5),
        title="Spatiotemporal Resampled", titleloc=:left,
    )
    c,m,e = binmeans(sim_org.age, sim_org.d13c_org, xmin, xmax, nbins)
    plot!(c, m, yerror=2e, 
        label="Observed", 
        color=isocolors.org_dark, lcolor=isocolors.org_dark, msc=:auto, 
        markershape=:circle,
    )
    c,m,e = binmeans(sim_org.age, corrected_sim, xmin, xmax, nbins)
    plot!(c, m, yerror=2e, 
        label="Corrected", 
        color=isocolors.ct_dark, lcolor=isocolors.ct_dark, msc=:auto, 
        markershape=:circle,
    )
    c,m,e = binmeans(sim_carb.age, sim_carb.d13c_carb, xmin, xmax, nbins, relbinwidth=2)
    plot!(c, m, yerror=2e, 
        label="Carbonate [200 Ma. avg.]", 
        color=isocolors.carb_dark, lcolor=isocolors.carb_dark, msc=:auto, 
        markershape=:circle,
    )
    display(h_sim)

    # Observed data, no resampling
    h_obs = plot(
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="d13c",
        fg_color_legend=:white,
        legend=:bottomleft,
        size=(600,1200),
        ylims=(-50, 5),
        title="Observed", titleloc=:left,
    )
    c,m,e = binmeans(carbon.age, carbon.d13c_org, xmin, xmax, nbins)
    plot!(c[.!isnan.(m)], m[.!isnan.(m)], yerror=2e, 
        label="Observed", 
        color=isocolors.org_dark, lcolor=isocolors.org_dark, msc=:auto, 
        markershape=:circle,
    )
    c,m,e = binmeans(carbon.age, corrected_obs, xmin, xmax, nbins)
    plot!(c[.!isnan.(m)], m[.!isnan.(m)], yerror=2e, 
        label="Corrected", 
        color=isocolors.ct_dark, lcolor=isocolors.ct_dark, msc=:auto, 
        markershape=:circle,
    )
    t = @. !isnan(carbon.d13c_org) & !isnan(carbon.hc)
    c,m,e = binmeans(carbon.age[t], corrected_min, xmin, xmax, nbins)
    plot!(c[.!isnan.(m)], m[.!isnan.(m)], yerror=2e, 
        label="Corrected [Obs. H/C]", 
        color=:purple, lcolor=:purple, msc=:auto, 
        markershape=:circle,
    )
    c,m,e = binmeans(carbon.age, carbon.d13c_carb, xmin, xmax, nbins, relbinwidth=2)
    plot!(c[.!isnan.(m)], m[.!isnan.(m)], yerror=2e, 
        label="Carbonate [200 Ma. avg.]", 
        color=isocolors.carb_dark, lcolor=isocolors.carb_dark, msc=:auto, 
        markershape=:circle,
    )
    display(h_obs)

    # Together!
    h = plot(h_sim, h_obs, layout=(1, 2), size=(1200,1200))
    display(h)
    savefig("$filepath/carbon_isotope_record.pdf")


## --- [PLOT] Inorganic carbonate record 
    h = plot(
        framestyle=:box,
        xlabel="Age [Ma.]", ylabel="d13c",
        fg_color_legend=:white,
    )
    # t = rand(1:length(carbon.age), 5_000)     # Random selection of data
    t = trues(length(carbon.age))
    plot!(carbon.age[t], carbon.d13c_carb[t], 
        label="Observed Record",
        color=isocolors.carb_light, msc=:auto, 
        seriestype=:scatter,
        markersize=1
    )
    c,m,e = binmeans(sim_carb.age, sim_carb.d13c_carb, xmin, xmax, nbins)
    plot!(c, m, yerror=2e, 
        label="Resampled Means", 
        color=isocolors.carb_dark, lcolor=isocolors.carb_dark, msc=:auto, 
        markershape=:circle,
    )
    display(h)
    savefig("$filepath/carbon_inorganic_record.pdf")


## --- [PLOT] Fraction of carbon buried as organic 
    # Define mantle and carbonate values 
    mantle = -5.5
    c,carbonate,e = binmeans(sim_carb.age, sim_carb.d13c_carb, xmin, xmax, nbins, relbinwidth=2)

    h = plot(
        xlabel="Age [Ma.]", ylabel="Fraction Buried as Organic",
        framestyle=:box,
        fontfamily=:Helvetica,
        ylims=(0,0.3),
        size=(400,500),
        legend=:bottomleft,
        fg_color_legend=:white
    );

    # Resampled uncorrected
    c,m,e = binmeans(sim_org.age, sim_org.d13c_org, xmin, xmax, nbins, relbinwidth=1)
    frog = (mantle .- carbonate) ./ (m .- carbonate)
    plot!(c, frog, label="Uncorrected",
        color=isocolors.org_dark, lcolor=isocolors.org_dark, msc=:auto, 
        markershape=:circle,
    )

    # Resampled corrected
    c,m,e = binmeans(sim_org.age, corrected_sim, xmin, xmax, nbins, relbinwidth=1)
    frog = (mantle .- carbonate) ./ (m .- carbonate)
    plot!(c, frog, label="Corrected",
        color=isocolors.ct_dark, lcolor=isocolors.ct_dark, msc=:auto, 
        markershape=:circle,
    )

    # # Observed uncorrected
    # c,m,e = binmeans(carbon.age, carbon.d13c_org, xmin, xmax, nbins, relbinwidth=1)
    # frog = (mantle .- carbonate) ./ (m .- carbonate)
    # plot!(c, frog, label="Uncorrected",
    #     color=isocolors.org_dark, lcolor=isocolors.org_dark, msc=:auto, 
    #     markershape=:utriangle,
    # )

    # # Observed corrected
    # c,m,e = binmeans(carbon.age, corrected_obs, xmin, xmax, nbins, relbinwidth=1)
    # frog = (mantle .- carbonate) ./ (m .- carbonate)
    # plot!(c, frog, label="Corrected",
    #     color=isocolors.ct_dark, lcolor=isocolors.ct_dark, msc=:auto, 
    #     markershape=:utriangle,
    # )

    # # Des Marais curve 
    # des_marais = (;
    #     age=[2650.0, 2495.5516180173463, 2047.2862072181294, 1949.6316329385436, 
    #         1853.1940688240234, 1747.3141844633037, 1646.8618856663252, 1553.2220460691974, 
    #         1451.8744754266536, 1350.582859274457, 1251.9162470079896, 1051.5664666811736, 
    #         957.5075512657118, 850.7471608277119, 756.0325287284863, 656.6550224336061],
    #     forg=[0.08958593677142582, 0.10020889676396522, 0.1494628368926606, 0.19049706238925668, 
    #         0.17004010071808257, 0.11977338431409115, 0.12975286766763028, 0.14035064813951315, 
    #         0.14039261400727404, 0.14105567471789607, 0.16009238967121547, 0.1497514947102282, 
    #         0.19076697804304346, 0.2210044219590288, 0.21133867232428727, 0.14991501911778413]
    # )
    # plot!(des_marais.age, des_marais.forg, label="Des Marais",
    #     markershape=:circle, linestyle=:dash,
    #     color=:black, msc=:auto,
    # )

    display(h)
    savefig("$filepath/f_org.pdf")


## --- End of file 