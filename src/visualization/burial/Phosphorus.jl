## --- Set up 
    # Potential correlation between the ratio of phosphorus to alkalininty and 
    # carbon burial 
    
    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots, Colors

    # Save figures to: 
    filepath = "results/figures/burial"

    # Script-wide definitions 
    nsims = Int(1e6)
    xmin, xmax, nbins = 0, 3800, 38
    age_error = 0.05                   # Minimum age error (%)
    age_error_abs = 50                 # Minimum age error (Ma)


## --- Load data 
    # Carbon isotope data
    carbon = importdataset("data/carbonisotope/compilation.csv", ',', importas=:Tuple)

    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

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


## --- Temporal resample P/Alk ratio 
    # L + ratio + we can't resample everything and then take a ratio, because that's bad 
    # for some reason. So calculate the ratio [actually, calculate X/(X+Y)], and then 
    # resample that before re-converting back into ratio. It's true! Google Keller and Schoene 
    # 2012 extended methods for more info

    # Definitions
    P2O5_err = 2.0      # Error set as 10x the error given in 
    alk_err = 1.0       # volcanic.mat (Keller et al., 2015)   

    # Get sample ages and uncertainties
    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
        uncert_rel=age_error, uncert_abs=age_error_abs
    )

    # Preallocate 
    alkalinity = Array{Float64}(undef, length(mbulk.SiO2), 1)
    phosphorus = Array{Float64}(undef, length(mbulk.SiO2), 1)

    # target = (:sed, :volc, :plut)
    target = (:plut, :volc, :sed)
    simout = NamedTuple{target}(Array{Float64}(undef, nbins, 4) for _ in target)

    # Conversion factors from wt.% element oxide to moles per 100g sample
    CaO_to_Ca =   (molarmass["Ca"]   + molarmass["O"]  )
    MgO_to_Mg =   (molarmass["Mg"]   + molarmass["O"]  )
    K2O_to_K =    (molarmass["K"] *2 + molarmass["O"]  ) * 2   # 2 mol K / 1 mol K₂O
    Na2O_to_Na =  (molarmass["Na"]*2 + molarmass["O"]  ) * 2   # 2 mol Na / 1 mol Na₂O
    P2O5_to_mol = (molarmass["P"] *2 + molarmass["O"]*5)

    # Moles of alkalinity (charge), phosphorus, 
    for i in eachindex(alkalinity)
        Ca²⁺ = mbulk.CaO[i] * CaO_to_Ca * 2     # +2 change
        Mg²⁺ = mbulk.MgO[i] * MgO_to_Mg * 2     # +2 change
        K⁺ = mbulk.K2O[i] * K2O_to_K
        Na⁺ = mbulk.Na2O[i] * Na2O_to_Na
        alkalinity[i] = Ca²⁺ + Mg²⁺ + K⁺ + Na⁺

        phosphorus[i] = mbulk.P2O5[i] * P2O5_to_mol
    end
    alkalinity = vec(alkalinity)
    phosphorus = vec(phosphorus)

    # Calculate temporal weights and resample ratio
    for key in target 
        t = @. match_cats[key] & !isnan(phosphorus) & !isnan(alkalinity)
        t .&= .!(match_cats.phosphorite) .| (mbulk.P2O5 .< 4); 

        k = invweight_age(sampleage[t])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        c,m,el,eu = bin_bsr_ratios(sampleage[t], phosphorus[t], alkalinity[t], 
            xmin, xmax, nbins,
            x_sigma = ageuncert[t],
            num_sigma = fill(P2O5_err, count(t)),
            denom_sigma = fill(alk_err, count(t)),
            p = p
        )
        simout[key] .= [c m el eu]
    end

    # Sanity check plot 
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(simout))
    c,m,el,eu = 1,2,3,4;
    h = plot(
        # xlabel="Age [Ma.]", ylabel="P / Alk [mol. ratio]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topleft,
    );
    for i in eachindex(target)
        key = target[i]
        h1 = deepcopy(h)
        plot!(h1, simout[key][:,c], simout[key][:,m], 
            yerror=(2*simout[key][:,el], 2*simout[key][:,eu]), 
            label="$key", 
            color=colors[key], lcolor=colors[key], msc=:auto, 
            markershape=:circle,
        )
        figs[i] = h1
    end
    
    h2 = plot(figs..., layout=(1,3), size=(1800,400))
    display(h2)
    savefig(h2, "$filepath/p_alk.pdf")


## --- Plot nonresampled data because what the fuck is going on 
    ratio = phosphorus ./ alkalinity;
    t = match_cats.sed .& (mbulk.P2O5 .< 1)
    c,m,e = binmeans(sampleage[t], ratio[t], 0,3800,38);
    plot(c, m, yerror=2e, markershape=:circle)

    s = (3100 .< sampleage .< 3300) .& match_cats.sed;
    count(s)
    plot(mbulk.P2O5[match_cats.sed .& s], alkalinity[match_cats.sed .& s], label="",
        xlabel="P wt.%", ylabel="Alk mols",
        seriestype=:scatter, 
    )
    plot!([0, 1],[0, P2O5_to_mol/0.1], label="0.1")
    plot!([0, 1],[0, P2O5_to_mol/0.04], label="0.04")
    plot!([0, 1],[0, P2O5_to_mol/0.02], label="0.02")

    # histogram(alkalinity[match_cats.sed .& s])


## --- variation alk vs p over time??
    h = plot(
        framestyle=:box,
    )
    c,m,e = binmeans(sampleage[match_cats.sed], phosphorus[match_cats.sed], 0,3800,38);
    plot!(c, m, label="", color=:darkorange, markershape=:circle, ylabel="P")

    c,m,e = binmeans(sampleage[match_cats.sed], alkalinity[match_cats.sed], 0,3800,38);
    plot!(twinx(), c, m, label="", color=:blue, markershape=:circle, ylabel="Alk")


## --- Phosphorus abundance over time: why are things weird? 
    h = plot(
        xlabel="Age [Ma.]", ylabel="P [wt.%]",
        framestyle=:box,
        fontfamily=:Helvetica,
        legend=:topright
    )
    c,m,e = binmeans(mbulk.Age[match_cats.sed], mbulk.P2O5[match_cats.sed], 0,3800,38)
    plot!(h, c, m, yerror=2e, 
        label="All Seds",
        color=colors.sed, lcolor=colors.sed, msc=:auto, markershape=:circle
    )
    t = match_cats.sed .& .!match_cats.phosphorite .& (mbulk.P2O5 .< 3)
    c,m,e = binmeans(mbulk.Age[t], mbulk.P2O5[t], 0,3800,38)
    plot!(h, c, m, yerror=2e, 
        label="Seds w/o Phosphorites",
        color=colors.carb, lcolor=colors.carb, msc=:auto, markershape=:circle
    )
    display(h)
    savefig(h, "$filepath/p_unresampled.pdf")

    # The cutoff for excluding phosphorites is super sensitive (2 vs. 4), which makes me 
    # feel like there's perhaps some outliers in there 
    t = 450 .< sampleage .< 650;
    t .&= match_cats.sed .& .!match_cats.phosphorite
    c,n = bincounts(mbulk.P2O5[t], 0, 6, 6)
    h = plot(c,n,label="", 
        seriestype=:bar, 
        yaxis=:log10,
        xlabel="P [wt.%]", ylabel="Count",
        framestyle=:box,
    )


## --- Spatiotemporal resample d13C data 
    # # Assign age (5% or 50 Ma.) and isotope (0.02 ‰) uncertainties to samples without them 
    # carb_uncert = Array{Float64}(undef, length(carbon.d13c_carb), 1)
    # org_uncert = Array{Float64}(undef, length(carbon.d13c_carb), 1)
    # ageuncert = Array{Float64}(undef, length(carbon.d13c_carb), 1)
    
    # for i in eachindex(carb_uncert)
    #     carb_uncert[i] = ifelse(isnan(carbon.d13c_carb_uncert[i]), 0.02, carbon.d13c_carb_uncert[i])
    #     org_uncert[i] = ifelse(isnan(carbon.d13c_org_uncert[i]), 0.02, carbon.d13c_org_uncert[i])
    #     ageuncert[i] = nanmaximum([carbon.age_uncert[i], carbon.age[i]*age_error, age_error_abs])
    # end
    
    # # Carbonate 
    # t = .!isnan.(carbon.d13c_carb)
    # k = invweight(carbon.lat[t], carbon.lon[t], carbon.age[t])
    # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    # sim_carb = bsresample(
    #     [carbon.d13c_carb[t] carbon.age[t]], 
    #     [carbon.d13c_carb_uncert[t] carbon.age_uncert[t]], 
    #     nsims, p
    # )

    # # Organic 
    # t = .!isnan.(carbon.d13c_org)
    # k = invweight(carbon.lat[t], carbon.lon[t], carbon.age[t])
    # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    # sim_org = bsresample(
    #     [carbon.d13c_org[t] carbon.age[t]], 
    #     [carbon.d13c_org_uncert[t] carbon.age_uncert[t]], 
    #     nsims, p
    # )

    # # Sanity check plot 
    # h = plot(
    #     xlabel="Age [Ma.]", ylabel="δ13C",
    #     framestyle=:box,
    #     fontfamily=:Helvetica
    # )
    # c,m,e = binmeans(sim_carb[:,2], sim_carb[:,1], xmin, xmax, nbins)   # Doesn't smooth CIE
    # plot!(c[.!isnan.(m)],m[.!isnan.(m)], yerror=2e[.!isnan.(m)], markershape=:circle, 
    #     label="Carbonate", 
    #     color=:blue, lcolor=:blue, msc=:auto,
    # )

    # c,m,e = binmeans(sim_org[:,2], sim_org[:,1], xmin, xmax, nbins) 
    # plot!(c[.!isnan.(m)],m[.!isnan.(m)], yerror=2e[.!isnan.(m)], markershape=:circle, 
    #     label="Organic",
    #     color=:darkorange, lcolor=:darkorange, msc=:auto,
    # )

    # display(h)
    # savefig(h, "$filepath/carbonisotope_binned.pdf")


## --- TO DO: 
    # Assign H/C ratios for correction to all samples 
    # Correct samples for diagenesis 
    # d13c_org_corr = data.d13c_org .- (r₀(hc_assigned, params) .-  r₀(1.5, params))

    # Create F(org) curve and plot next to P/Alk ratio     



## --- Plot 
    # p = Plots.palette(colorpalette, 5)
    # p2 = Plots.palette(:isoluminant_cm_70_c39_n256, 3, rev=true)
    
    # fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target));
    # figname = ("Whole Earth", "Sedimentary", "Igneous", "Volcanic", "Plutonic")
    # for i in eachindex(target)
    #     k = target[i] 

    #     # Initialize plot
    #     h = Plots.plot(
    #         framestyle=:box,
    #         fontfamily=:Helvetica, 
    #         # xlabel="Age [Ma.]", ylabel="P / Alk [mol. ratio]",
    #         # fg_color_legend=:white, legend=:topright, legendfontsize=14,
    #         labelfontsize=labelfontsize,
    #         legend=false,
    #         xlims=(-50, 3850)
    #     )

    #     # Get (future) ylims 
    #     lim_l = round(minimum(simout_m[k] .- simout_e[k][:,1])*0.9, sigdigits=4)
    #     lim_u = round(maximum(simout_m[k] .+ simout_e[k][:,2])*1.1, sigdigits=4)

    #     # Color GTS eons 
    #     Plots.plot!(h, [2500, 3850, 3850, 2500], [lim_l, lim_l, lim_u, lim_u],
    #         seriestype=:shape, color=p2[1], alpha=0.15, lcolor=:match
    #     )
    #     Plots.plot!(h, [541, 2500, 2500, 541], [lim_l, lim_l, lim_u, lim_u],
    #         seriestype=:shape, color=p2[2], alpha=0.15, lcolor=:match
    #     )
    #     Plots.plot!(h, [-50, 541, 541, -50], [lim_l, lim_l, lim_u, lim_u],
    #         seriestype=:shape, color=p2[3], alpha=0.15, lcolor=:match
    #     )

    #     # Plot data and reset y limits
    #     Plots.plot!(h, c, simout_m[k], yerror=(simout_e[k][:,1], simout_e[k][:,2],), 
    #         left_margin=(15,:px), bottom_margin=(15,:px), 
    #         label="$(figname[i])", color=p[i], lcolor=p[i], msc=:auto,
    #         seriestype=:scatter)
    #     Plots.ylims!(lim_l, lim_u)
    #     Plots.annotate!(((0.9, 0.97), (figname[i] * "\nn = $(count(class[k]))", 
    #         labelfontsize, :right, :top)))
    #     fig[i] = h
    # end

    # # Assemble plots
    # xlabel!(fig[end], "Age [Ma.]")
    # ylabel!(fig[1], "P / Alk [mol. ratio]")
    # ylabel!(fig[end-1], "P / Alk [mol. ratio]")
    # h = Plots.plot(fig..., layout=(2, 3), size=(1600, 1000), 
    #     left_margin=(30,:px), bottom_margin=(30,:px)
    # )


## --- End of file 