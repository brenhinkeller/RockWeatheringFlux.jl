## --- Set up 
    # Changes in the ratio of phosphorus to alkalininty over geologic time 

    # Load data and base packages
    include("Definitions.jl");

    # Definitions
    age_error = 0.05                         # Minimum age error (%)
    age_error_abs = 50                       # Minimum age error (Ma)
    P2O5_err = 2.0                           # Error set as 10x the error given in 
    alk_err = 1.0                            # volcanic.mat (Keller et al., 2015)   


## --- Resampling setup: calculate alkalinity and phosphorus for every sample
    # L + ratio + we can't resample everything and then take a ratio, because that's bad 
    # for some reason. So calculate the ratio [actually, calculate X/(X+Y)], and then 
    # resample that before re-converting back into ratio. It's true! Google Keller et al., 
    # 2012 extended methods for more info

    # Conversion factor from wt.% element oxide to moles per 100g sample
    CaO_to_Ca =   (molarmass["Ca"]   + molarmass["O"]  )
    MgO_to_Mg =   (molarmass["Mg"]   + molarmass["O"]  )
    K2O_to_K =    (molarmass["K"] *2 + molarmass["O"]  ) * 2   # 2 mol K / 1 mol K₂O
    Na2O_to_Na =  (molarmass["Na"]*2 + molarmass["O"]  ) * 2   # 2 mol Na / 1 mol Na₂O
    P2O5_to_mol = (molarmass["P"] *2 + molarmass["O"]*5)

    # Calculate moles of alkalinity (charge) for every sample
    alk = Array{Float64}(undef, length(mbulk.SiO2), 1)
    for i in eachindex(alk)
        Ca²⁺ = mbulk.CaO[i] * CaO_to_Ca * 2     # +2 change
        Mg²⁺ = mbulk.MgO[i] * MgO_to_Mg * 2
        K⁺ = mbulk.K2O[i] * K2O_to_K
        Na⁺ = mbulk.Na2O[i] * Na2O_to_Na
        alk[i] = Ca²⁺ + Mg²⁺ + K⁺ + Na⁺
    end

    # Calculate moles of phosphorus for every sample 
    phosphorus = mbulk.P2O5 .* P2O5_to_mol

    # Calculate sample age
    # Use the sample age, unless the sample doesn't have an age: then use map age
    t = @. !isnan(mbulk.Age);
    sampleage = copy(mbulk.Age);
    ageuncert = nanadd.(mbulk.Age_Max, .- mbulk.Age_Min) ./ 2;
    sampleage[t] .= macrostrat.age[t]
    ageuncert[t] .= nanadd.(macrostrat.agemax[t], .- macrostrat.agemin[t]) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = max(sampleage[i]*age_error, ageuncert[i], age_error_abs)
    end


## --- Resample, temporal weights 
    # Look at bulk rock, sedimentary, igneous, and separate volcanic / plutonic 
    # May be worth it to separate basalts from granites, since we don't really want to 
    # get into weird shit with fractional crystallization and remelting 

    # Preallocate 
    xmin, xmax, nbins = 0, 3800, 38
    c = cntr(xmin:(xmax-xmin)/nbins:xmax)
    
    target = (:bulk, :sed, :ign, :volc, :plut)
    simout_m = NamedTuple{target}(Array{Float64}(undef, nbins) for _ in eachindex(target))
    simout_e = NamedTuple{target}(Array{Float64}(undef, nbins, 2) for _ in eachindex(target))

    # Ratio of phosphorus (convert to moles) to alkalininty, binned by age
    t = @. !isnan.(sampleage);
    for key in target 
        s = t .& class[key]
        k = invweight_age(sampleage[s])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        c, m, el, eu = bin_bsr_ratios(sampleage[s], mbulk.P2O5[s].*P2O5_to_mol, alk[s], 
            xmin, xmax, nbins,
            x_sigma=ageuncert[s], p=p,
            num_sigma = fill(alk_err, count(s)), denom_sigma = fill(P2O5_err, count(s)),
        )

        simout_m[key] .= m
        simout_e[key] .= [eu el]
    end


## --- Plot 
    p = Plots.palette(colorpalette, 5)
    p2 = Plots.palette(:isoluminant_cm_70_c39_n256, 3, rev=true)
    
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target));
    figname = ("Whole Earth", "Sedimentary", "Igneous", "Volcanic", "Plutonic")
    for i in eachindex(target)
        k = target[i] 

        # Initialize plot
        h = Plots.plot(
            framestyle=:box,
            fontfamily=:Helvetica, 
            # xlabel="Age [Ma.]", ylabel="P / Alk [mol. ratio]",
            # fg_color_legend=:white, legend=:topright, legendfontsize=14,
            labelfontsize=labelfontsize,
            legend=false,
            xlims=(-50, 3850)
        )

        # Get (future) ylims 
        lim_l = round(minimum(simout_m[k] .- simout_e[k][:,1])*0.9, sigdigits=4)
        lim_u = round(maximum(simout_m[k] .+ simout_e[k][:,2])*1.1, sigdigits=4)

        # Color GTS eons 
        Plots.plot!(h, [2500, 3850, 3850, 2500], [lim_l, lim_l, lim_u, lim_u],
            seriestype=:shape, color=p2[1], alpha=0.15, lcolor=:match
        )
        Plots.plot!(h, [541, 2500, 2500, 541], [lim_l, lim_l, lim_u, lim_u],
            seriestype=:shape, color=p2[2], alpha=0.15, lcolor=:match
        )
        Plots.plot!(h, [-50, 541, 541, -50], [lim_l, lim_l, lim_u, lim_u],
            seriestype=:shape, color=p2[3], alpha=0.15, lcolor=:match
        )

        # Plot data and reset y limits
        Plots.plot!(h, c, simout_m[k], yerror=(simout_e[k][:,1], simout_e[k][:,2],), 
            left_margin=(15,:px), bottom_margin=(15,:px), 
            label="$(figname[i])", color=p[i], lcolor=p[i], msc=:auto,
            seriestype=:scatter)
        Plots.ylims!(lim_l, lim_u)
        Plots.annotate!(((0.9, 0.97), (figname[i] * "\nn = $(count(class[k]))", 
            labelfontsize, :right, :top)))
        fig[i] = h
    end

    # Assemble plots
    xlabel!(fig[end], "Age [Ma.]")
    ylabel!(fig[1], "P / Alk [mol. ratio]")
    ylabel!(fig[end-1], "P / Alk [mol. ratio]")
    h = Plots.plot(fig..., layout=(2, 3), size=(1600, 1000), 
        left_margin=(30,:px), bottom_margin=(30,:px)
    )


## --- End of file 