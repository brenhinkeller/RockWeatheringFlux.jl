## --- Set up 
    # Compare phosphorus / alkalinity to fraction of carbon buried as organic 
    
    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots

    # Save figures to: 
    filepath = "results/figures/burial"

    # Definitions 
    xmin, xmax, nbins = 0, 3800, 38
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];
    seds = (minorsed..., :sed)

    colors2 = (
        sed = :mediumturquoise,
        carb = :dodgerblue,
        shale = :slategrey,
        siliciclast = :blueviolet,
        ign = :crimson,
        plut = :mediumvioletred,
        volc = :deeppink,
        basalt = :pink,
    )


## --- Load data
    # Carbon isotope data 
    fid = h5open("src/burial/resampled_carbon.h5", "r")
    
    head = keys(fid["vars"]["carb"])
    sim_carb = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["carb"]["$k"]) for k in head)
    
    head = keys(fid["vars"]["org"])
    sim_org = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["org"]["$k"]) for k in head)

    head = keys(fid["vars"]["corrected"])
    sim_corr = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["corrected"]["$k"]) for k in head)
    close(fid)

    # Resampled geochemical data (phosphorus / alkalinity ratios)
    fid = h5open("src/burial/resampled_geochem.h5", "r")
    head = read(fid["vars"]["ratio"]["head_ratio"])
    head_class = keys(fid["vars"]["ratio"]["data"])
    sim_ratio = NamedTuple{Tuple(Symbol.(head_class))}(NamedTuple{Tuple(Symbol.(head))}(
            read(fid["vars"]["ratio"]["data"]["$k"])[:,i] for i in eachindex(head)
        ) for k in head_class
    )
    # close(fid)

    # Resampled geochemical data (phosphorus / alkalinity mole abundance)
    head_data = read(fid["vars"]["mole"]["head_data"])
    head_cats = read(fid["vars"]["mole"]["head_cats"])
    data = read(fid["vars"]["mole"]["data"]["data"])
    cats = read(fid["vars"]["mole"]["cats"]["filter_cats"])
    cats = @. cats > 0
    sim_mol = (;
        data = NamedTuple{Tuple(Symbol.(head_data))}([data[:,i] for i in eachindex(head_data)]),
        cats = NamedTuple{Tuple(Symbol.(head_cats))}([cats[:,i] for i in eachindex(head_cats)])
    )

    # Resampled geochemical data (phosphorus)
    head_data = read(fid["vars"]["wt"]["head_data"])
    head_cats = read(fid["vars"]["wt"]["head_cats"])
    data = read(fid["vars"]["wt"]["data"]["data"])
    cats = read(fid["vars"]["wt"]["cats"]["match_cats"])
    cats = @. cats > 0
    sim_wt = (;
        data = NamedTuple{Tuple(Symbol.(head_data))}([data[:,i] for i in eachindex(head_data)]),
        cats = NamedTuple{Tuple(Symbol.(head_cats))}([cats[:,i] for i in eachindex(head_cats)])
    )

    close(fid)


## --- Calculate fraction of carbon buried as organic 
    mantle = -5
    carb = sim_carb.m .± nanmean([sim_carb.el sim_carb.eu], dims=2)
    # org = sim_org.m .± nanmean([sim_org.el sim_org.eu], dims=2)
    org_corrected = sim_corr.m .± nanmean([sim_corr.el sim_corr.eu], dims=2)

    frog = (mantle .- carb) ./ (org_corrected .- carb)
    frog = (;
        val = Measurements.value.(frog),
        err = Measurements.uncertainty.(frog),
    )


## --- [PLOT] Fraction buried as organic with phosphorus / alkalinity ratios 
    h = Plots.plot(
        xlabel="Age [Ma.]", ylabel="P / Alk [mol. ratio]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );

    # target = keys(sim_ratio)
    target = (
        :ign, :volc, :basalt, :plut, 
        :sed, :siliciclast, :shale, :carb
    )
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(figs)
        # Remove Archean sedimentary values because uncertainties are so large they 
        # take out any useful younger signal 
        s = trues(length(sim_ratio[target[i]].c))
        if target[i] in seds
            s[25:38] .= false
        end

        # Plot data 
        hᵢ = deepcopy(h)
        Plots.plot!(hᵢ, sim_ratio[target[i]].c[s], sim_ratio[target[i]].m[s],
            yerror=(2*sim_ratio[target[i]].el[s], 2*sim_ratio[target[i]].eu[s]),
            ribbon=(2*sim_ratio[target[i]].el[s], 2*sim_ratio[target[i]].eu[s]),
            fillalpha=0.25,
            label="",
            title="$(target[i])",
            color=colors2[target[i]], lcolor=colors2[target[i]], msc=:auto,
                y_foreground_color_border=colors2[target[i]],
                y_foreground_color_text=colors2[target[i]],
                y_foreground_color_axis=colors2[target[i]],
                y_guidefontcolor=colors2[target[i]],
            markershape=:circle,
            seriestype=:scatter,
            linewidth=2,
        )
        Plots.plot!(twinx(), c, frog.val, 
            ribbon=2*frog.err,
            label="", ylabel="Fraction Buried as Organic",
            color=:seagreen,
        )
        figs[i] = hᵢ
    end

    h = Plots.plot(figs..., layout=(2,4), size=(2000,900))
    display(h)
    savefig("$filepath/primaryproduction.pdf")

    # # Check y limits
    # for i in eachindex(figs)
    #     println(ylims(figs[i])[2] - ylims(figs[i])[1])
    # end


## --- [PLOT] Phosphourus / alkalinity breakdown (who's driving the bus?)
    h = Plots.plot(
        xlabel="Age [Ma.]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(figs)
        # Remove Archean sedimentary values because uncertainties are so large they 
        # take out any useful younger signal 
        s = trues(nbins)
        if target[i] in seds
            s[25:38] .= false
        end

        # Set up 
        hᵢ = deepcopy(h)
        f = sim_mol.cats[target[i]]

        # Phosphorus 
        c,m,e = binmeans(sim_mol.data.Age[f], sim_mol.data.P[f], xmin, xmax, nbins)
        Plots.plot!(hᵢ, c[s], m[s],
            yerror=2e[s],
            ribbon=2e[s],
            fillalpha=0.25,
            label="",
            ylabel="Phosphorus [mol.]",
            title="$(target[i])",
            color=colors2[target[i]], lcolor=colors2[target[i]], msc=:auto,
                y_foreground_color_border=colors2[target[i]],
                y_foreground_color_text=colors2[target[i]],
                y_foreground_color_axis=colors2[target[i]],
                y_guidefontcolor=colors2[target[i]],
            markershape=:circle,
            seriestype=:scatter,
            linewidth=2,
            xlims=(xmin, xmax),
        )

        # Alkalinity 
        c,m,e = binmeans(sim_mol.data.Age[f], sim_mol.data.Alk[f], xmin, xmax, nbins)
        Plots.plot!(twinx(), c[s], m[s],
            yerror=2e[s],
            ribbon=2e[s],
            fillalpha=0.25,
            label="",
            ylabel="Alkalinity [mol.]",
            color=:black, lcolor=:black, msc=:auto,
            # color=colors2[target[i]], lcolor=colors2[target[i]], msc=:auto,
                # y_foreground_color_border=colors2[target[i]],
                # y_foreground_color_text=colors2[target[i]],
                # y_foreground_color_axis=colors2[target[i]],
                # y_guidefontcolor=colors2[target[i]],
            markershape=:circle,
            seriestype=:scatter,
            linewidth=2,
            xlims=(xmin, xmax),
        )

        # Save to array 
        figs[i] = hᵢ
    end
    h = Plots.plot(figs..., layout=(2,4), size=(2000,900))
    display(h)
    savefig(h, "$filepath/p_alk_breakdown.pdf")
    

## --- End of file 