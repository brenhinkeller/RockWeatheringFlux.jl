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


## --- Load data
    # Carbon isotope data 
    fid = h5open("src/visualization/burial/resampled_carbon.h5", "r")
    
    head = keys(fid["vars"]["carb"])
    sim_carb = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["carb"]["$k"]) for k in head)
    
    head = keys(fid["vars"]["org"])
    sim_org = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["org"]["$k"]) for k in head)

    head = keys(fid["vars"]["corrected"])
    sim_corr = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["corrected"]["$k"]) for k in head)
    close(fid)

    # Resampled geochemical data (phosphorus / alkalinity ratios)
    fid = h5open("src/visualization/burial/resampled_geochem.h5", "r")
    head = read(fid["vars"]["ratio"]["head_ratio"])
    head_class = keys(fid["vars"]["ratio"]["data"])
    sim_ratio = NamedTuple{Tuple(Symbol.(head_class))}(NamedTuple{Tuple(Symbol.(head))}(
            read(fid["vars"]["ratio"]["data"]["$k"])[:,i] for i in eachindex(head)
        ) for k in head_class
    )
    # close(fid)

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

    target = keys(sim_ratio)
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
            color=colors[target[i]], lcolor=colors[target[i]], msc=:auto,
                y_foreground_color_border=colors[target[i]],
                y_foreground_color_text=colors[target[i]],
                y_foreground_color_axis=colors[target[i]],
                y_guidefontcolor=colors[target[i]],
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

    h = Plots.plot(figs..., layout=(2,3), size=(1800,800))
    display(h)
    savefig("$filepath/primaryproduction.pdf")


## --- Fraction buried as organic and phosphorus abundance 
    h = plot(
        xlabel="Age [Ma.]", ylabel="Phosphorus [wt.]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );

    target = keys(sim_ratio)
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(figs)
        # Remove Archean sedimentary values because uncertainties are so large they 
        # take out any useful younger signal 
        s = trues(length(sim_ratio[target[i]].c))
        f = sim_wt.cats[target[i]]
        # if target[i] in seds
        #     s[25:38] .= false
        # end

        # Plot data 
        hᵢ = deepcopy(h)
        c,m,e = binmeans(sim_wt.data.Age[f], sim_wt.data.P2O5[f], xmin, xmax, nbins)
        plot!(hᵢ, c[s], m[s],
            yerror=2e[s],
            ribbon=2e[s],
            fillalpha=0.25,
            label="",
            title="$(target[i])",
            color=colors[target[i]], lcolor=colors[target[i]], msc=:auto,
                y_foreground_color_border=colors[target[i]],
                y_foreground_color_text=colors[target[i]],
                y_foreground_color_axis=colors[target[i]],
                y_guidefontcolor=colors[target[i]],
            markershape=:circle,
            seriestype=:scatter,
            linewidth=2,
        )
        plot!(twinx(), c, frog.val, 
            yerror=2*frog.err,
            label="", ylabel="Fraction Buried as Organic",
            color=:black,
        )
        figs[i] = hᵢ
    end

    h = plot(figs..., layout=(2,3), size=(1800,800))
    display(h)
    savefig("$filepath/primaryproduction_Pwt.pdf")


## --- End of file 