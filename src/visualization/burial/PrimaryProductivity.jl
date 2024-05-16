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
    close(fid)

    # Resampled geochemical data (phosphorus / alkalinity ratios)
    fid = h5open("src/visualization/burial/resampled_geochem.h5", "r")
    head = read(fid["vars"]["ratio"]["head_ratio"])
    head_class = keys(fid["vars"]["ratio"]["data"])
    sim_ratio = NamedTuple{Tuple(Symbol.(head_class))}(NamedTuple{Tuple(Symbol.(head))}(
            read(fid["vars"]["ratio"]["data"]["$k"])[:,i] for i in eachindex(head)
        ) for k in head_class
    )
    close(fid)


## --- Calculate fraction of carbon buried as organic 
    mantle = -5.5
    
    c,m,e = binmeans(sim_carb.age, sim_carb.d13c_carb, xmin, xmax, nbins, relbinwidth=2)
    carbonate = m .± e
    
    c,m,e = binmeans(sim_org.age, sim_org.d13c_org_corrected, xmin, xmax, nbins)
    organic = m .± e

    frog = (mantle .- carbonate) ./ (organic .- carbonate)
    frog = (;
        val = Measurements.value.(frog),
        err = Measurements.uncertainty.(frog),
    )


## --- [PLOT] Fraction buried as organic with phosphorus / alkalinity ratios 
    h = plot(
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
        # Remove Archean sedimentary uncertainties because they're so large they erase the 
        # younger signal we can constrain 
        el, eu = 2*sim_ratio[target[i]].el, 2*sim_ratio[target[i]].eu
        if target[i] in seds
            el[25:38] .= sim_ratio[target[i]].m[25:38] * 0.05
            eu[25:38] .= sim_ratio[target[i]].m[25:38] * 0.05
        end

        # Plot data 
        hᵢ = deepcopy(h)
        plot!(hᵢ, sim_ratio[target[i]].c, sim_ratio[target[i]].m,
            yerror=(el,eu),
            ribbon=(el,eu),
            fillalpha=0.25,
            label="",
            title="$(target[i])",
            color=colors[target[i]], lcolor=colors[target[i]], msc=:auto,
                y_foreground_color_border=colors[target[i]],
                y_foreground_color_text=colors[target[i]],
                y_foreground_color_axis=colors[target[i]],
                y_guidefontcolor=colors[target[i]],
            markershape=:circle,
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
    savefig("$filepath/primaryproduction.pdf")


## --- End of file 