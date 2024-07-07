## --- Set up 
    # Better characterize the interaction between siliciclastic and volcanic rocks 
    
    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots

    # Save figures to: 
    filepath = "results/figures/burial"

    # Definitions 
    xmin, xmax, nbins = 0, 3800, 38


## --- Load data 
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

## --- [PLOT] Siliciclastic / volcanic P/Alk covariance 
    h = plot(
        xlabel="Age [Ma.]", 
        ylabel="P / Alk [mol. ratio]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        # left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
        size=(600,800),
        # ylims=(0,2.6),
    );
    k = :siliciclast
    plot!(h, sim_ratio[k].c, sim_ratio[k].m,
        yerror=(2*sim_ratio[k].el, 2*sim_ratio[k].eu),
        ribbon=(2*sim_ratio[k].el, 2*sim_ratio[k].eu),
        fillalpha=0.25,
        label="",
        ylabel="Siliciclastic P/Alk [mol.]",
        color=colors2[k], lcolor=colors2[k], msc=:auto,
            y_foreground_color_border=colors2[k],
            y_foreground_color_text=colors2[k],
            y_foreground_color_axis=colors2[k],
            y_guidefontcolor=colors2[k],
        markershape=:circle,
        linewidth=2,
    )
    k = :volc
    plot!(twinx(), sim_ratio[k].c, sim_ratio[k].m,
        yerror=(2*sim_ratio[k].el, 2*sim_ratio[k].eu),
        ribbon=(2*sim_ratio[k].el, 2*sim_ratio[k].eu),
        fillalpha=0.25,
        label="",
        ylabel="Volcanic P/Alk [mol.]",
        color=colors2[k], lcolor=colors2[k], msc=:auto,
            y_foreground_color_border=colors2[k],
            y_foreground_color_text=colors2[k],
            y_foreground_color_axis=colors2[k],
            y_guidefontcolor=colors2[k],
        markershape=:circle,
        linewidth=2,
        ylims=(0.002, 0.01),
    )
    display(h)
    savefig(h, "$filepath/phos_siliciclast_volc.pdf")

## --- TiO2?
    mafic = @. (43 < sim_wt.data.SiO2 < 51) .& sim_wt.cats.volc;
    c,m,e = binmeans(sim_wt.data.Age[mafic], sim_wt.data.TiO2[mafic], xmin,xmax,nbins)
    plot(c,m,yerror=2e, 
        label="",
        color=:red, lcolor=:red, msc=:auto, 
            y_foreground_color_border=:red,
            y_foreground_color_text=:red,
            y_foreground_color_axis=:red,
            y_guidefontcolor=:red,
        ylabel="Mafic Volcanic TiO2",
        linewidth=2,
        xlabel="Age [Ma.]",
        size=(600,600),
    )

    f = sim_wt.cats.siliciclast
    c,m,e = binmeans(sim_wt.data.Age[f], sim_wt.data.TiO2[f], xmin,xmax,nbins)
    plot!(twinx(), c,m,yerror=2e, 
        label="True age sed", 
        legend=:topright,
        color=:blue, lcolor=:blue, msc=:auto, 
        linestyle=:dash,
        y_foreground_color_border=:seagreen,
        y_foreground_color_text=:seagreen,
        y_foreground_color_axis=:seagreen,
        y_guidefontcolor=:seagreen,
        # ylabel="Siliciclastic TiO2"
        # linewidth=2,
    )
    
    plot!(twinx(), c.+250,m,yerror=2e, 
        label="250 Ma. offset sed", 
        legend=:top,
        color=:seagreen, lcolor=:seagreen, msc=:auto, 
            y_foreground_color_border=:seagreen,
            y_foreground_color_text=:seagreen,
            y_foreground_color_axis=:seagreen,
            y_guidefontcolor=:seagreen,
        ylabel="Siliciclastic TiO2",
        linewidth=2,
    )



## --- End of file