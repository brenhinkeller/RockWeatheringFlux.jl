## --- Set up 
    # DZ and F-org
    
    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots

    # Save figures to: 
    filepath = "results/figures/burial"

    # Definitions 
    xmin, xmax, nbins = 0, 3800, 38
    c = cntr(xmin:100:xmax)
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];
    seds = (minorsed..., :sed)


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

    dz = importdataset("data/dz_voice_2011.csv", ',', importas=:Tuple)


## --- Process data
    # Fraction buried as organic in 100 Myr. bins
    mantle = -5
    carb = sim_carb.m .± nanmean([sim_carb.el sim_carb.eu], dims=2)
    org_corrected = sim_corr.m .± nanmean([sim_corr.el sim_corr.eu], dims=2)

    frog = (mantle .- carb) ./ (org_corrected .- carb)
    frog = (;
        val = Measurements.value.(frog),
        err = Measurements.uncertainty.(frog),
    )

    # DZ in 25 Myr. bins
    c₁,n = bincounts(dz.best_age_2_Ma, xmin, xmax, nbins*4)


## --- Plot 
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        foreground_color_border=:seagreen,
        foreground_color_text=:seagreen,
        foreground_color_axis=:seagreen,
        guidefontcolor=:seagreen,
        xlabel="Carbon Age [Ma.]"
    );
    # plot!(c₁ .+ 350, n, 
    #     seriestype=:bar,
    #     label="Offset DZ Age",
    #     ylabel="Global DZ Abundance",
    #     # xlabel="DZ Age [Ma.]",
    #     color=mineralcolors["zircon"], lcolor=mineralcolors["zircon"], msc=:auto,
    #         y_foreground_color_border=mineralcolors["zircon"],
    #         y_foreground_color_text=mineralcolors["zircon"],
    #         y_foreground_color_axis=mineralcolors["zircon"],
    #         y_guidefontcolor=mineralcolors["zircon"],
    #     # alpha=0.25, lalpha=0.25,
    # )
    plot!(c₁, n, 
    seriestype=:bar,
    label="Measured DZ Age",
    ylabel="Global DZ Abundance",
    # xlabel="DZ Age [Ma.]",
    color=mineralcolors["zircon"], lcolor=mineralcolors["zircon"], msc=:auto,
        y_foreground_color_border=mineralcolors["zircon"],
        y_foreground_color_text=mineralcolors["zircon"],
        y_foreground_color_axis=mineralcolors["zircon"],
        y_guidefontcolor=mineralcolors["zircon"],
    alpha=0.25, linealpha=0,
)
    # vline!([1000])
    plot!(twinx(), c, frog.val, 
        ribbon=2*frog.err,
        label="",
        ylabel="Fraction of Carbon Buried as Organic",
        color=:seagreen,
            foreground_color_border=:seagreen,
            foreground_color_text=:seagreen,
            foreground_color_axis=:seagreen,
            guidefontcolor=:seagreen,
    )
    display(h)
    savefig(h, "$filepath/dz.pdf")


## --- End of File 