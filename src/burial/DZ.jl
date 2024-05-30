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
        xlabel="Age [Ma.]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
    );
    plot!(c₁, n, 
        seriestype=:bar,
        label="",
        ylabel="Global DZ Abundance",
        color=mineralcolors["zircon"], lcolor=mineralcolors["zircon"], msc=:auto,
            y_foreground_color_border=mineralcolors["zircon"],
            y_foreground_color_text=mineralcolors["zircon"],
            y_foreground_color_axis=mineralcolors["zircon"],
            y_guidefontcolor=mineralcolors["zircon"],
    )
    plot!(twinx(), c, frog.val, 
        ribbon=2*frog.err,
        label="",
        ylabel="Fraction of Carbon Buried as Organic",
        color=:seagreen,
            y_foreground_color_border=:seagreen,
            y_foreground_color_text=:seagreen,
            y_foreground_color_axis=:seagreen,
            y_guidefontcolor=:seagreen,
    )
    display(h)
    savefig(h, "$filepath/dz.pdf")


## --- End of File 