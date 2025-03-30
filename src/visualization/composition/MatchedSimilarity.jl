## --- Set up 
    # Compare the age and locations of mapped lithologies to matched geochemical samples.

    # Load data and base packages
    include("Definitions.jl");


## --- Age 
    agediff = abs.(macrostrat.age .- mbulk.Age);
    c, n = bincounts(agediff, 0, 3800, 38)
    # n = float(n) ./ nansum(float(n) .* step(c))
    m = nanmedian(agediff)

    xticklabels = 0:1000:3800
    h1 = plot(c[n .> 0], n[n .> 0], label="",
        seriestype=:bar, barwidths = 100,
        seriescolor=colors_contrast.a, lcolor=:match, alpha=0.5,

        ylabel="Sample Count",
        xlabel="Absolute Age Difference [Myr.]", labelfontsize=14,

        xlims=(0, 3800),
        xticks=(xticklabels, xticklabels), xminorticks=2,

        ylims=(1e-1, maximum(n)*2),
        yscale=:log10, 
        yminorticks=5, tickfontsize=12,

        title="A. Temporal Difference", titlefont=16, titlepos=:left,
        framestyle=:box, 
        fontfamily=:Helvetica, 
    )

    vline!([m], label="Median = $(round(m, sigdigits=3)) Myr.",
        color=:black, linestyle=:dot, linewidth=2, 
        fg_color_legend=:transparent, bg_color_legend=:transparent,
        legendfontsize=12,
    )
    display(h1)
    # savefig(h1, "$filepath/diff_age.pdf")


## --- Location 
    locdiff = haversine.(mbulk.Latitude, mbulk.Longitude, macrostrat.rocklat, 
        macrostrat.rocklon)
    c, n = bincounts(locdiff, 0, 180, 45)
    # n = float(n) ./ nansum(float(n) .* step(c))
    n = float(n) ./ 10^4
    m = nanmedian(locdiff)

    xticklabels = 0:45:180
    h2 = plot(c, n, label="",
        seriestype=:bar, barwidths = 4,
        seriescolor=colors_contrast.b, lcolor=:match, alpha=0.5,

        ylabel="Sample Count [x 10⁴]",
        xlabel="Distance [arc degrees]", labelfontsize=14,

        xlims=(0, 180),
        xticks=(xticklabels, xticklabels), xminorticks=2,

        ylims=(0, maximum(n)*1.05),
        yminorticks=5, tickfontsize=12,

        title="B. Spatial Difference", titlefont=16, titlepos=:left,
        framestyle=:box, 
        fontfamily=:Helvetica, 
    )
    vline!([m], label="Median = $(round(m, sigdigits=3)) arc degrees",
        color=:black, linestyle=:dot, linewidth=2, 
        fg_color_legend=:transparent, bg_color_legend=:transparent,
        legendfontsize=12,
    )
    display(h2)
    # savefig(h2, "$filepath/diff_loc.pdf")


## --- Collect plots
    h = Plots.plot(h1, h2, layout=(1,2), size=(1200, 400),
        left_margin=(40,:px),
        bottom_margin=(40, :px),
    )
    display(h)
    savefig(h, "$filepath/matchingdifference.pdf")
    

## --- End of File