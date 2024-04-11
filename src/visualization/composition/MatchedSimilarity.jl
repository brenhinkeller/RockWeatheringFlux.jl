## --- Set up 
    # Compare the age and locations of mapped lithologies to matched geochemical samples.

    # Load data and base packages
    include("Definitions.jl");


## --- Age 
    agediff = abs.(macrostrat.age .- mbulk.Age);
    c, n = bincounts(agediff, 0, 3800, 38*2)
    n = float(n) ./ nansum(float(n) .* step(c))
    m = nanmedian(agediff)

    h1 = plot(
        framestyle=:box, 
        fontfamily=:Helvetica, 
        xlabel="Absolute Age Difference [Myr.]", ylabel="Logarithmic Relative Abundance",
        ylims=(9e-9,0.02),
        yticks=false,
        xlims=(-100, 3800),
        left_margin=(15,:px)
    )
    Plots.plot!(h1, c[n .> 0], n[n .> 0], label="",
        seriestype=:bar, 
        yscale=:log10,
        color=colors_contrast.a, lcolor=:match,
        barwidths = 100,
    )

    vline!([m], color=:black, linestyle=:dot, linewidth=2, label="")
    Plots.annotate!([(m*8, 2e-8, text("Median = $(round(m, sigdigits=3))", 9, :left, 
        :bottom, :black, rotation=90))]
    )

    jump = 500
    xticks!((xlims(h1)[1]jump)*jump:jump:(xlims(h1)[2]÷jump)*jump)

    display(h1)
    savefig(h1, "$filepath/diff_age.pdf")


## --- Location 
    locdiff = haversine.(mbulk.Latitude, mbulk.Longitude, macrostrat.rocklat, 
        macrostrat.rocklon
    )
    c, n = bincounts(locdiff, 0, 180, 45)
    n = float(n) ./ nansum(float(n) .* step(c))
    m = nanmedian(locdiff)

    h2 = plot(
        framestyle=:box, 
        fontfamily=:Helvetica, 
        xlabel="Distance [arc degrees]", # ylabel="Abundance",
        # ylims=(0.1, round(Int, maximum(n)*1.05)),
        ylims=(-1e-4,0.026),
        yticks=false,
        xlims=(-2, 180),
    )
    Plots.plot!(h2, c, n, label="",
        seriestype=:bar, 
        color=colors_contrast.b, lcolor=:match,
        barwidths = 4,
    )
    vline!([m], color=:black, linestyle=:dot, linewidth=2, label="")
    Plots.annotate!([(m*0.9,0.0026, text("Median = $(round(m, sigdigits=3))", 9, :left, 
        :bottom, :black, rotation=90))]
    )
    display(h2)
    savefig(h2, "$filepath/diff_loc.pdf")


## --- Collect plots for LaTeX placeholder

    h = Plots.plot(h1, h2, layout=(1,2), size=(1200,400))
    display(h)
    savefig(h, "$filepath/diff_combo.pdf")

## --- End of File