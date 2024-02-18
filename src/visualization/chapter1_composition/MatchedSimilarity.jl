## --- Set up 
    # Compare the age and locations of mapped lithologies to matched geochemical samples.

    # Load data and base packages
    include("Definitions.jl");


## --- Age 
    agediff = macrostrat.age .- mbulk.Age;
    c, n = bincounts(agediff, -3800, 3800, 38*2)
    m = nanmedian(agediff)

    h = plot(
        framestyle=:box, 
        fontfamily=:Helvetica, 
        xlabel="Age Difference [Myr.]", ylabel="Abundance",
        ylims=(0.1, round(Int, maximum(n)*1.5)), yticks=10.0.^(0:2:5),
        xlims=(-1100, 3000),
    )
    Plots.plot!(h, c[n .> 0], n[n .> 0], label="",
        seriestype=:bar, 
        yscale=:log10,
        color=colors.evap, lcolor=:match,
        barwidths = 100,
    )
    vline!([0], linestyle=:solid, color=:black, linewidth=1, label="")
    Plots.annotate!(((0.03, 0.97), ("Geochemical\nage is older", 12, :left, :top)))
    Plots.annotate!(((0.97, 0.97), ("Mapped age\nis older", 12, :right, :top)))

    vline!([m], color=:black, linestyle=:dash, linewidth=2, label="")
    Plots.annotate!([(m*40, 0.25, text("Median = $(round(m, sigdigits=3))", 9, :left, 
        :bottom, :black, rotation=90))]
    )

    jump = 500
    xticks!((xlims(h)[1]jump)*jump:jump:(xlims(h)[2]÷jump-1)*jump)

    display(h)
    savefig(h, "$filepath/agediff.pdf")


## --- Location 
    locdiff = haversine.(mbulk.Latitude, mbulk.Longitude, macrostrat.rocklat, 
        macrostrat.rocklon
    )
    c, n = bincounts(locdiff, 0, 180, 90)
    m = nanmedian(locdiff)

    h = plot(
        framestyle=:box, 
        fontfamily=:Helvetica, 
        xlabel="Distance [arc degrees]", ylabel="Abundance",
        ylims=(0.1, round(Int, maximum(n)*1.05)),
        xlims=(-2, 180),
    )
    Plots.plot!(h, c, n, label="",
        seriestype=:bar, 
        color=colors.sed, lcolor=:match,
        barwidths = 2,
    )
    vline!([m], color=:black, linestyle=:dash, linewidth=2, label="")
    Plots.annotate!([(m*0.9, 500, text("Median = $(round(m, sigdigits=3))", 9, :left, 
        :bottom, :black, rotation=90))]
    )
    display(h)
    savefig(h, "$filepath/locdiff.pdf")


## --- End of File