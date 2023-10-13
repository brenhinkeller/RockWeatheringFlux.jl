## --- Closer look at metamorphic rocks (specifically met and metaign)
## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using DelimitedFiles

    using StatsPlots

    # Local utilities
    using Static, LoopVectorization, Measurements
    include("../utilities/Utilities.jl")

    # Settings and local definitions
    SiO2min, SiO2max = 0, 100
    nbins = 200
    nbins_matched = Int(nbins/2)

    types = [:metaign, :met,]
    labels = string.(types)


## --- Load data
    # Resampled EarthChem data
    fid = h5open("output/resampled/resampled.h5", "r")
    header = read(fid["vars"]["header"])
    i = findfirst(x -> x=="SiO2", header)
    rocktypes = keys(fid["vars"]["data"])

    bsrsilica = NamedTuple{Tuple(Symbol.(rocktypes))}([read(fid["vars"]["data"][r]["data"])[:,i]
        for r in rocktypes
    ])
    close(fid)

    # Matched EarthChem data
    fid = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(fid)

    # Macrostrat data
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    # For this analysis, we want metamorphic rocks to include metaigns, but not metaseds
    macro_cats.met .|= macro_cats.metaign


## --- Show distribution of matched metamorphic rocks
    c, n = bincounts(mbulk.SiO2[macro_cats.met], SiO2min, SiO2max, nbins_matched)
    n₁ = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n₂, seriestype=:bar, color=colors.met, linecolor=colors.met,
        label="Matched met", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims = (0,0.09)
    )

    c, n = bincounts(mbulk.SiO2[macro_cats.metaign], SiO2min, SiO2max, nbins_matched)
    n₂ = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n₂, seriestype=:bar, color=colors.metaign, linecolor=colors.metaign,
        label="Matched metaign", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims = (0,0.09)
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        legend=:topleft, 
    )
    display(h)


## --- Resample matched sample distributions
    # Universal filters and uncertainties
    t = @. !isnan(macrostrat.rocklat) & !isnan(macrostrat.rocklat) & !isnan(macrostrat.age)
    ageuncert = (macrostrat.agemax .- macrostrat.agemin) ./ 2
    for i in eachindex(ageuncert)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), macrostrat.age[i]*0.05, ageuncert[i])
    end

    # Define order of data (and resampled output)
    latc, lonc, agec, SiO2c = 1,2,3,4

    # Metamorphic rocks
    test = t .& macro_cats.met
    data = [macrostrat.rocklat[test] macrostrat.rocklon[test] macrostrat.age[test] mbulk.SiO2[test]]
    uncertainty = [zeros(count(test)) zeros(count(test)) ageuncert[test] fill(0.01, count(test))]
    simmet = bsresample(data, uncertainty, Int(1e6), ones(count(test)))

    # Metaigneous rocks 
    test = t .& macro_cats.metaign
    data = [macrostrat.rocklat[test] macrostrat.rocklon[test] macrostrat.age[test] mbulk.SiO2[test]]
    uncertainty = [zeros(count(test)) zeros(count(test)) ageuncert[test] fill(0.01, count(test))]
    simmetaign = bsresample(data, uncertainty, Int(1e6), ones(count(test)))


## --- Silica distributions of resampled data
    c, n = bincounts(simmet[:,SiO2c], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        
        label="Resampled Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    c, n = bincounts(simmetaign[:,SiO2c], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Resampled Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        legend=:topleft, fg_color_legend=:white,
    )


## --- Distribution of samples over time
    # Archean
    t = @. 4000 > simmet[:,agec] >= 2500;
    c, n = bincounts(simmet[:,SiO2c][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Archean Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 4000 > simmetaign[:,agec] >= 2500;
    c, n = bincounts(simmetaign[:,SiO2c][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Archean Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    # Proterozoic
    t = @. 2500 > simmet[:,agec] >= 541;
    c, n = bincounts(simmet[:,SiO2c][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h3 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Proterozoic Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 2500 > simmetaign[:,agec] >= 541;
    c, n = bincounts(simmetaign[:,SiO2c][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h4 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Proterozoic Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    # Phanerozoic
    t = @. 541 > simmet[:,agec] >= 0;
    c, n = bincounts(simmet[:,SiO2c][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h5 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Phanerozoic Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 541 > simmetaign[:,agec] >= 0;
    c, n = bincounts(simmetaign[:,SiO2c][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h6 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Phanerozoic Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    # All together
    h = Plots.plot(h1, h3, h5, h2, h4, h6, layout=(2, 3), size=(1800,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        legend=:topleft, fg_color_legend=:white,
    )


## --- Distribution of samples over time
    # Metamorphic
    archean = @. 4000 > simmet[:,agec] >= 2500;
    proterozoic = @. 2500 > simmet[:,agec] >= 541;
    phanerozoic = @. 541 > simmet[:,agec] >= 0;

    h1 = Plots.plot(simmet[:,SiO2c][archean], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins_matched/2),
        linecolor=:red, linewidth=3, 
        label="Archean", barwidths = ((SiO2max-SiO2min)/nbins_matched),
    )
    Plots.plot!(simmet[:,SiO2c][proterozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins_matched/2),
        linecolor=:darkorange, linewidth=3, 
        label="Proterozoic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
    )
    Plots.plot!(simmet[:,SiO2c][phanerozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins_matched/2),
        linecolor=:forestgreen, linewidth=3, 
        label="Phanerozoic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
    )
    display(h1)

    # Metaigneous
    archean = @. 4000 > simmetaign[:,agec] >= 2500;
    proterozoic = @. 2500 > simmetaign[:,agec] >= 541;
    phanerozoic = @. 541 > simmetaign[:,agec] >= 0;

    h1 = Plots.plot(simmetaign[:,SiO2c][archean], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins_matched/2),
        linecolor=:red, linewidth=3, 
        label="Archean", barwidths = ((SiO2max-SiO2min)/nbins_matched),
    )
    Plots.plot!(simmetaign[:,SiO2c][proterozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins_matched/2),
        linecolor=:darkorange, linewidth=3, 
        label="Proterozoic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
    )
    Plots.plot!(simmetaign[:,SiO2c][phanerozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins_matched/2),
        linecolor=:forestgreen, linewidth=3, 
        label="Phanerozoic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
    )
    display(h1)


## --- Silica content timeseries
    c,m,e = binmeans(simmet[:,agec], simmet[:,SiO2c], 0,3800, 38)
    h1 = Plots.plot(c, m, yerror=e, label="Resampled Metamorphic", 
        markershape=:circle,
        color=colors.met, linecolor=colors.met, msc=colors.met, 
    )

    c,m,e = binmeans(simmetaign[:,agec], simmetaign[:,SiO2c], 0,3800, 38)
    h2 = Plots.plot(c, m, yerror=e, label="Resampled Metaigneous", 
        markershape=:circle,
        color=colors.metaign, linecolor=colors.metaign, msc=colors.metaign, 
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
    )


## --- Archean metamorphic rocks in Australia
    cont = find_geolcont(simmet[:,latc], simmet[:,lonc])
    aus = cont .== 5

    t = @. 4000 > simmet[:,agec] >= 2500;
    c, n = bincounts(simmet[:,SiO2c][t .& aus], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Archean + Australian Met", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 4000 > simmetaign[:,agec] >= 2500;
    c, n = bincounts(simmetaign[:,SiO2c][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Archean + Australian Metaign", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )


## --- Relative contributions of geologic eons to each bin 
    # Metaigneous
    archean = @. 4000 > simmetaign[:,agec] >= 2500;
    proterozoic = @. 2500 > simmetaign[:,agec] >= 541;
    phanerozoic = @. 541 > simmetaign[:,agec] >= 0;

    c, nar = bincounts(simmetaign[:,SiO2c][archean], SiO2min, SiO2max, nbins_matched)
    c, npr = bincounts(simmetaign[:,SiO2c][proterozoic], SiO2min, SiO2max, nbins_matched)
    c, nph = bincounts(simmetaign[:,SiO2c][phanerozoic], SiO2min, SiO2max, nbins_matched)

    h = groupedbar([nar npr nph], bar_position=:stack, 
        label=["Archean" "Proterozoic" "Phanerozoic"], xlabel="SiO₂ [wt.%]", ylabel="Weight",
        barwidths = ((SiO2max-SiO2min)/nbins_matched),
        framestyle=:box, legend=:topright,
        color=[:red :darkorange :forestgreen],
    )

    x = SiO2min:10:SiO2max
    Plots.xticks!(x, string.(x))


## --- Resampled matched sample density 
    # histogram2d(simmet[:,SiO2c], simmet[:,agec], ylims=(0,3800), xlims=(40,80), bins=100)
    # histogram2d(mbulk.SiO2, macrostrat.age, ylims=(0,3800), xlims=(40,80), bins=100)

    # For each age bin, need to get and normalize a histogram of silica content. Then normalize
    # This isn't doing what I want it to and part of that is the x/y distinction
    xmin, xmax, xbins = 40, 80, 40
    # xedges = xmin:(xmax-xmin)/xbins:xmax

    ymin, ymax, ybins = 0, 3800, 38
    yedges = ymin:(ymax-ymin)/ybins:ymax

    out = zeros(xbins, ybins)
    for i = 1:ybins
        # Filter for samples in this age bin
        t = @. yedges[i] <= macrostrat.age[macro_cats.ign] < yedges[i+1]

        # Count and normalize distribution of silica
        c, n = bincounts(mbulk.SiO2[macro_cats.ign][t], xmin, xmax, xbins)
        n = float(n) ./ nansum(float(n) .* step(c))

        # Put in array for that age bin
        out[:,i] .= n
    end

    zeronan!(out)
    histogram2d(out)

## --- End of file 