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

    ages = (
        archean = (4000,2500),
        proterozoic = (2500,541),
        phanerozoic = (541,0),
    )


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
    c, n = bincounts(mbulk.SiO2[macro_cats.metaign], SiO2min, SiO2max, nbins_matched)
    n₁ = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n₁, seriestype=:bar, color=colors.metaign, linecolor=colors.metaign,
        label="Matched metaign", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims = (0,0.09)
    )

    c, n = bincounts(mbulk.SiO2[macro_cats.met], SiO2min, SiO2max, nbins_matched)
    n₂ = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n₂, seriestype=:bar, color=colors.met, linecolor=colors.met,
        label="Matched met", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims = (0,0.09)
    )
    
    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
    )
    display(h)


## --- Resample matched sample distributions
    # Universal filters and uncertainties
    t = @. !isnan(macrostrat.rocklat) & !isnan(macrostrat.rocklat) & !isnan(macrostrat.age)
    ageuncert = (macrostrat.agemax .- macrostrat.agemin) ./ 2
    for i in eachindex(ageuncert)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), macrostrat.age[i]*0.05, ageuncert[i])
    end

    # Metamorphic rocks
    test = t .& macro_cats.met
    k = invweight(macrostrat.rocklat[test], macrostrat.rocklon[test],macrostrat.age[test])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    data = [macrostrat.age[test] mbulk.SiO2[test]]
    uncertainty = [ageuncert[test] fill(0.01, count(test))]

    simmet = bsresample(data, uncertainty, Int(1e6), p)

    # Metaigneous rocks 
    test = t .& macro_cats.metaign
    k = invweight(macrostrat.rocklat[test], macrostrat.rocklon[test],macrostrat.age[test])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    data = [macrostrat.age[test] mbulk.SiO2[test]]
    uncertainty = [ageuncert[test] fill(0.01, count(test))]

    simmetaign = bsresample(data, uncertainty, Int(1e6), p)


## --- Silica distributions of resampled data
    c, n = bincounts(simmet[:,2], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        
        label="Resampled Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    c, n = bincounts(simmetaign[:,2], SiO2min, SiO2max, nbins_matched)
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
    t = @. 4000 > simmet[:,1] >= 2500;
    c, n = bincounts(simmet[:,2][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Archean Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )
    t = @. 4000 > simmetaign[:,1] >= 2500;
    c, n = bincounts(simmetaign[:,2][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Archean Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    # Proterozoic
    t = @. 2500 > simmet[:,1] >= 541;
    c, n = bincounts(simmet[:,2][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h3 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Proterozoic Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 2500 > simmetaign[:,1] >= 541;
    c, n = bincounts(simmetaign[:,2][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h4 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Proterozoic Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    # Phanerozoic
    t = @. 541 > simmet[:,1] >= 0;
    c, n = bincounts(simmet[:,2][t], SiO2min, SiO2max, nbins_matched)
    n = float(n) ./ nansum(float(n) .* step(c))
    h5 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Phanerozoic Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins_matched),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 541 > simmetaign[:,1] >= 0;
    c, n = bincounts(simmetaign[:,2][t], SiO2min, SiO2max, nbins_matched)
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


## --- Silica content timeseries
    c,m,e = binmeans(simmet[:,1], simmet[:,2], 0,3800, 38)
    h1 = Plots.plot(c, m, yerror=e, label="Resampled Metamorphic", 
        markershape=:circle,
        color=colors.met, linecolor=colors.met, msc=colors.met, 
    )

    c,m,e = binmeans(simmetaign[:,1], simmetaign[:,2], 0,3800, 38)
    h2 = Plots.plot(c, m, yerror=e, label="Resampled Metaigneous", 
        markershape=:circle,
        color=colors.metaign, linecolor=colors.metaign, msc=colors.metaign, 
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
    )


## --- End of file 