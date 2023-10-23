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
    nbins = 100


## --- Load data
    # Matched EarthChem data
    fid = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])

    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] 
        for i in eachindex(header)])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Bulk rock name, type, and material
    header = read(fid["bulktext"]["sampledata"]["header"])
    index = read(fid["bulktext"]["sampledata"]["index"])
    target = ["Rock_Name", "Type", "Material"]
    targetind = [findall(==(i), header)[1] for i in target]

    mbulktext = NamedTuple{Tuple(Symbol.(target))}(
        [lowercase.(read(fid["bulktext"]["sampledata"]["elements"][target[i]]))[
            index[:,targetind[i]]][bulkidx[t]] for i in eachindex(target)]
    )
    bulktext = NamedTuple{Tuple(Symbol.(target))}(
        [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] 
            for i in eachindex(target)]
    )
    close(fid)

    # Macrostrat data
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        rocktype = read(fid["vars"]["rocktype"])[t],
        rockname = read(fid["vars"]["rockname"])[t],
        rockdescrip = read(fid["vars"]["rockdescrip"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] 
        for i in eachindex(header)])
    close(fid)

    # For this analysis, we want metamorphic rocks to include metaigns, but not metaseds
    macro_cats.met .|= macro_cats.metaign


## --- Get Macrostrat rock names for the Archean felsic mode 
    archean = @. 4000 > macrostrat.age >= 2500;
    felsic = @. mbulk.SiO2 > 60;
    old_metaigns = macrostrat.rocktype[archean .& macro_cats.metaign .& felsic]
    felsicnames = unique(old_metaigns)

    target = [count(==(i), old_metaigns) for i in felsicnames]
    
    p = reverse(sortperm(target))
    # display([felsicnames[p][1:15] target[p][1:15]])

    n = length(old_metaigns)
    t = @. target[p] > (ceil(Int, n * 0.01));
    display([felsicnames[p][t] target[p][t]])


## --- Show distribution of matched metamorphic rocks
    c, n = bincounts(mbulk.SiO2[macro_cats.metased], SiO2min, SiO2max, nbins)
    n₁ = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n₁, seriestype=:bar, color=colors.met, linecolor=colors.met,
        label="Matched metamorphic", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims = (0,0.09)
    )

    c, n = bincounts(mbulk.SiO2[macro_cats.metaign], SiO2min, SiO2max, nbins)
    n₂ = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n₂, seriestype=:bar, color=colors.metaign, linecolor=colors.metaign,
        label="Matched metaigneous", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims = (0,0.09)
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        legend=:topleft, fg_color_legend=:white,
    )
    display(h)

    # Modal sample analysis for metamorphic distributions
    geochemkeys, = get_elements()
    sameindex(bulkidx[t][macro_cats.met], geochemkeys, (25,100,75), bulk, bulktext);
    sameindex(bulkidx[t][macro_cats.metaign], geochemkeys, (25,100,75), bulk, bulktext);
    

## --- Resample matched sample distributions (defacto spatial)
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
    c, n = bincounts(simmet[:,SiO2c], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        
        label="Resampled Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    c, n = bincounts(simmetaign[:,SiO2c], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Resampled Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        legend=:topleft, fg_color_legend=:white,
    )


## --- Resample matched sample distributions (temporal)
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

    k = invweight_age(macrostrat.age[test])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    simmet_t = bsresample(data, uncertainty, Int(1e6), p)

    # Metaigneous rocks 
    test = t .& macro_cats.metaign
    data = [macrostrat.rocklat[test] macrostrat.rocklon[test] macrostrat.age[test] mbulk.SiO2[test]]
    uncertainty = [zeros(count(test)) zeros(count(test)) ageuncert[test] fill(0.01, count(test))]

    k = invweight_age(macrostrat.age[test])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    simmetaign_t = bsresample(data, uncertainty, Int(1e6), p)


## --- [TEMPORAL] Distribution of samples over time
    # Archean
    t = @. 4000 > simmet_t[:,agec] >= 2500;
    c, n = bincounts(simmet_t[:,SiO2c][t], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Archean Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 4000 > simmetaign_t[:,agec] >= 2500;
    c, n = bincounts(simmetaign_t[:,SiO2c][t], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h2 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Archean Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    # Proterozoic
    t = @. 2500 > simmet_t[:,agec] >= 541;
    c, n = bincounts(simmet_t[:,SiO2c][t], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h3 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Proterozoic Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 2500 > simmetaign_t[:,agec] >= 541;
    c, n = bincounts(simmetaign_t[:,SiO2c][t], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h4 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Proterozoic Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    # Phanerozoic
    t = @. 541 > simmet_t[:,agec] >= 0;
    c, n = bincounts(simmet_t[:,SiO2c][t], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h5 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.met, linecolor=colors.met,
        label="Phanerozoic Metamorphic", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    t = @. 541 > simmetaign_t[:,agec] >= 0;
    c, n = bincounts(simmetaign_t[:,SiO2c][t], SiO2min, SiO2max, nbins)
    n = float(n) ./ nansum(float(n) .* step(c))
    h6 = Plots.plot(c, n, seriestype=:bar, 
        color=colors.metaign, linecolor=colors.metaign,
        label="Phanerozoic Metaigneous", barwidths = ((SiO2max-SiO2min)/nbins),
        ylims=(0,maximum(n)+0.01)
    )

    # All together
    h = Plots.plot(h1, h3, h5, h2, h4, h6, layout=(2, 3), size=(1800,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        legend=:topleft, fg_color_legend=:white,
    )


## --- [TEMPORAL] Alternate visualization of distribution of samples over time
    # Metamorphic
    archean = @. 4000 > simmet_t[:,agec] >= 2500;
    proterozoic = @. 2500 > simmet_t[:,agec] >= 541;
    phanerozoic = @. 541 > simmet_t[:,agec] >= 0;

    h1 = Plots.plot(simmet_t[:,SiO2c][archean], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins/2),
        linecolor=:red, linewidth=3, 
        label="Archean", barwidths = ((SiO2max-SiO2min)/nbins),
    )
    Plots.plot!(simmet_t[:,SiO2c][proterozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins/2),
        linecolor=:darkorange, linewidth=3, 
        label="Proterozoic", barwidths = ((SiO2max-SiO2min)/nbins),
    )
    Plots.plot!(simmet_t[:,SiO2c][phanerozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins/2),
        linecolor=:forestgreen, linewidth=3, 
        label="Phanerozoic", barwidths = ((SiO2max-SiO2min)/nbins),
        legendtitle="Metamorphic"
    )

    # Metaigneous
    archean = @. 4000 > simmetaign_t[:,agec] >= 2500;
    proterozoic = @. 2500 > simmetaign_t[:,agec] >= 541;
    phanerozoic = @. 541 > simmetaign_t[:,agec] >= 0;

    h2 = Plots.plot(simmetaign_t[:,SiO2c][archean], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins/2),
        linecolor=:red, linewidth=3, 
        label="Archean", barwidths = ((SiO2max-SiO2min)/nbins),
    )
    Plots.plot!(simmetaign_t[:,SiO2c][proterozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins/2),
        linecolor=:darkorange, linewidth=3, 
        label="Proterozoic", barwidths = ((SiO2max-SiO2min)/nbins),
    )
    Plots.plot!(simmetaign_t[:,SiO2c][phanerozoic], 
        seriestype=:stephist, normalize=:pdf, nbins=Int(nbins/2),
        linecolor=:forestgreen, linewidth=3, 
        label="Phanerozoic", barwidths = ((SiO2max-SiO2min)/nbins),
        legendtitle="Metaigneous"
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        legend=:topleft, fg_color_legend=:white,
    )

## --- [TEMPORAL] Silica content timeseries
    c,m,e = binmeans(simmet_t[:,agec], simmet_t[:,SiO2c], 0,3800, 38)
    h1 = Plots.plot(c, m, yerror=e, label="Resampled Metamorphic", 
        markershape=:circle,
        color=colors.met, linecolor=colors.met, msc=colors.met, 
    )

    c,m,e = binmeans(simmetaign_t[:,agec], simmetaign_t[:,SiO2c], 0,3800, 38)
    h2 = Plots.plot(c, m, yerror=e, label="Resampled Metaigneous", 
        markershape=:circle, 
        color=colors.metaign, linecolor=colors.metaign, msc=colors.metaign, 
    )

    h = Plots.plot(h1, h2, layout=(2,1), size=(600,800),
        framestyle=:box, left_margin=(30,:px), ylabel="Weight", xlabel="SiO₂ [wt.%]",
        fg_color_legend=:white,
    )


## --- Resampled matched sample density 
    # For each age bin, get and normalize a histogram of silica content.
    xmin, xmax, xbins = 40, 80, 80
    xedges = xmin:(xmax-xmin)/xbins:xmax

    ymin, ymax, ybins = 0, 3800, 38*2
    yedges = ymin:(ymax-ymin)/ybins:ymax

    out = zeros(ybins, xbins)
    for i = 1:ybins
        # Filter for samples in this age bin
        t = @. yedges[i] <= simmet_t[:,agec] < yedges[i+1]

        # Count and normalize distribution of silica
        c, n = bincounts(simmet_t[:,SiO2c][t], xmin, xmax, xbins)
        # n = float(n) ./ nansum(float(n) .* step(c))
        n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))

        # Put in array for that age bin
        out[i,:] .= n
    end

    zeronan!(out)
    heatmap(out, xlabel="SiO₂ [wt.%]", ylabel="Age [Ma]", framestyle=:box,
        colorbar_title="Relative Sample Density", size=(700,400),
        color=c_gradient, left_margin=(25,:px), bottom_margin=(15,:px))
    x = 5:20:xmax-5
    y = 0:10:ybins
    z = 
    xticks!(x, string.(xmin+5:10:xmax-5))
    yticks!(y, string.(collect(0:5:38)*100))
    

## --- Compare to prior EarthChem sample density
    # Load data
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    close(fid)

    bulk_cats.met .|= bulk_cats.metaign

    out = zeros(ybins, xbins)
    for i = 1:ybins
        # Filter for samples in this age bin
        t = @. yedges[i] <= bulk.Age[bulk_cats.met] < yedges[i+1]

        # Count and normalize distribution of silica
        c, n = bincounts(bulk.SiO2[bulk_cats.met][t], xmin, xmax, xbins)
        n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))

        # Put in array for that age bin
        out[i,:] .= n
    end

    zeronan!(out)
    heatmap(out, xlabel="SiO₂ [wt.%]", ylabel="Age [Ma]", framestyle=:box,
        colorbar_title="Relative Sample Density", size=(700,400),
        color=c_gradient, left_margin=(25,:px), bottom_margin=(15,:px))
    x = 5:20:xmax-5
    y = 0:10:ybins
    z = 
    xticks!(x, string.(xmin+5:10:xmax-5))
    yticks!(y, string.(collect(0:5:38)*100))


## --- End of file 