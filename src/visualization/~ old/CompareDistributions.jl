## --- Compare resampled and matched silica distributions
    # Igneous silica distributions (resampled from bulk, and matched) are compared to the
    # resampled distributions from Keller et al., 2015

## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using MAT
    using DelimitedFiles

    using Plots
    using StatsPlots

    # Local utilities
    using Static
    using LoopVectorization
    using Measurements

    include("../utilities/Utilities.jl")

    # Parameters for VP resampling
    SiO2_err = 1.0
    nsims = Int(1e7)


## --- Run Monte Carlo simulations for VP Plutonic data
    plutonic = matread("data/volcanicplutonic/plutonic.mat")["plutonic"];
    plutonic = NamedTuple{Tuple(Symbol.(keys(plutonic)))}(values(plutonic));

    # Get resampling weights
    tₚ = @. !isnan(plutonic.Latitude) & !isnan(plutonic.Longitude) & (plutonic.Elevation .> -30)
    k = invweight_location(plutonic.Latitude[tₚ], plutonic.Longitude[tₚ])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Get data and uncertainty
    data = plutonic.SiO2[tₚ]
    uncertainty = fill(SiO2_err, length(data))

    # Run simulation
    simplutonic = bsresample(data, uncertainty, nsims, p)


## --- Run Monte Carlo simulations for VP Volcanic data
    volcanic = matread("data/volcanicplutonic/volcanic.mat")["volcanic"];
    volcanic = NamedTuple{Tuple(Symbol.(keys(volcanic)))}(values(volcanic));

    # Get resampling weights
    tᵥ = @. !isnan(volcanic.Latitude) & !isnan(volcanic.Longitude) & (volcanic.Elevation .> -30)
    k = invweight_location(volcanic.Latitude[tᵥ], volcanic.Longitude[tᵥ])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Get data and uncertainty
    data = volcanic.SiO2[tᵥ]
    uncertainty = fill(SiO2_err, length(data))

    # Run simulation
    simvolcanic = bsresample(data, uncertainty, nsims, p)


## --- Run Monte Carlo simulations for VP Igneous data (spatiotemporal)
    # Get resampling weights
    tₚ .&= .!isnan.(plutonic.Age)
    tᵥ .&= .!isnan.(volcanic.Age)

    lat_ign = [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]]
    lon_ign = [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]]
    age_ign = [plutonic.Age[tₚ]; volcanic.Age[tᵥ]]

    k = invweight(lat_ign, lon_ign, age_ign)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Get data and uncertainty
    data = [plutonic.SiO2[tₚ]; volcanic.SiO2[tᵥ]]
    uncertainty = fill(SiO2_err, length(data))

    # Run simulation
    simigneous = bsresample(data, uncertainty, nsims, p)


## --- Load resampled EarthChem data
    fid = h5open("output/resampled/resampled.h5", "r")
    header = read(fid["vars"]["header"])
    i = findfirst(x -> x=="SiO2", header)
    rocktypes = keys(fid["vars"]["data"])

    bsrsilica = NamedTuple{Tuple(Symbol.(rocktypes))}(
        [read(fid["vars"]["data"][r]["data"])[:,i] for r in rocktypes])
    close(fid)


## --- Load matched Earthchem data 
    fid = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    bulktype = string.(vec(fid[:,2]))
    macro_cats = match_rocktype(bulktype[t])    # Majors inclusive

    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] 
        for i in eachindex(header)]
    )
    close(fid)


## --- Load unmatched EarthChem data
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    close(fid)

    # Major types include minor types
    minorsed, minorign, minormet = get_minor_types()
    for type in minorsed
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end
    for type in minormet
        bulk_cats.met .|= bulk_cats[type]
    end


## --- Plot igneous data
    SiO2min, SiO2max = 40, 80
    nbins = 160
    nbins_matched = Int(nbins/2)
    types = [:volc, :plut, :ign]
    labels = ["volcanic", "plutonic", "all igneous"]
    simVP = [simvolcanic, simplutonic, simigneous]

    # Volcanic, Plutonic
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(types))
    for i in eachindex(fig)
        h = Plots.plot(framestyle=:box, # xlabel="SiO2 [wt.%]", ylabel="Weight", 
            xlims=(SiO2min, SiO2max), left_margin=(30, :px), bottom_margin=(30, :px), 
        )

        i==1 && ylabel!("Abundance")
        i==2 && xlabel!("SiO2 [wt.%]")

        # Matched samples
        c, n = bincounts(mbulk.SiO2[macro_cats[types[i]]], SiO2min, SiO2max, nbins_matched)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, seriestype=:bar, color=colors[types[i]], linecolor=:match,
            label="Matched $(labels[i])", barwidths = ((SiO2max-SiO2min)/nbins_matched)
        )

        # Resampled EarthChem
        c, n = bincounts(bsrsilica[types[i]], SiO2min, SiO2max, nbins)
        n₂ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₂, seriestype=:path, color=:blue, linecolor=:blue, linewidth=2,
            label="Resampled EarthChem",
        )

        # Prior Earthchem
        c, n = bincounts(bulk.SiO2[bulk_cats[types[i]]], SiO2min, SiO2max, Int(nbins/2))
        n₃ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₃, seriestype=:path, color=:red, linecolor=:red, linewidth=2,
            label="EarthChem prior",
        )

        # Keller et al., 2015
        c, n = bincounts(simVP[i], SiO2min, SiO2max, nbins)
        n₄ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₄, seriestype=:path, color=:black, linecolor=:black, linewidth=4,
            label="Keller et al., 2015",
        )

        # Plots.ylims!(0, round(maximum([n₁; n₂; n₃; n₄]), digits=2)+0.01)
        Plots.ylims!(0, round(maximum([n₁; n₄]), digits=2)+0.01)
        fig[i] = h
    end

    # h = Plots.plot(fig..., layout=(2, 2), size=(1000,800), legend=false)
    h = Plots.plot(fig..., layout=(1, 3), size=(1800,500), legend=false,
        left_margin=(35,:px), bottom_margin=(35,:px), top_margin=(35,:px),
        title =["Volcanic" "Plutonic" "All Igneous"], titleloc=:left, titlefont = font(15)
    )    
    display(h)
    # savefig(h, "results/figures/dist_ign.png")
    savefig(h, "results/figures/dist_ign.pdf")


## --- Plot other rock types 
    SiO2min, SiO2max = 0, 100
    nbins = 200
    nbins_matched = Int(nbins/2)
    types = [:metased, :metaign, :met, :siliciclast, :shale, :carb, :chert, :sed]
    labels = string.(types)

    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(types))
    for i in eachindex(fig)
        h = Plots.plot(framestyle=:box, xlabel="SiO2 [wt.%]", ylabel="Weight", 
            xlims=(SiO2min, SiO2max), left_margin=(30, :px), bottom_margin=(30, :px),
            legend=:topleft    
        )

        # Matched samples
        c, n = bincounts(mbulk.SiO2[macro_cats[types[i]]], SiO2min, SiO2max, nbins_matched)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, seriestype=:bar, color=colors[types[i]], linecolor=colors[types[i]],
            label="Matched $(labels[i])", barwidths = ((SiO2max-SiO2min)/nbins_matched)
        )

        # Prior Earthchem
        c, n = bincounts(bulk.SiO2[bulk_cats[types[i]]], SiO2min, SiO2max, nbins)
        n₃ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₃, seriestype=:path, color=:red, linecolor=:red, linewidth=3,
            label="EarthChem prior",
        )

        # Resampled EarthChem
        c, n = bincounts(bsrsilica[types[i]], SiO2min, SiO2max, nbins)
        n₂ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₂, seriestype=:path, color=:black, linecolor=:black, linewidth=3,
            label="Resampled EarthChem",
        )

        Plots.ylims!(0, round(maximum([n₁; n₂; n₃]), digits=2)+0.01)
        fig[i] = h
    end

    h = Plots.plot(fig..., layout=(3, 3), size=(2000,1400))
    display(h)
    savefig(h, "results/figures/dist_sedmet.png")


## --- Plot metamorphic and metaigneous rocks for thesis proposal
    h1 = Plots.plot(framestyle=:box, xlabel="SiO2 [wt.%]", ylabel="Abundance", 
        xlims=(SiO2min, SiO2max), ylims=(0, 0.1),
        left_margin=(25, :px), right_margin=(25, :px), bottom_margin=(25, :px),
        legend=:topleft, fg_color_legend=:white,)

    c, n = bincounts(mbulk.SiO2[macro_cats.metaign], SiO2min, SiO2max, nbins_matched)
    n₁ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(c, n₁, seriestype=:bar, color=colors.metaign, linecolor=colors.metaign,
        label="Matched metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched))
    
    c, n = bincounts(bulk.SiO2[bulk_cats[types[i]]], SiO2min, SiO2max, nbins)
    n₂ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(c, n₂, seriestype=:path, color=:red, linecolor=:red, linewidth=3,
        label="EarthChem prior",)
    
    c, n = bincounts(bsrsilica.metaign, SiO2min, SiO2max, nbins)
    n₃ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(c, n₃, seriestype=:path, color=:black, linecolor=:black, linewidth=3,
        label="Resampled EarthChem",
    )

    h2 = Plots.plot(framestyle=:box, xlabel="SiO2 [wt.%]", ylabel="Abundance", 
        xlims=(SiO2min, SiO2max), ylims=(0, 0.1),
        left_margin=(25, :px), right_margin=(25, :px), bottom_margin=(25, :px),
        legend=:topleft, fg_color_legend=:white,)
    
    c, n = bincounts(mbulk.SiO2[macro_cats.metaign], SiO2min, SiO2max, nbins_matched)
    n₁ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(c, n₁, seriestype=:bar, color=colors.metaign, linecolor=colors.metaign,
        label="Matched metaigneous", barwidths = ((SiO2max-SiO2min)/nbins_matched))
    c, n = bincounts(bulk.SiO2[bulk_cats[types[i]]], SiO2min, SiO2max, nbins)
    n₂ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(c, n₂, seriestype=:path, color=:red, linecolor=:red, linewidth=3,
        label="EarthChem prior",)
    c, n = bincounts(bsrsilica.metaign, SiO2min, SiO2max, nbins)
    n₃ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(c, n₃, seriestype=:path, color=:black, linecolor=:black, linewidth=3,
        label="Resampled EarthChem",)

    h = Plots.plot(h1, h2, layout=(1,2))
    display(h)
    savefig(h, "results/figures/dist_metaign.pdf")


## --- Plot CaO data for carbonates
    # Load resampled CaO
    fid = h5open("output/resampled/resampled.h5", "r")
    header = read(fid["vars"]["header"])
    i = findfirst(x -> x=="CaO", header)
    rocktypes = ["carb"]

    bsrCaO = NamedTuple{Tuple(Symbol.(rocktypes))}([read(fid["vars"]["data"][r]["data"])[:,i]
        for r in rocktypes
    ])
    close(fid)

    # Set up
    CaOmin, CaOmax = 0, 100
    nbins = 100
    nbins_matched = 100
    types = Symbol.(rocktypes)
    labels = string.(types)

    # Make plots
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(types))
    for i in eachindex(fig)
        h = plot(framestyle=:box, xlabel="CaO [wt.%]", ylabel="Weight", 
            xlims=(CaOmin, CaOmax), left_margin=(30, :px), bottom_margin=(30, :px),
            legend=:topleft    
        )

        # Matched samples
        c, n = bincounts(mbulk.CaO[macro_cats[types[i]]], CaOmin, CaOmax, nbins_matched)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, seriestype=:bar, color=colors[types[i]], linecolor=colors[types[i]],
            label="Matched $(labels[i])", barwidths = ((CaOmax-CaOmin)/nbins_matched)
        )

        # Prior Earthchem
        c, n = bincounts(bulk.CaO[bulk_cats[types[i]]], CaOmin, CaOmax, nbins)
        n₃ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₃, seriestype=:path, color=:red, linecolor=:red, linewidth=3,
            label="EarthChem prior",
        )

        # Resampled EarthChem
        c, n = bincounts(bsrCaO[types[i]], CaOmin, CaOmax, nbins)
        n₂ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₂, seriestype=:path, color=:black, linecolor=:black, linewidth=2,
            label="Resampled EarthChem",
        )

        Plots.ylims!(0, round(maximum([n₁; n₂]), digits=2)+0.01)
        fig[i] = h
    end

    h = Plots.plot(fig..., layout=(1, 2), size=(600,400))
    display(h)
    savefig(h, "results/figures/dist_cao.png")


## --- End of file