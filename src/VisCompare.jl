## --- Compare resampled and matched silica distributions
    # Igneous silica distributions (resampled from bulk, and matched) are compared to the
    # resampled distributions from Keller et al., 2015

## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using MAT
    using Plots

    # Local utilities
    using Static
    using LoopVectorization
    include("utilities/Utilities.jl")

    # Parameters for VP resampling
    SiO2_err = 1.0
    nsims = Int(1e7)


## --- Run Monte Carlo simulations for VP Plutonic data
    plutonic = matread("data/volcanicplutonic/plutonic.mat")["plutonic"];
    plutonic = NamedTuple{Tuple(Symbol.(keys(plutonic)))}(values(plutonic));

    # Get resampling weights
    t = @. !isnan(plutonic.Latitude) & !isnan(plutonic.Longitude) & (plutonic.Elevation .> -30)
    k = invweight_location(plutonic.Latitude[t], plutonic.Longitude[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Get data and uncertainty
    data = plutonic.SiO2[t]
    uncertainty = fill(SiO2_err, length(data))

    # Run simulation
    simplutonic = bsresample(data, uncertainty, nsims, p)


## --- Run Monte Carlo simulations for VP Volcanic data
    volcanic = matread("data/volcanicplutonic/volcanic.mat")["volcanic"];
    volcanic = NamedTuple{Tuple(Symbol.(keys(volcanic)))}(values(volcanic));

    # Get resampling weights
    t = @. !isnan(volcanic.Latitude) & !isnan(volcanic.Longitude) & (volcanic.Elevation .> -30)
    k = invweight_location(volcanic.Latitude[t], volcanic.Longitude[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Get data and uncertainty
    data = volcanic.SiO2[t]
    uncertainty = fill(SiO2_err, length(data))

    # Run simulation
    simvolcanic = bsresample(data, uncertainty, nsims, p)


## --- Run Monte Carlo simulations for VP Igneous data (spatial)
    # Get resampling weights
    tₚ = @. !isnan(plutonic.Latitude) & !isnan(plutonic.Longitude) & (plutonic.Elevation .> -30)
    tᵥ = @. !isnan(volcanic.Latitude) & !isnan(volcanic.Longitude) & (volcanic.Elevation .> -30)

    lat_ign = [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]]
    lon_ign = [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]]

    k = invweight_location(lat_ign, lon_ign)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Get data and uncertainty
    data = [plutonic.SiO2[t₁]; volcanic.SiO2[t₂]]
    uncertainty = fill(SiO2_err, length(data))

    # Run simulation
    simigneous1 = bsresample(data, uncertainty, nsims, p)

    c, n = bincounts(simigneous1, SiO2min, SiO2max, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=:dimgrey, linecolor=:dimgrey,
        label="", ylabel="Abundance", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01), xlims=(SiO2min, SiO2max)
    )
    display(h)


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
    data = [plutonic.SiO2[t₁]; volcanic.SiO2[t₂]]
    uncertainty = fill(SiO2_err, length(data))

    # Run simulation
    simigneous2 = bsresample(data, uncertainty, nsims, p)

    c, n = bincounts(simigneous2, SiO2min, SiO2max, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=:dimgrey, linecolor=:dimgrey,
        label="", ylabel="Abundance", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01), xlims=(SiO2min, SiO2max)
    )
    display(h)


## --- Load resampled EarthChem data
    fid = h5open("output/resampled/resampled.h5", "r")
    header = read(fid["vars"]["header"])
    i = findfirst(x -> x=="SiO2", header)
    rocktypes = keys(fid["vars"]["data"])

    bsrsilica = NamedTuple{Tuple(Symbol.(rocktypes))}([read(fid["vars"]["data"][r]["data"])[:,i]
        for r in rocktypes
    ])
    close(fid)


## --- Load matched Earthchem data 
    data = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(data[:,1]))
    t = @. bulkidx != 0

    bulktype = string.(vec(data[:,2]))
    macro_cats = match_rocktype(bulktype[t])

    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] 
        for i in eachindex(header)]
    )
    close(fid)


## --- Plot igneous data
    SiO2min, SiO2max = 40, 80
    nbins = 160
    types = [:volc, :plut]
    labels = ["volcanic", "plutonic"]
    simVP = [simvolcanic, simplutonic]

    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(types))
    for i in eachindex(fig)
        h = plot(framestyle=:box, xlabel="SiO2 [wt.%]", ylabel="Weight", xlims=(SiO2min, SiO2max))

        # Matched samples
        c, n = bincounts(mbulk.SiO2[macro_cats[types[i]]], SiO2min, SiO2max, Int(nbins/2))
        n = float(n) ./ nansum(float(n) .* step(c))
        plot!(c, n, seriestype=:bar, color=colors[types[i]], linecolor=colors[types[i]],
            label="Matched $(labels[i])", ylims=(0, round(maximum(n), digits=2)+0.01), 
        )

        # Resampled EarthChem
        c, n = bincounts(bsrsilica[types[i]], SiO2min, SiO2max, nbins)
        n = float(n) ./ nansum(float(n) .* step(c))
        plot!(c, n, seriestype=:path, color=:black, linecolor=:black, linewidth=3,
            label="Resampled EarthChem",
        )

        # Keller et al., 2015
        c, n = bincounts(simVP[i], SiO2min, SiO2max, nbins)
        n = float(n) ./ nansum(float(n) .* step(c))
        plot!(c, n, seriestype=:path, color=:mediumblue, linecolor=:mediumblue, linewidth=3,
            label="Keller et al., 2015",
        )

        fig[i] = h
    end

    h = plot(fig..., layout=(1, 2), size=(1000,400))
    display(h)
    

## --- Plot other rock types 


## --- End of file