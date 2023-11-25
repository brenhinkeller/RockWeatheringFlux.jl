# Make figures for AGU poster 11-15 Dec. 2023
# This will replicate some code that already exists, but I'll be able to run it through
# with no other files involved

## --- Set up
    using StatGeochem
    using HDF5
    using Plots

    # Local utilities and definitions
    using Measurements, LoopVectorization, Static
    include("../utilities/Utilities.jl")
    filepath = "results/figures/AGU"


## --- Slope vs. erosion rate 
    using Isoplot
    include("../utilities/yorkfit.jl")

    octopusdata = importdataset("output/octopusdata.tsv",'\t', importas=:Tuple)
    basin_srtm = importdataset("output/basin_srtm15plus_avg_maxslope.tsv", '\t', importas=:Tuple)

    # Collect erosion rate and slope data for Be and Al data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(basin_srtm.err)
    t_be = t .& .!isnan.(octopusdata.ebe_mmkyr) .& .!isnan.(octopusdata.ebe_err)
    t_al = t .& .!isnan.(octopusdata.eal_mmkyr) .& .!isnan.(octopusdata.eal_err)

    x = (
        v = [basin_srtm.avg_slope[t_be]; basin_srtm.avg_slope[t_al]],
        e = [basin_srtm.err[t_be]; basin_srtm.err[t_al]]
    )
    y = (
        v = [octopusdata.ebe_mmkyr[t_be]; octopusdata.eal_mmkyr[t_al]],
        e = [octopusdata.ebe_err[t_be]; octopusdata.eal_err[t_al]]
    )

    # Get mean in regular space for bins with equal numbers of points
    c, m, ex, ey, ey_bound = binmeans_percentile(x.v, y.v, step=5, bounderror=true)

    # Update percentiles to work for plotting
    l, u = untupleify(ey_bound)
    l = m .- l      # Lower bound
    u = u .- m      # Upper bound
    ey_bound = (l, m)

    # Fit slope to means
    fobj = yorkfit(collect(c), ex, log.(m), log.(ey))
    emmkyr(slp) = exp(slp * (fobj.slope) + (fobj.intercept))
    model, = unmeasurementify(emmkyr.(1:600))
    
    # Build plot
    h = Plots.plot(xlabel="SRTM15+ Slope [m/km]", ylabel="Erosion rate [mm/kyr]",
        framestyle=:box, legend=:bottomright, fg_color_legend=:white, 
        grid=false,
        fontfamily=:Helvetica,
        yscale=:log10,
        ylims=(10^-1, 10^4.5),
        yticks=10.0.^(0:2:4),
    )
    Plots.scatter!(basin_srtm.avg_slope,octopusdata.ebe_mmkyr, label="Be-10", 
        msc=:auto, color=:darkseagreen, 
        markersize = 2,
    )
    Plots.scatter!(h, basin_srtm.avg_slope,octopusdata.eal_mmkyr, label="Al-26", 
        msc=:auto, color=:royalblue, 
        markersize = 2,
    )
    Plots.scatter!(h, c, m, yerror=ey*2, label="",
        msc=:auto, linecolor=:black, linewidth=2, markersize=0,
    )
    Plots.scatter!(h, c, m, label="Binned Means ± 2 SEM", 
        msc=:auto, color=:black,
    )
    Plots.plot!(h, 1:length(model), model, label="exp(slope * 0.0098) + 2.97)", 
        color=:black, linewidth=3
    )
    
    display(h)
    savefig("$filepath/erosionslope.pdf")


## --- Igneous silica distributions 
    using MAT 
    using DelimitedFiles
    nsims = Int(1e7)
    SiO₂_error = 1.0

    # Volcanic (spatial)
    volcanic = matread("data/volcanicplutonic/volcanic.mat")["volcanic"];
    volcanic = NamedTuple{Tuple(Symbol.(keys(volcanic)))}(values(volcanic));
    tᵥ = @. !isnan(volcanic.Latitude) & !isnan(volcanic.Longitude) & (volcanic.Elevation .> -140)
    k = invweight_location(volcanic.Latitude[tᵥ], volcanic.Longitude[tᵥ])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = volcanic.SiO2[tᵥ]
    uncertainty = fill(SiO₂_error, length(data))
    simvolcanic = bsresample(data, uncertainty, nsims, p)
    
    # Plutonic (spatial)
    plutonic = matread("data/volcanicplutonic/plutonic.mat")["plutonic"];
    plutonic = NamedTuple{Tuple(Symbol.(keys(plutonic)))}(values(plutonic));
    tₚ = @. !isnan(plutonic.Latitude) & !isnan(plutonic.Longitude) & (plutonic.Elevation .> -140)
    k = invweight_location(plutonic.Latitude[tₚ], plutonic.Longitude[tₚ])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = plutonic.SiO2[tₚ]
    uncertainty = fill(SiO₂_error, length(data))
    simplutonic = bsresample(data, uncertainty, nsims, p)
    
    # All igneous (spatiotemporal; spatial commented out)
    tₚ .&= .!isnan.(plutonic.Age)
    tᵥ .&= .!isnan.(volcanic.Age)
    k = invweight(
        [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]], 
        [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]], 
        [plutonic.Age[tₚ]; volcanic.Age[tᵥ]]
    )
    # k = invweight_location([plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]], 
    #     [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]]
    # )
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [plutonic.SiO2[tₚ]; volcanic.SiO2[tᵥ]]
    uncertainty = fill(SiO₂_error, length(data))
    simigneous = bsresample(data, uncertainty, nsims, p)

    # Load matched EarthChem data
    fid = readdlm(matchedbulk_io)
    matches = Int.(vec(fid[:,1]))
    t = @. matches != 0

    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][matches[t]] 
        for i in eachindex(header)])
    close(fid)

    # Set up rock types to be inclusive of all subtypes
    fid = h5open("$macrostrat_io", "r")
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    sample_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorvolc
        sample_cats.volc .|= sample_cats[type]
    end
    for type in minorplut
        sample_cats.plut .|= sample_cats[type]
    end
    for type in minorign
        sample_cats.ign .|= sample_cats[type]
    end

    # Build plot
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    fig_types = (:volc, :plut, :ign)
    fig_names = ("Volcanic", "Plutonic", "Igneous")
    simVP = [simvolcanic, simplutonic, simigneous]
    for i in eachindex(fig)
        h = Plots.plot(
            framestyle=:box, grid = false,
            fontfamily=:Helvetica, 
            # title=fig_names[i],
            xlims=(40,80),
            xticks=false,
            yticks=false
        )

        c, n = bincounts(mbulk.SiO2[sample_cats[fig_types[i]]], 40, 80, 40)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, 
            seriestype=:bar, barwidths=1,
            color=colors[fig_types[i]], linecolor=:match,
            # label="Matched samples",
            label=""
        )

        c, n = bincounts(simVP[i], 40, 80, 40)
        n₂ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₂, 
            seriestype=:path, linewidth=4,
            color=:black, linecolor=:black,
            # label="Keller et al., 2015",
            label=""
        )

        Plots.ylims!(0, round(maximum([n₁; n₂]), digits=2)+0.01)
        # Plots.annotate!(42, ylims()[2]-0.005, fig_names[i])
        Plots.annotate!(((0.1, 0.92), (fig_names[i], 20)))
        fig[i] = h
    end

    # Make a common legend
    Plots.plot!(fig[1], legendfontsize = 15, fg_color_legend=:white, legend=:topright)
    Plots.plot!(fig[1], [0],[0], color=colors.volc, linecolor=:match, seriestype=:bar, 
        label=" ",)
    Plots.plot!(fig[1], [0],[0], color=colors.plut, linecolor=:match, seriestype=:bar, 
        label=" This study")
    Plots.plot!(fig[1], [0],[0], color=colors.ign, linecolor=:match, seriestype=:bar, 
        label=" ")
    Plots.plot!(fig[1], [0],[0], color=:white, linecolor=:match, label=" ")
    Plots.plot!(fig[1], [0], [0], linewidth=5, color=:black, label=" Keller et al., 2015")
    
    # Set axes
    ylabel!(fig[2], "Abundance")
    xticks!(fig[3], (40:10:80, string.(40:10:80)))
    xlabel!(fig[3], "SiO₂ [wt.%]")

    # Assemble subplots
    h = Plots.plot(fig..., layout=(3, 1), size=(900, 1800), 
        # legend=false,
        left_margin=(45,:px), right_margin=(25,:px), # bottom_margin=(15,:px),
        tickfontsize=15,
        # titleloc=:center, titlefont = font(18),
        labelfontsize=20
    )
    display(h)
    savefig(h, "$filepath/ignsilica.pdf")
    

## --- P content in igneous rocks over time 
    # Load EarthChem data
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    # Condense igneous rocks
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorvolc
        bulk_cats.volc .|= bulk_cats[type]
    end
    for type in minorplut
        bulk_cats.plut .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end

    # Get age uncertainty. If unknonwn, set at 5%
    age_uncert = Array{Float64}(undef, length(bulk.Age))
    age_uncert .= (bulk.Age_Max .- bulk.Age_Min)/2
    t = isnan.(age_uncert)
    age_uncert[t] .= bulk.Age[t] .* 0.05

    # Resample mafic and felsic silica ranges separately to distinguish mantle melting and
    # subsequent evolution and fractionation (after Keller and Schoene, 2012)
    SiO₂_range = ((43,51),(51,62),(62,74),(74,80)) 
    nsims = Int(1e7)
    # P₂O₅_error = 0.02       # From volcanic.mat (Keller et al., 2015)
    P₂O₅_error = 1.0
    simout = [Array{Float64}(undef, nsims, 2) for _ in 1:length(SiO₂_range)]

    for i in eachindex(SiO₂_range)
        # Restrict data and calculate resampling weights
        t = SiO₂_range[i][1] .< bulk.SiO2 .<= SiO₂_range[i][2]
        t .&= bulk_cats.ign
        t .&= .!isnan.(bulk.Latitude) .& .!isnan.(bulk.Longitude) .& .!isnan.(bulk.Age)
        k = invweight(bulk.Latitude[t], bulk.Longitude[t], bulk.Age[t])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

        # Resample
        data = [bulk.P2O5[t] bulk.Age[t]]
        uncertainty = [fill(P₂O₅_error, length(bulk.SiO2[t])) age_uncert[t]]
        simout[i] .= bsresample(data, uncertainty, nsims, p)
        # c, m, e = bin_bsr_means(bulk.Age[t], bulk.P2O5[t], 0, 3800, 38, 
        #     x_sigma = bulk.Age[t], y_sigma = fill(P₂O₅_error, length(bulk.SiO2[t])), 
        #     p=p, nresamplings=nsims
        # )
    end

    # Plot!
    labels=["Mafic", "Intermediate", "Felsic", "High-Silica"]
    colors=[:darkred, :red, :deeppink, :darkorange]
    h = plot(
        xlabel="Age [Ma]",
        ylabel="P₂O₅ [wt.%]",
        grid=false,
        framestyle=:box,
        fg_color_legend=:white, legendfontsize = 12,
        fontfamily=:Helvetica, 
        tickfontsize=12,
        labelfontsize=15, 
        xlims=(0,3800),
        ylims=(0,0.6),
        yticks=0:0.1:0.6,
        size=(600,600),
        bottom_margin=(15,:px), left_margin=(15,:px)
    )
    for i in eachindex(simout)
        c, m, e = binmeans(simout[i][:,2], simout[i][:,1], 0, 3800, 38)
        plot!(h, c, m, yerror=2*e, 
            linewidth=2,
            label=labels[i], 
            color=colors[i], linecolor=colors[i], msc=colors[i],
        )
    end
    display(h)
    savefig(h, "$filepath/ign_phos_timeseries.pdf")


## --- End of file 