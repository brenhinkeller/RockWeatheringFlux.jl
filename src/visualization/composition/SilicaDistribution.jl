## --- Set up 
    # Silica distributions by rock class. Compare igneous distributions to Keller et al., 
    # 2015 (10.1038/nature14584). Compare all distributions to resampled bulk dataset
    # distributions.

    # Unique packages
    using KernelDensity, MAT 

    # Load data and base packages
    include("Definitions.jl");

    # Filter NaNs from my samples 
    here = @. !isnan(mbulk.SiO2);

    # Definitions
    nsims = Int(1e7)
    SiO₂_error = 1.0        # Given error in volcanic.mat is 0.01


## --- Load and resample volcanic / plutonic data from Keller et al., 2015
    # Volcanic (spatial)
    volcanic = matread("data/volcanicplutonic/volcanic.mat")["volcanic"];
    volcanic = NamedTuple{Tuple(Symbol.(keys(volcanic)))}(values(volcanic));
    tᵥ = @. !isnan(volcanic.Latitude) & !isnan(volcanic.Longitude) & (volcanic.Elevation .> -140)
    # k = invweight_location(volcanic.Latitude[tᵥ], volcanic.Longitude[tᵥ])
    k = vec(volcanic.k)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = volcanic.SiO2[tᵥ]
    uncertainty = fill(SiO₂_error, length(data))
    simvolcanic = bsresample(data, uncertainty, nsims, p)
    
    # Plutonic (spatial)
    plutonic = matread("data/volcanicplutonic/plutonic.mat")["plutonic"];
    plutonic = NamedTuple{Tuple(Symbol.(keys(plutonic)))}(values(plutonic));
    tₚ = @. !isnan(plutonic.Latitude) & !isnan(plutonic.Longitude) & (plutonic.Elevation .> -140)
    # k = invweight_location(plutonic.Latitude[tₚ], plutonic.Longitude[tₚ])
    k = vec(plutonic.k)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = plutonic.SiO2[tₚ]
    uncertainty = fill(SiO₂_error, length(data))
    simplutonic = bsresample(data, uncertainty, nsims, p)

    # All igneous (spatial; spatiotemporal commented out)
    tₚ .&= .!isnan.(plutonic.Age)
    tᵥ .&= .!isnan.(volcanic.Age)
    k = invweight(
        [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]], 
        [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]], 
        [plutonic.Age[tₚ]; volcanic.Age[tᵥ]]
    )
    # k = invweight_location(
    #     [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]], 
    #     [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [plutonic.SiO2[tₚ]; volcanic.SiO2[tᵥ]]
    uncertainty = fill(SiO₂_error, length(data))
    simigneous = bsresample(data, uncertainty, nsims, p)


## --- Resample (spatial) all silica distributions 
    # Preallocate 
    simout = NamedTuple{keys(bulk_cats)}(Array{Float64}(undef, nsims) for _ in keys(bulk_cats))

    # Restrict to samples with data and resample 
    t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude);
    for key in keys(match_cats) 
        s = t .& bulk_cats[key]
        k = invweight_location(bulk.Latitude[s], bulk.Longitude[s])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        data = bulk.SiO2[s]
        uncertainty = fill(SiO₂_error, count(s))
        simout[key] .= bsresample(data, uncertainty, nsims, p)
    end


## --- Igneous classes
    # Preallocate
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    fig_types = (:volc, :plut, :ign)
    fig_names = ("A. Volcanic", "B. Plutonic", "C. Igneous")
    simVP = [simvolcanic, simplutonic, simigneous]

    # Build plots
    for i in eachindex(fig)
        h = plot(
            framestyle=:box, 
            grid = false,
            fontfamily=:Helvetica, 
            xlims=(40,80),
            xticks=(40:10:80, string.(40:10:80)),
            yticks=false,
            fg_color_legend=:white,
            labelfontsize=18, titlefont=20, tickfontsize=16,legendfontsize=18,
        )

        # Raw data
        c, n = bincounts(mbulk.SiO2[match_cats[fig_types[i]]], 40, 80, 80)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, 
            seriestype=:bar, barwidths=0.5,
            color=colors[fig_types[i]], linecolor=:match, alpha=0.25,
            # label="Matched samples",
            label=""
        )

        # Kernel density estimate 
        u = kde(mbulk.SiO2[match_cats[fig_types[i]] .& here])
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[fig_types[i]], linecolor=colors[fig_types[i]],
            # label="Kernel density estimate",
            label=""
        )

        # Keller et al., 2015
        c, n = bincounts(simVP[i], 40, 80, 80)
        n₂ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₂, 
            seriestype=:path, linewidth=2,
            color=:black, linecolor=:black,
            linestyle=:dot,
            # label="Keller et al., 2015",
            label=""
        )

        # Final formatting
        Plots.ylims!(0, round(maximum([n₁; n₂; u.density]), digits=2)+0.02)
        npoints = count(match_cats[fig_types[i]])
        Plots.annotate!(((0.03, 0.97), ("n = $npoints", labelfontsize, :left, :top)))
        fig[i] = h
    end

    # Axis labels
    # ylabel!(fig[1], "Relative Abundance")
    # xlabel!(fig[2], "SiO2 [wt.%]")

    # Shared legend
    Plots.plot!(fig[1], [0],[0], color=:white, linecolor=:match, label=" ")
    Plots.plot!(fig[1], [0],[0], color=colors.volc, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" This study",)
    Plots.plot!(fig[1], [0], [0], linewidth=2, color=:black, linestyle=:dot,
        label=" Keller et al., 2015")

    # Assemble plots
    h = Plots.plot(fig..., layout=(1, 3), size=(2000, 500), 
        left_margin=(25,:px), right_margin=(25,:px), bottom_margin=(55,:px),
        labelfontsize=18, titlefont=20, tickfontsize=16,
        legendfontsize=18, fg_color_legend=:white, legend=:topright,
        xlabel="SiO2 [wt.%]"
    )
    display(h)
    savefig(h, "$filepath/histogram_ign.pdf")

    savefig(fig[1], "$filepath/histogram_volc.pdf")
    savefig(fig[2], "$filepath/histogram_plut.pdf")
    savefig(fig[3], "$filepath/histogram_ign_isolate.pdf")


## --- All 
    # Preallocate 
    fig_types = collect(keys(match_cats))
    fig_names =  ("A. Siliciclastic", "B. Shale", "C. Carbonate", "D. Evaporite", 
        "E. Chert", "F. Phosphorite", "G. Coal", "H. Sedimentary", "I. Komatiite", 
        "J. Basalt", "K. Andesite", "L. Dacite", "M. Rhyolite", "N. Alkaline Volcanic", 
        "O. Volcaniclastic", "P. Volcanic", "Q. Peridotite", "R. Pyroxenite", "S. Gabbro", 
        "T. Diorite", "U. Trondhjemite", "V. Tonalite", "W. Granodiorite", "X. Granite", 
        "Y. Alkaline Plutonic", "Z. Plutonic", "AA. Carbonatite", "AB. Igneous", 
        "AC. All Metamorphic"
    )
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(fig_names) + 1)

    # Build plots 
    for i in eachindex(fig_types)
        h = plot(
            framestyle=:box, 
            grid = false,
            fontfamily=:Helvetica, 
            xlims=(0,100),
            xticks=(0:20:100, string.(0:20:100)),
            yticks=false
        )

        # Make sure there's actually data
        if count(match_cats[fig_types[i]]) == 0
            Plots.annotate!(((0.03, 0.97), (fig_names[i] * "\nNo Data", 18, :left, :top)))
            fig[i] = h
            continue
        end

        # Raw data
        c, n = bincounts(mbulk.SiO2[match_cats[fig_types[i]]], 0, 100, 100)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, 
            seriestype=:bar, barwidths=0.5,
            color=colors[fig_types[i]], linecolor=:match, alpha=0.25,
            # label="Matched samples",
            label=""
        )

        # Kernel density estimate 
        u = kde(mbulk.SiO2[match_cats[fig_types[i]] .& here])
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[fig_types[i]], linecolor=colors[fig_types[i]],
            # label="Kernel density estimate",
            label=""
        )

        # Resampled distribution
        c, n = bincounts(simout[fig_types[i]], 0, 100, 100)
        n₂ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₂, 
            seriestype=:path, linewidth=2,
            color=:black, linecolor=:black,
            linestyle=:dot,
            label=""
        )

        # Final formatting
        Plots.ylims!(0, round(maximum([n₁; n₂; u.density]), digits=2)+0.01)
        npoints = count(match_cats[fig_types[i]])
        Plots.annotate!(((0.03, 0.97), (fig_names[i] * "\nn = $npoints", 18, :left, :top)))
        fig[i] = h
    end
    
    # Axis labels
    ylabel!(fig[11], "Relative Abundance")
    xlabel!(fig[28], "SiO2 [wt.%]")

    # Common legend, as it's own plot
    h = Plots.plot(
        framestyle=:none, grid = false,
        fontfamily=:Helvetica,
        xlims = (1, 10), ylims = (1, 10),
        xticks=false, yticks=false
    )
    h = Plots.plot!(h, legendfontsize = 24, fg_color_legend=:white, legend=:inside)
    Plots.plot!(h, [0],[0], color=colors.evap, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" Matched Samples")

    Plots.plot!(h, [0],[0], color=:white, linecolor=:match, label=" ")
    Plots.plot!(h, [0],[0], color=colors.evap, linewidth=5, 
        label=" Kernel Density Estimate")

    Plots.plot!(h, [0],[0], color=:white, linecolor=:match, label=" ")
    Plots.plot!(h, [0],[0], color=:black, linestyle=:dot, linewidth=5,
        label=" Spatially Resampled\n Geochemical Data")
    fig[30] = h

    # Assemble plots
    # Size: 500px for each row, 600 px for each column
    h = Plots.plot(fig..., layout=(6, 5), size=(3000, 3000), 
        # legend=false,
        left_margin=(75,:px), right_margin=(25,:px), bottom_margin=(45,:px),
        tickfontsize=12,
        # titleloc=:center, titlefont = font(18),
        labelfontsize=24
    )
    display(h)
    savefig(h, "$filepath/histogram_all_classes.pdf")

    
## --- End of file