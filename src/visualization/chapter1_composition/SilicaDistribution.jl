## --- Set up 
    # Silica distributions by rock class. Compare igneous distributions to Keller et al., 
    # 2015 (10.1038/nature14584)

    # Unique packages
    using KernelDensity, MAT 

    # Load data and base packages
    include("Definitions.jl");

    # Filter NaNs from my samples 
    here = @. !isnan(mbulk.SiO2);


## --- Load and resample volcanic / plutonic data from Keller et al., 2015
    # Set up
    nsims = Int(1e7)
    SiO₂_error = 1.0

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
    # k = invweight(
    #     [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]], 
    #     [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]], 
    #     [plutonic.Age[tₚ]; volcanic.Age[tᵥ]]
    # )
    k = invweight_location(
        [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]], 
        [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [plutonic.SiO2[tₚ]; volcanic.SiO2[tᵥ]]
    uncertainty = fill(SiO₂_error, length(data))
    simigneous = bsresample(data, uncertainty, nsims, p)


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
            yticks=false
        )

        # Raw data
        c, n = bincounts(mbulk.SiO2[macro_cats[fig_types[i]]], 40, 80, 80)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, 
            seriestype=:bar, barwidths=0.5,
            color=colors[fig_types[i]], linecolor=:match, alpha=0.25,
            # label="Matched samples",
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

        # Kernel density estimate 
        u = kde(mbulk.SiO2[macro_cats[fig_types[i]] .& here])
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[fig_types[i]], linecolor=colors[fig_types[i]],
            # label="Kernel density estimate",
            label=""
        )

        # Final formatting
        Plots.ylims!(0, round(maximum([n₁; n₂; u.density]), digits=2)+0.01)
        npoints = count(macro_cats[fig_types[i]])
        Plots.annotate!(((0.03, 0.97), (fig_names[i] * "\nn = $npoints", 12, :left, :top)))
        fig[i] = h
    end

    # Axis labels
    ylabel!(fig[1], "Relative Abundance")
    xlabel!(fig[2], "SiO2 [wt.%]")

    # Shared legend
    Plots.plot!(fig[1], legendfontsize = 12, fg_color_legend=:white, legend=:topright)
    Plots.plot!(fig[1], [0],[0], color=colors.volc, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" ",)
    Plots.plot!(fig[1], [0],[0], color=colors.plut, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" This study")
    Plots.plot!(fig[1], [0],[0], color=colors.ign, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" ")

    Plots.plot!(fig[1], [0],[0], color=:white, linecolor=:match, label=" ")
    Plots.plot!(fig[1], [0],[0], color=colors.volc, linewidth=5, label=" ",)
    Plots.plot!(fig[1], [0],[0], color=colors.plut, linewidth=5, 
        label=" Kernel Density Estimate")
    Plots.plot!(fig[1], [0],[0], color=colors.ign, linewidth=5, label=" ")

    Plots.plot!(fig[1], [0],[0], color=:white, linecolor=:match, label=" ")
    Plots.plot!(fig[1], [0], [0], linewidth=2, color=:black, linestyle=:dot,
        label=" Keller et al., 2015")

    # Assemble plots
    h = Plots.plot(fig..., layout=(1, 3), size=(1800, 500), 
        # legend=false,
        left_margin=(75,:px), right_margin=(25,:px), bottom_margin=(45,:px),
        tickfontsize=12,
        # titleloc=:center, titlefont = font(18),
        labelfontsize=14
    )
    display(h)
    savefig(h, "$filepath/silica_ign.pdf")


## --- All 
    # Preallocate 
    fig_types = collect(keys(macro_cats))
    fig_names =  ("A. Siliciclastic", "B. Shale", "C. Carbonate", "D. Evaporite", 
        "E. Chert", "F. Phosphorite", "G. Coal", "H. Sedimentary", "I. Komatiite", 
        "J. Basalt", "K. Andesite", "L. Dacite", "M. Rhyolite", "N. Alkaline Volcanic", 
        "O. Volcaniclastic", "P. Volcanic", "Q. Peridotite", "R. Pyroxenite", "S. Gabbro", 
        "T. Diorite", "U. Trondhjemite", "V. Tonalite", "W. Granodiorite", "X. Granite", 
        "Y. Alkaline Plutonic", "Z. Plutonic", "AA. Carbonatite", "AB. Igneous", 
        "AC. Unspecified Metamorphic"
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

        # Raw data
        c, n = bincounts(mbulk.SiO2[macro_cats[fig_types[i]]], 0, 100, 100)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, 
            seriestype=:bar, barwidths=0.5,
            color=colors[fig_types[i]], linecolor=:match, alpha=0.25,
            # label="Matched samples",
            label=""
        )

        # Kernel density estimate 
        u = kde(mbulk.SiO2[macro_cats[fig_types[i]] .& here])
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[fig_types[i]], linecolor=colors[fig_types[i]],
            # label="Kernel density estimate",
            label=""
        )

        # Final formatting
        Plots.ylims!(0, round(maximum([n₁; u.density]), digits=2)+0.01)
        npoints = count(macro_cats[fig_types[i]])
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
    savefig(h, "$filepath/silica_all.pdf")

    

## --- More in depth look at igneous rocks 
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];
    mplut = (:plut, minorplut...)   # 10 (3x4)
    mvolc = (:volc, minorvolc...)   # 8  (3x3)
    exclude_minor!(macro_cats)

    division = (mplut, mvolc)
    for d in division
        fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(d))
        fig_names = string.(d)
        for i in eachindex(fig)
            h = plot(
                framestyle=:box, 
                grid = false,
                fontfamily=:Helvetica, 
                xlims=(40,80),
                xticks=(40:10:80, string.(40:10:80)),
                yticks=false
            )

            # Raw data
            c, n = bincounts(mbulk.SiO2[macro_cats[d[i]]], 40, 80, 80)
            n = float(n) ./ nansum(float(n) .* step(c))
            Plots.plot!(c, n, 
                seriestype=:bar, barwidths=0.5,
                color=colors[d[i]], linecolor=:match, alpha=0.25, label=""
            )

            # Kernel density estimate 
            u = kde(mbulk.SiO2[macro_cats[d[i]] .& here])
            Plots.plot!(u.x, u.density, 
                seriestype=:path, linewidth=4,
                color=colors[d[i]], linecolor=colors[d[i]], label=""
            )

            # Final formatting
            Plots.ylims!(0, round(maximum([n; u.density]), digits=2)+0.01)
            npoints = count(macro_cats[d[i]])
            Plots.annotate!(((0.03, 0.97), (fig_names[i] * "\nn = $npoints", 18, :left, :top)))
            fig[i] = h
        end
        ncols = ceil(Int, length(d)/3)
        h = Plots.plot(fig..., layout = (3,ncols), size=(600*ncols, 400*3), 
            title="Matched Samples [$geochem_fid]")
        display(h)
    end


## --- Does the same trend happen in the unmatched data?
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];
    mplut = (:plut, minorplut...)   # 10 (3x4)
    mvolc = (:volc, minorvolc...)   # 8  (3x3)
    exclude_minor!(bulk_cats)
    here = @. !isnan(bulk.SiO2);

    division = (mplut, mvolc)
    for d in division
        fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(d))
        fig_names = string.(d)
        for i in eachindex(fig)
            h = plot(
                framestyle=:box, 
                grid = false,
                fontfamily=:Helvetica, 
                xlims=(40,80),
                xticks=(40:10:80, string.(40:10:80)),
                yticks=false
            )

            # Raw data
            c, n = bincounts(bulk.SiO2[bulk_cats[d[i]]], 40, 80, 80)
            n = float(n) ./ nansum(float(n) .* step(c))
            Plots.plot!(c, n, 
                seriestype=:bar, barwidths=0.5,
                color=colors[d[i]], linecolor=:match, alpha=0.25, label=""
            )

            # Kernel density estimate 
            u = kde(bulk.SiO2[bulk_cats[d[i]] .& here])
            Plots.plot!(u.x, u.density, 
                seriestype=:path, linewidth=4,
                color=colors[d[i]], linecolor=colors[d[i]], label=""
            )

            # Final formatting
            Plots.ylims!(0, round(maximum([n; u.density]), digits=2)+0.01)
            npoints = count(bulk_cats[d[i]])
            Plots.annotate!(((0.03, 0.97), (fig_names[i] * "\nn = $npoints", 18, :left, :top)))
            fig[i] = h
        end
        ncols = ceil(Int, length(d)/3)
        h = Plots.plot(fig..., layout = (3,ncols), size=(600*ncols, 400*3), 
            title="Bulk Geochemical Data [$geochem_fid]"
        )
        display(h)
    end


## --- Silica distribution in volc / plut rocks if we exclude undifferentiated samples 
    # For igneous rocks, plot volcanic / plutonic silica distribution, where the silica
    # distribution is made from the unmatched sample distribution weighted by the relative
    # abundance of each rock class

    # Get proportional distribution on Earth's surface
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    include_minor!(macro_cats)
    include_minor!(bulk_cats)
    here = @. !isnan(mbulk.SiO2);

    ptype = (
        sed = float.([count(macro_cats[i]) for i in minorsed]),
        volc = float.([count(macro_cats[i]) for i in minorvolc]),
        plut = float.([count(macro_cats[i]) for i in minorplut]),
        ign = float.([count(macro_cats[i]) for i in minorign]),
    )
    ptype.sed ./= nansum(ptype.sed)
    ptype.volc ./= nansum(ptype.volc)
    ptype.plut ./= nansum(ptype.plut)
    ptype.ign ./= nansum(ptype.ign)

    # Make distributions
    minorclasses = (volc=minorvolc, plut=minorplut)
    for class in (:volc, :plut)
        dist = zeros(80)
        for i in eachindex(minorclasses[class])
            c, n = bincounts(bulk.SiO2[bulk_cats[minorclasses[class][i]]], 40, 80, 80)
            n = float(n) ./ nansum(float(n) .* step(c))
            dist .+= (n*ptype[class][i])
        end

        # Get matched KDE 
        u = kde(mbulk.SiO2[macro_cats[class] .& here])

        # Assemble plot
        h = plot(
            framestyle=:box, 
            grid = false,
            fontfamily=:Helvetica, 
            xlims=(40,80), ylims=(0, maximum(dist)+0.01),
            xticks=(40:10:80, string.(40:10:80)),
            yticks=false,
            xlabel="SiO₂ [wt.%]", ylabel="Relative Abundance",
            title="$class samples [$geochem_fid]",
            left_margin=(20,:px)
        )
        Plots.plot!(h, cntr(40:0.5:80), dist, 
            seriestype=:bar, barwidths=0.5,
            color=colors[class], linecolor=:match, alpha=0.5, label="Weighted Bulk"
        )
        Plots.plot!(h, u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[class], linecolor=colors[class], label="Matched KDE"
        )
        display(h)
    end


## --- OK so as above, but if we resample the distributions for each rock type.
    # The above code assumes that the silica distribution of each rock class in the
    # geochemical dataset is representative of the "true" silica distribution. 
    # 
    # To avoid making that assumption, what if we spatially resample (NOT temporal, after
    # Keller et al., 2015) each constituent sub-class before getting its distribution and 
    # computing the total volcanic / plutonic distribution

    # Set up 
    include_minor!(macro_cats)
    include_minor!(bulk_cats)
    here = @. !isnan(mbulk.SiO2);

    nsims = Int(1e7)
    SiO₂_error = 1.0

    # Proportional distribution of igneous rock types on Earth's surface
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    ptype = (
        volc = float.([count(macro_cats[i]) for i in minorvolc]),
        plut = float.([count(macro_cats[i]) for i in minorplut]),
        ign = float.([count(macro_cats[i]) for i in minorign]),
    )
    ptype.volc ./= nansum(ptype.volc)
    ptype.plut ./= nansum(ptype.plut)
    ptype.ign ./= nansum(ptype.ign)

    # Make plots
    minorclasses = (volc=minorvolc, plut=minorplut)
    for class in (:volc, :plut)
        # Compute distributions
        dist = zeros(80)
        for i in eachindex(minorclasses[class])
            # Get samples in this class
            s = bulk_cats[minorclasses[class][i]]
            
            # Resample
            k = invweight_location(bulk.Latitude[s], bulk.Longitude[s])
            p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
            data = bulk.SiO2[s]
            uncertainty = fill(SiO₂_error, length(data))
            simout = bsresample(data, uncertainty, nsims, p)

            # Get distribution of the simulation and add to the total distribution
            c, n = bincounts(simout, 40, 80, 80)
            n = float(n) ./ nansum(float(n) .* step(c))
            dist .+= (n*ptype[class][i])
        end

        # Compare resampled distribution of the entire class 
        k = invweight_location(bulk.Latitude[bulk_cats[class]], bulk.Longitude[bulk_cats[class]])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        data = bulk.SiO2[bulk_cats[class]]
        uncertainty = fill(SiO₂_error, length(data))
        simout = bsresample(data, uncertainty, nsims, p)
        c₁, n₁ = bincounts(simout, 40, 80, 80)
        n₁ = float(n₁) ./ nansum(float(n₁) .* step(c₁))

        # Get matched KDE 
        u = kde(mbulk.SiO2[macro_cats[class] .& here])

        # Assemble plot
        h = plot(
            framestyle=:box, 
            grid = false,
            fontfamily=:Helvetica, 
            xlims=(40,80), ylims=(0, maximum(dist)+0.01),
            xticks=(40:10:80, string.(40:10:80)),
            yticks=false,
            xlabel="SiO₂ [wt.%]", ylabel="Relative Abundance",
            title="$class samples [$geochem_fid]",
            left_margin=(20,:px)
        )
        Plots.plot!(h, cntr(40:0.5:80), dist, 
            seriestype=:bar, barwidths=0.5,
            color=colors[class], linecolor=:match, alpha=0.5, label="Weighted Bulk"
        )
        Plots.plot!(h, c₁, n₁, 
            seriestype=:path, linewidth=4, linestyle=:dot,
            color=:grey, linecolor=:grey, label="Resampled $class"
        )
        Plots.plot!(h, u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[class], linecolor=colors[class], label="Matched KDE"
        )
        display(h)
    end


## --- End of file