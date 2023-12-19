# Figures for the submitted paper

    # Packages
    using StatGeochem
    using HDF5
    using MAT
    using DelimitedFiles
    using KernelDensity
    using Plots
    using Colors

    # Local utilities
    using Measurements, Static
    using LoopVectorization: @turbo
    include("../utilities/Utilities.jl")
    filepath = "results/figures/CompPaper"


## --- Load data
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

    # Filter all NaNs
    here = @. !isnan(mbulk.SiO2);

    # Set up rock types to be inclusive of all subtypes
    fid = h5open("$macrostrat_io", "r")
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    sample_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorsed
        sample_cats.sed .|= sample_cats[type]
    end
    for type in minorvolc
        sample_cats.volc .|= sample_cats[type]
    end
    for type in minorplut
        sample_cats.plut .|= sample_cats[type]
    end
    for type in minorign
        sample_cats.ign .|= sample_cats[type]
    end

    # I hate cover.
    sample_cats = delete_cover(sample_cats)


## --- Apply a kernel density estimate to the matched igneous distributions
    # Load and resample Keller et al., 2015 data 
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

    # All igneous (spatiotemporal; spatial commented out)
    tₚ .&= .!isnan.(plutonic.Age)
    tᵥ .&= .!isnan.(volcanic.Age)
    k = invweight(
        [plutonic.Latitude[tₚ]; volcanic.Latitude[tᵥ]], 
        [plutonic.Longitude[tₚ]; volcanic.Longitude[tᵥ]], 
        [plutonic.Age[tₚ]; volcanic.Age[tᵥ]]
    )
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [plutonic.SiO2[tₚ]; volcanic.SiO2[tᵥ]]
    uncertainty = fill(SiO₂_error, length(data))
    simigneous = bsresample(data, uncertainty, nsims, p)

    # Plot igneous samples
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    fig_types = (:volc, :plut, :ign)
    fig_names = ("Volcanic", "Plutonic", "Igneous")
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
        c, n = bincounts(mbulk.SiO2[sample_cats[fig_types[i]]], 40, 80, 80)
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
        u = kde(mbulk.SiO2[sample_cats[fig_types[i]] .& here])
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[fig_types[i]], linecolor=colors[fig_types[i]],
            # label="Kernel density estimate",
            label=""
        )

        # Final formatting
        Plots.ylims!(0, round(maximum([n₁; n₂; u.density]), digits=2)+0.01)
        npoints = count(sample_cats[fig_types[i]])
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
    
    
## --- All rock types (supplemental figure)
    # Set up. 
    # Some rock types removed if silica content has no pattern and sample density is low
    fig_types = collect(keys(sample_cats))
    # deleteat!(fig_types, findall(x->x==:coal,fig_types))            # Poor silica content
    # deleteat!(fig_types, findall(x->x==:phosphorite,fig_types))     # Poor silica content
    # deleteat!(fig_types, findall(x->x==:carbonatite,fig_types))     # Low sample density

    # fig_names =  ("Siliciclastic", "Shale", "Carbonate", "Evaporite", "Chert",
    #     "Sedimentary", "Komatiite", "Basalt", "Andesite", "Dacite", 
    #     "Rhyolite", "Alkaline Volcanic", "Volcaniclastic", "Volcanic", "Peridotite", 
    #     "Pyroxenite", "Gabbro", "Diorite", "Trondhjemite", "Tonalite", "Granodiorite", 
    #     "Granite", "Alkaline Plutonic", "Plutonic", "Igneous", "Metamorphic"
    # )
    fig_names =  ("Siliciclastic", "Shale", "Carbonate", "Evaporite", "Chert", 
        "Phosphorite", "Coal", "Sedimentary", "Komatiite", "Basalt", "Andesite", "Dacite", 
        "Rhyolite", "Alkaline Volcanic", "Volcaniclastic", "Volcanic", "Peridotite", 
        "Pyroxenite", "Gabbro", "Diorite", "Trondhjemite", "Tonalite", "Granodiorite", 
        "Granite", "Alkaline Plutonic", "Plutonic", "Carbonatite", "Igneous", 
        "Unspecified Metamorphic"
    )
    
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(fig_names) + 1)

    # Note that some rock types, particularly carbonate and chert, will have apparently
    # odd silica distribution. This is because these rock types often occur together
    # (e.g. chert nodules in silica) but are mapped as a single unit. The anomalous 
    # distributions are the impact of this "cross-contamination"

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
        c, n = bincounts(mbulk.SiO2[sample_cats[fig_types[i]]], 0, 100, 100)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, 
            seriestype=:bar, barwidths=0.5,
            color=colors[fig_types[i]], linecolor=:match, alpha=0.25,
            # label="Matched samples",
            label=""
        )

        # Kernel density estimate 
        u = kde(mbulk.SiO2[sample_cats[fig_types[i]] .& here])
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[fig_types[i]], linecolor=colors[fig_types[i]],
            # label="Kernel density estimate",
            label=""
        )

        # Final formatting
        Plots.ylims!(0, round(maximum([n₁; u.density]), digits=2)+0.01)
        npoints = count(sample_cats[fig_types[i]])
        Plots.annotate!(((0.03, 0.97), (fig_names[i] * "\nn = $npoints", 12, :left, :top)))
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


## --- Compare trace elements to Rudnick and Gao, 2014 (10.1016/B978-0-08-095975-7.00301-6)
    # Load all estimates from Rudnick and Gao. 
    rg_all = importdataset("data/rudnick_gao_2014_table1-2.csv",  ',', importas=:Tuple)
    ucc = importdataset(ucc_out, '\t', importas=:Tuple)

    # rg = Dict(zip(rg.Element, rg.Percent))
    ucc = Dict(zip(ucc.element, ucc.bulk))                              # My estimate
    rg = Dict(zip(rg_all.Element, rg_all.This_Study))                   # Rudnick and Gao
    tm = Dict(zip(rg_all.Element, rg_all.Taylor_and_McLennan_1985))     # Taylor and McLennan
    sw = Dict(zip(rg_all.Element, rg_all.Shaw_et_al__1967))             # Shaw et al.
    kc = Dict(zip(rg_all.Element, rg_all.Condie_1993))                  # Condie
    sg = Dict(zip(rg_all.Element, rg_all.Gao_et_al_1998))               # Gao et al.

    estimates = [rg, tm, sw, kc, sg]

    # Convert units to percent for normalization
    units = Dict(zip(rg_all.Element, rg_all.Units))
    for e in eachindex(estimates)
        for k in keys(estimates[e])
            if units[k] == "percent"
                continue
            elseif units[k] == "ppm"
                estimates[e][k] = estimates[e][k] / 10_000
            elseif units[k] == "ppb"
                estimates[e][k] = estimates[e][k] / 10_000_000
            end
        end
    end

    # Remove volatiles from my estimate and normalize all estimates to 100%
    delete!(ucc_recast, "Volatiles")
    ucc = Dict(zip(keys(ucc), normalize!(collect(values(ucc)))))
    for e in eachindex(estimates)
        estimates[e] = Dict(
            zip(keys(estimates[e]), normalize!(collect(values(estimates[e]))))
        )
    end
    
    # Get the elements we want to analyze
    cnorm = get_chondrite_norm()
    REEs = keys(cnorm)

    # Convert estimates into normalized REE space, recasting to mg/g
    ucc_REE = NamedTuple{Tuple(REEs)}([ucc[string(i)] / cnorm[i] * 10000 for i in REEs])

    # Spider diagram, add an empty value for Pm
    all_REEs = get_REEs()
    i = findfirst(x->x==:Pm, all_REEs)
    x = collect([1:i-1; i+1:length(all_REEs)])

    all_REEs = string.(all_REEs)
    all_REEs[i] = ""

    labels = ["Rudnick and Gao, 2014", "Taylor and McLennan, 1985/1995", 
        "Shaw et al., 1967", "Condie, 1993", "Gao et al., 1998"
    ]
    estcolors = [:blue, :green, :hotpink, :purple, :red]
    shapes = [:utriangle, :dtriangle, :star5, :diamond, :x]

    h = Plots.plot(
        ylabel="Chondrite Normalized",
        fg_color_legend=:white,
        framestyle=:box,
        grid=false,
        yaxis=:log10,
        ylims=(10^0, 10^3),
        yticks=(10.0.^(0:3), ("1", "10", "100", "1000")),
        xticks=(x, string.(REEs)),
        yminorticks=log.(1:10),
    )
    for e in eachindex(estimates)
        # Normalize REEs, recasting into mg/g
        REE_i = NamedTuple{Tuple(REEs)}([
            estimates[e][string(i)] / cnorm[i] * 10000 for i in REEs
        ])

        # Plot
        Plots.plot!(h, x, collect(values(REE_i)),
            markershape=shapes[e], color=estcolors[e], msc=estcolors[e],
            label=labels[e],
        )
    end
    Plots.plot!(h, x, collect(values(ucc_REE)),
        markershape=:circle, color=:darkorange, msc=:darkorange,
        label="This study",
    )

    display(h)
    savefig("$filepath/spidergram.pdf")


## --- Slope vs. erosion rate
    using Isoplot
    using StatsBase
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
    ey_bound = [(l[i]*2, m[i]*2) for i in eachindex(l)]

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
        xlims=(0, 675),
        tickfontsize=12,
        labelfontsize=12,
        legendfontsize=12,
        size=(600,400),
        # left_margin=(10,:px), bottom_margin=(10,:px)
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
        markersize = 5,
    )
    Plots.plot!(h, 1:length(model), model, label="",
        # label="exp(slope * 0.0098) + 2.97)", 
        color=:black, linewidth=3
    )
    
    display(h)
    savefig("$filepath/erosionslope.pdf")


## --- End of file