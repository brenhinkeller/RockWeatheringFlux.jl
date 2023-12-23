## --- Set up 
    # Experiments with silica distributions and rock classes of Archean rocks

    # Unique packages
    using StatsPlots

    # Load data and base packages
    if !@isdefined(filepath)
        include("Definitions.jl")
    end


## --- Rock class distribution in the Archean
    # Filter for just Archean samples
    archean = macrostrat.age .>= 2500;

    # Define felsic / intermediate / mafic / ultramafic types
    fels = (:rhyolite, :trondhjemite, :tonalite, :granodiorite, :granite,)
    intr = (:andesite, :dacite, :diorite,) 
    mafc = (:basalt, :gabbro,)
    umaf = (:komatiite, :peridotite, :pyroxenite,)
    alkn = (:alk_volc, :alk_plut)

    types = [alkn, fels, intr, mafc, umaf]
    counts = zeros(Int, length(types))
    for t in eachindex(types)
        for k in types[t]
            counts[t] += count(macro_cats[k] .& archean)
        end
    end

    # Visualize!
    xlabels = 
    a = Plots.plot(
        framestyle=:none,
        fontfamily=:Helvetica,
        grid=false,
        xticks=false,
        yticks=false,
        ylims=(0.5, 1), 
        legend=false,
        size=(800,200)
    )
    StatsPlots.groupedbar!(a, counts', 
        bar_position = :stack, orientation=:h, bar_width=1,
        group = ["Alkaline", "Felsic", "Intermediate", "Mafic", "Ultramafic"],
        color = [colors.alk_plut, colors.rhyolite, colors.andesite, colors.basalt, 
            colors.peridotite], 
        linecolor=:match
    )

    b = plot(
        framestyle=:none,
        xticks=false, yticks=false,
        ylims=(1,2), xlims=(1,2),
        fg_color_legend=:white, 
        legendfontsize = 12, 
        legend=:outerbottom, 
        legendcolumns=5,
        legendtitle="",
        size=(1200,200),
    )
    plot!(b, [0], [0], label=" Ultramafic", color=colors.peridotite, seriestype=:bar)
    plot!(b, [0], [0], label=" Mafic", color=colors.basalt, seriestype=:bar)
    plot!(b, [0], [0], label=" Intermediate", color=colors.andesite, seriestype=:bar)
    plot!(b, [0], [0], label=" Felsic", color=colors.rhyolite, seriestype=:bar)
    plot!(b, [0], [0], label=" Alkaline", color=colors.alk_plut, seriestype=:bar)

    l = @layout [
        a{0.95h} 
        b{0.05h}
    ]
    h = plot(a, b, layout=l)

    display(h)
    savefig("$filepath/archeanrockclass.pdf")


## --- Silica distribution across the Archean
    # Preallocate
    nsims = Int(1e6)                         # 1 M simulations
    SiO2_error = 1.0                         # Large SiO₂ error to smooth data
    age_error = 0.05                         # Assumed % age error, if unknown

    # Resample (spatiotemporal weights) unmatched EarthChem samples
    notsed = bulk_cats.ign .| bulk_cats.met
    t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude) & !isnan(bulk.Age) & notsed

    ageuncert = Array{Float64}(undef, count(t), 1)
    ageuncert .= (bulk.Age_Max[t] .- bulk.Age_Min[t])/2
    s = vec(isnan.(ageuncert))
    ageuncert[s] .= bulk.Age[t][s] .* age_error

    k = invweight(bulk.Latitude[t], bulk.Longitude[t], bulk.Age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2[t] bulk.Age[t]]
    uncertainty = [fill(SiO2_error, count(t)) ageuncert]
    simbulk = bsresample(data, uncertainty, nsims, p)

    # Resample (temporal weights) matched EarthChem samples 
    notsed = macro_cats.ign .| macro_cats.met
    t = @. !isnan(macrostrat.rocklat) & !isnan(macrostrat.rocklon) & !isnan(macrostrat.age) & notsed

    ageuncert = Array{Float64}(undef, count(t), 1)
    ageuncert .= (macrostrat.agemax[t] .- macrostrat.agemin[t])/2
    s = vec(isnan.(ageuncert))
    ageuncert[s] .= macrostrat.age[t][s] .* age_error

    k = invweight_age(macrostrat.age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [mbulk.SiO2[t] macrostrat.age[t]]
    uncertainty = [fill(SiO2_error, count(t)) ageuncert]
    simmacro = bsresample(data, uncertainty, nsims, p)

## --
    # Filter archean samples
    old_bulk = simbulk[:,2] .<= 2500;
    old_macro = simmacro[:,2] .<= 2500;

    # Build plot
    h = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        grid=false,
        xlabel="SiO2 [wt.%]", 
        ylabel="Relative Abundance",
        xlims=(40,80),
        legend=:topright, 
        fg_color_legend=:white, 
    )

    # Resampled unmatched EarthChem
    c, n = bincounts(simbulk[:,1][old_bulk], 40, 80, 80)
    n₁ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(h, c, n₁, 
        seriestype=:bar, linewidth=2,
        color=:grey, linecolor=:match, alpha=0.5,
        label="EarthChem"
    )

    # Resampled matched EarthChem
    c, n = bincounts(simmacro[:,1][old_macro], 40, 80, 80)
    n₂ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(h, c, n₂, 
        seriestype=:bar, linewidth=2,
        color=:steelblue, linecolor=:match, alpha=0.5,
        label="Matched Samples"
    )

    ylims!(0, round(maximum([n₁; n₂]), digits=2)+0.01,)
    display(h)

    
## --- End of file