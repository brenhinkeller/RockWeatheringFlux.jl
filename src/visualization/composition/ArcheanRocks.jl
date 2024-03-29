## --- Set up 
    # Predicted (bulk geochemical data) and observed (matched data) silica distributions 
    # in Archean rocks

    # Load data and base packages
    include("Definitions.jl")

    # Preallocate / Local definitions
    nsims = Int(1e7)                         # 10 M simulations
    SiO2_error = 1.0                         # Assumed SiO₂ error
    age_error = 0.05                         # Minimum age error (%)
    age_error_abs = 50                       # Minimum age error (Ma)

    age_ceil = 2500                          # Youngest age sample to consider


## --- Resample bulk geochemical dataset 
    # Preallocate 
    simbulk = Array{Float64}(undef, nsims, 2)

    # Compute age uncertainties 
    ageuncert = nanadd.(bulk.Age_Max, .- bulk.Age_Min) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = max(bulk.Age[i]*age_error, ageuncert[i], age_error_abs)
    end

    # Restrict to igneous Archean samples with data
    t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude) & (bulk.Age .> age_ceil);
    t .&= bulk_cats.ign;

    # Resample
    k = invweight(bulk.Latitude[t], bulk.Longitude[t], bulk.Age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2[t] bulk.Age[t]]
    uncertainty = [fill(SiO2_error, count(t)) ageuncert[t]]
    simbulk .= bsresample(data, uncertainty, nsims, p)


## --- Resample matched geochemical dataset 
    # Preallocate 
    sim_mbulk = Array{Float64}(undef, nsims, 2)

    # Use the sample age, unless the sample doesn't have an age: then use map age
    t = @. !isnan(mbulk.Age);
    sampleage = copy(mbulk.Age);
    ageuncert = nanadd.(mbulk.Age_Max, .- mbulk.Age_Min) ./ 2;
    sampleage[t] .= macrostrat.age[t]
    ageuncert[t] .= nanadd.(macrostrat.agemax[t], .- macrostrat.agemin[t]) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = max(sampleage[i]*age_error, ageuncert[i], age_error_abs)
    end

    # Restrict to igneous Archean samples with data
    t = @. (sampleage > age_ceil) & match_cats.ign;

    # Resample
    k = invweight_age(sampleage[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [mbulk.SiO2[t] sampleage[t]]
    uncertainty = [fill(SiO2_error, count(t)) ageuncert[t]]
    sim_mbulk .= bsresample(data, uncertainty, nsims, p)


## --- Build plot
    # Bulk geochemical data 
    h1 = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        grid=false,
        xlabel="SiO2 [wt.%]", 
        ylabel="Relative Abundance",
        yticks=false,
        xlims=(40,80),
        legend=:topleft, 
        fg_color_legend=:white, 
        left_margin=(15,:px)
    )
    c, n = bincounts(simbulk[:,1], 40, 80, 80)
    n₁ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(h1, c, n₁, 
        seriestype=:bar, barwidths=1,
        color=colors.komatiite, linecolor=:match,
        title="A. Bulk Geochemical Data\n", label="",
        ylims=(0, round(maximum(n₁), digits=2)+0.01,)
    )

    # Matched data
    h2 = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        grid=false,
        xlabel="SiO2 [wt.%]", 
        # ylabel="Relative Abundance",
        yticks=false,
        xlims=(40,80),
        legend=:topleft, 
        fg_color_legend=:white, 
        left_margin=(15,:px)
    )
    c, n = bincounts(sim_mbulk[:,1], 40, 80, 80)
    n₂ = float(n) ./ nansum(float(n) .* step(c))
    Plots.plot!(h2, c, n₂, 
        seriestype=:bar, barwidths=1,
        color=colors.plut, linecolor=:match,
        title="B. Matched Samples\n", label="",
        ylims=(0, round(maximum(n₂), digits=2)+0.01,)
    )

    # Both plots together
    h = Plots.plot(h1, h2, layout=(1,2), size=(1200,450), 
        bottom_margin=(40,:px), left_margin=(50,:px), right_margin=(10,:px),
        labelfontsize=labelfontsize, titlefont=titlefontsize, tickfontsize=tickfontsize
    )
    display(h)
    savefig("$filepath/archeansilica.pdf")

    
## --- End of file