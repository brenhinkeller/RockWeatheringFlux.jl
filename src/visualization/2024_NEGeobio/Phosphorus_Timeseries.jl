## --- Set up 

    # Load data and base packages
    include("Definitions.jl")

    # Definitions 
    age_error = 0.05        # Minimum age error (%)
    age_error_abs = 50      # Minimum age error (Ma)
    P2O5_err = 0.2          # Error given in Keller et al., 2015
    xmin, xmax, nbins = 0, 3800, 38
    c = cntr(xmin:(xmax-xmin)/nbins:xmax)


## --- Resample phosphorus data using matched geochemical dataset
    # # Define intervals of interest
    # include_minor!(match_cats);
    # target = (
    #     sed = match_cats.sed,
    #     clastic = match_cats.siliciclast .| match_cats.shale,
    #     carb = match_cats.sed,
    #     ign = match_cats.ign,
    #     mafic = match_cats.ign .& (43 .< mbulk.SiO2 .<= 51),
    #     int = match_cats.ign .& (51 .< mbulk.SiO2 .<= 62),
    #     felsic = match_cats.ign .& (62 .< mbulk.SiO2 .<= 73),
    # )
    # simout = NamedTuple{keys(target)}(Array{Float64}(undef, nbins, 3) for _ in keys(target))

    # # Calculate sample age, use the sample age unless missing, then use map age
    # t = @. !isnan(mbulk.Age);
    # sampleage = copy(mbulk.Age);
    # ageuncert = nanadd.(mbulk.Age_Max, .- mbulk.Age_Min) ./ 2;
    # sampleage[t] .= macrostrat.age[t]
    # ageuncert[t] .= nanadd.(macrostrat.agemax[t], .- macrostrat.agemin[t]) ./ 2;
    # for i in eachindex(ageuncert)
    #     ageuncert[i] = max(sampleage[i]*age_error, ageuncert[i], age_error_abs)
    # end
    # t = @. !isnan.(sampleage);

    # # Abundance of phosphorus [wt.%] in 100 Ma age bins
    # for key in keys(target)
    #     s = t .& target[key]
    #     k = invweight_age(sampleage[s])
    #     p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    #     c,m,el,eu = bin_bsr(sampleage[s], mbulk.P2O5[s], xmin, xmax, nbins;
    #         x_sigma = ageuncert[s],
    #         y_sigma = fill(P2O5_err, count(s)),
    #         nresamplings = 10_000,
    #         sem = :pctile,
    #         p = p,
    #     )
    #     simout[key] .= [m eu el]
    # end

    # # Save data to a file 
    # fid = h5open("src/visualization/2024_NEGeobio/simout_phosphorus.h5", "w")
    # g = create_group(fid, "vars")
    # for key in keys(simout) 
    #     write(g, "$key", simout[key])
    # end
    # close(fid)

    
## --- Load data from file 
    fid = h5open("src/visualization/2024_NEGeobio/simout_phosphorus.h5", "r")
    target = Tuple(Symbol.(keys(fid["vars"])))
    simout = NamedTuple{target}(Array{Float64}(undef, nbins, 3) for _ in target)
    for key in target
        simout[key] .= read(fid["vars"]["$key"])
    end
    close(fid)


## --- Create base plot
    h = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlims=(-50, 3850),
        # ylims=(0.05,0.55),
        grid=false,
        fg_color_legend=:white,
        # xlabel="Age [Ma.]",
        # ylabel="P concentration [wt.%]",
        legendfont=11,
        tickfontsize=12, labelfontsize=12,
        legend=:topright
    );


## --- Sedimentary 
    h1 = deepcopy(h)

    # Plots.plot!(h, c, simout.carb[:,1], label="Carbonate",
    #     yerror=(simout.carb[:,2], simout.carb[:,3],), 
    #     color=pal[9], lcolor=pal[8], msc=:auto,
    #     seriestype=:scatter
    # )
    Plots.plot!(h1, c, simout.clastic[:,1], label="Siliciclastic",
        yerror=(simout.clastic[:,2], simout.clastic[:,3],), 
        color=pal[2], lcolor=pal[2], msc=:auto, alpha=0.5,
        seriestype=:scatter, markersize=4, linewidth=2
    )
    Plots.plot!(h1, c, simout.sed[:,1], label="All Seds",
        yerror=(simout.sed[:,2], simout.sed[:,3],), 
        color=pal[6], lcolor=pal[6], msc=:auto,
        seriestype=:scatter, markersize=4, linewidth=2
    )

    Plots.vline!(h1, [753.3], color=:black, linestyle=:dash, label="")
    Plots.annotate!([(803.3, 0.53, text("753.3 Ma.", 
        12, :right, :top, :black, rotation=90)
    )])

    Plots.vline!(h1, [2500], color=:black, linestyle=:dash, label="")
    Plots.annotate!([(2550, 0.53, text("2500 Ma.", 
        12, :right, :top, :black, rotation=90)
    )])

    ylims!(0.06, 0.56)
    yticks!(h2, 0.1:0.1:0.4, string.(0.1:0.1:0.4))

    display(h1)
    savefig(h1, "$filepath/phosphorus_sed.pdf")


## --- Igneous 
    h2 = deepcopy(h)

    Plots.plot!(h2, c, simout.mafic[:,1], label="Mafic",
        yerror=(simout.mafic[:,2], simout.mafic[:,3],), 
        color=pal[2], lcolor=pal[2], msc=pal[1], alpha=0.5,
        seriestype=:scatter, markersize=4, linewidth=2,
    )
    # Plots.plot!(h2, c, simout.int[:,1], label="Intermediate",
    #     yerror=(simout.int[:,2], simout.int[:,3],), 
    #     color=pal[3], lcolor=pal[3], msc=pal[3], alpha=0.5,
    #     seriestype=:scatter, markersize=4, linewidth=2
    # )
    Plots.plot!(h2, c, simout.felsic[:,1], label="Felsic",
        yerror=(simout.felsic[:,2], simout.felsic[:,3],), 
        color=pal[4], lcolor=pal[4], msc=pal[5], alpha=0.5,
        seriestype=:scatter, markersize=4, linewidth=2
    )
    Plots.plot!(h2, c, simout.ign[:,1], label="All Igneous",
        yerror=(simout.ign[:,2], simout.ign[:,3],), 
        color=pal[6], lcolor=pal[6], msc=:auto,
        seriestype=:scatter, markersize=4, linewidth=2
    )

    Plots.vline!(h2, [753.3], color=:black, linestyle=:dash, label="")
    Plots.annotate!([(803.3, 0.4, text("753.3 Ma.", 
        12, :right, :top, :black, rotation=90)
    )])

    Plots.vline!(h2, [2500], color=:black, linestyle=:dash, label="")
    Plots.annotate!([(2550, 0.4, text("2500 Ma.", 
        12, :right, :top, :black, rotation=90)
    )])

    ylims!(0.07, 0.43)
    yticks!(h2, 0.1:0.1:0.4, string.(0.1:0.1:0.4))

    display(h2)
    savefig(h2, "$filepath/phosphorus_ign.pdf")


## --- End of file