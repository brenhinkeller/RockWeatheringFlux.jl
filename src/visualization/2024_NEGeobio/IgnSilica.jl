## --- Set up 
    # Igneous silcia distributions, including time-variance. Allows a little more 
    # granular control (and to change the colors) than the code for the paper

    # Load data and base packages
    include("Definitions.jl")
    using MAT, KernelDensity

    # Definitions
    here = @. !isnan(mbulk.SiO2);            # Filter NaNs from my data
    nsims = Int(1e7)
    SiO₂_error = 1.0                         # Error in volcanic.mat is 0.01
    age_error = 0.05                         # Minimum age error (%)
    age_error_abs = 50                       # Minimum age error (Ma)

    xmin, xmax, xbins = 40, 80, 240          # Silica
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380         # Age
    yedges = ymin:(ymax-ymin)/ybins:ymax

    # Rock classes of interest
    target = (:volc, :plut, :ign)


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


## --- Create time-integrated histograms 
    # Preallocate
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    fig_colors = (volc=pal[2], plut=pal[4], ign=pal[6])
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
            tickfontsize=16, labelfontsize=18,
            left_margin=(15,:px), right_margin=(15,:px),
        )

        # Raw data
        c, n = bincounts(mbulk.SiO2[match_cats[target[i]]], 40, 80, 80)
        n₁ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₁, 
            seriestype=:bar, barwidths=0.5,
            color=fig_colors[target[i]], linecolor=:match, alpha=0.25,
            label=""
        )

        # Kernel density estimate 
        u = kde(mbulk.SiO2[match_cats[target[i]] .& here])
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=fig_colors[target[i]], linecolor=fig_colors[target[i]],
            label=""
        )

        # Keller et al., 2015
        c, n = bincounts(simVP[i], 40, 80, 80)
        n₂ = float(n) ./ nansum(float(n) .* step(c))
        Plots.plot!(c, n₂, 
            seriestype=:path, linewidth=2,
            color=:black, linecolor=:black,
            linestyle=:dot,
            label=""
        )

        # Final formatting
        Plots.ylims!(0, round(maximum([n₁; n₂; u.density]), digits=2)+0.02)
        display(h)
        savefig(h, "$filepath/histogram_$(target[i]).pdf")
    end


    # Shared legend
    h1 = Plots.plot(framestyle=:none, xlims=(1,2), ylims=(1,2), fg_color_legend=:white)
    Plots.plot!(h1, [0],[0], color=:white, linecolor=:match, label=" ")
    Plots.plot!(h1, [0],[0], color=fig_colors.volc, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" ",)
    Plots.plot!(h1, [0],[0], color=fig_colors.plut, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" ")
    Plots.plot!(h1, [0],[0], color=fig_colors.ign, linecolor=:match, seriestype=:bar, 
        alpha=0.25, label=" ")

    Plots.plot!(h1, [0],[0], color=:white, linecolor=:match, label="")
    Plots.plot!(h1, [0], [0], linewidth=2, color=:black, linestyle=:dot,
        label=" "
    )
    display(h1)
    savefig(h1, "$filepath/histogram_legend.pdf")


## --- Resample matched geochemical samples 
    # # Preallocate 
    # sim_mbulk = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

    # # Use the sample age, unless the sample doesn't have an age: then use map age
    # t = @. !isnan(mbulk.Age);
    # sampleage = copy(mbulk.Age);
    # ageuncert = nanadd.(mbulk.Age_Max, .- mbulk.Age_Min) ./ 2;
    # sampleage[t] .= macrostrat.age[t]
    # ageuncert[t] .= nanadd.(macrostrat.agemax[t], .- macrostrat.agemin[t]) ./ 2;
    # for i in eachindex(ageuncert)
    #     ageuncert[i] = max(sampleage[i]*age_error, ageuncert[i], age_error_abs)
    # end

    # # Restrict to only samples with data and resample 
    # t = @. !isnan.(sampleage);
    # for key in target 
    #     s = t .& match_cats[key]
    #     k = invweight_age(sampleage[s])
    #     p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    #     data = [mbulk.SiO2[s] sampleage[s]]
    #     uncertainty = [fill(SiO₂_error, count(s)) ageuncert[s]]
    #     sim_mbulk[key] .= bsresample(data, uncertainty, nsims, p)
    # end

    # # Sort matched data into time bins, normalize, and save to a file
    # fid = h5open("src/visualization/2024_NEGeobio/simout_ignsilica.h5", "w")
    # g = create_group(fid, "vars")
    # for i in eachindex(target)
    #     # Sort into time bins and normalize
    #     out_mbulk = zeros(ybins, xbins)
    #     for j = 1:ybins
    #         t = @. yedges[j] <= sim_mbulk[target[i]][:,2] < yedges[j+1]
    #         c, n = bincounts(sim_mbulk[target[i]][:,1][t], xmin, xmax, xbins)
    #         out_mbulk[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    #     end

    #     # Save data
    #     write(g, "$(target[i])", out_mbulk)
    # end
    # close(fid)


## --- Create 2D age / silica histogram 
    # Open normalized data file 
    fid = h5open("src/visualization/2024_NEGeobio/simout_ignsilica.h5", "r")
    simout = NamedTuple{target}(read(fid["vars"]["$key"]) for key in target)
    close(fid)

    # Base plot 
    h = Plots.plot(
        fontfamily=:Helvetica,
        framestyle=:box,
        ylims=(0,ybins),
        yticks=(0:ybins/7.6*2:ybins, string.(0:1000:3800)),
        yflip=true,
        xticks=(0:xbins/4:xbins, string.(40:10:80)),
        # xlabel="SiO2 [wt.%]", ylabel="Age [Ma.]",
        size=(650,500),
        tickfontsize=16,labelfontsize=16,
    )
    
    for i in eachindex(target)
        # Plot
        h1 = deepcopy(h)
        Plots.heatmap!(h1, simout[target[i]], 
            colorbar=false,
            color=reverse(pal),
        )

        display(h1)
        savefig(h1, "$filepath/heatmap_$(target[i]).pdf")
    end

    # Make a color bar legend
    h2 = Plots.heatmap(rand(2,2), clims=(0,1), 
        framestyle=:none, color=:managua, 
        colorbar_title="Relative Sample Density", # ticks=false,
        lims=(-1,0)
    )
    savefig(h2, "$filepath/heatmap_colorbar.pdf")


## --- End of file