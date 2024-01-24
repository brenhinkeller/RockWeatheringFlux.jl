## --- Set up 
    # 2D histograms of igneous and metamorphic samples as a function of age and silica 
    # content. After Figure 6.9 from Keller, 2016 (10.31237/osf.io/q7yra)

    # Load data and base packages
    include("Definitions.jl")

    # Preallocate / Local definitions
    nsims = Int(1e7)                         # 10 M simulations
    SiO2_error = 5.0                         # Large SiO₂ error to smooth data
    age_min_error = 0.15                     # Minimum 15% age error

    xmin, xmax, xbins = 40, 80, 240          # Silica
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380         # Age
    yedges = ymin:(ymax-ymin)/ybins:ymax


## --- EarthChem 
    notsed = bulk_cats.ign .| bulk_cats.met
    t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude) & !isnan(bulk.Age) & notsed

    s = vec(isnan.(bulk.Age))
    ageuncert = bulk.Age .* age_min_error
    for i in eachindex(ageuncert)
        calc_uncert = (bulk.Age_Max[i] .- bulk.Age_Min[i])/2
        ageuncert[i] = ifelse(calc_uncert > ageuncert[i], calc_uncert, ageuncert[i])
    end

    k = invweight(bulk.Latitude[t], bulk.Longitude[t], bulk.Age[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2[t] bulk.Age[t]]
    uncertainty = [fill(SiO2_error, count(t)) ageuncert[t]]
    simbulk = bsresample(data, uncertainty, nsims, p)

    # Make a 2d-histogram / heatmap. Normalize each time step to between 0 and 1
    out_bulk = zeros(ybins, xbins)                # Preallocate
    for i = 1:ybins
        # Filter for samples in this age bin
        t = @. yedges[i] <= simbulk[:,2] < yedges[i+1]

        # Count and normalize distribution of silica
        c, n = bincounts(simbulk[:,1][t], xmin, xmax, xbins)
        n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
        # n = float(n) ./ nansum(float(n) .* step(c))

        # Put the output in the output array
        out_bulk[i,:] .= n
    end


## --- Matched samples
    notsed = macro_cats.ign .| macro_cats.met
    t = @. !isnan(macrostrat.rocklat) & !isnan(macrostrat.rocklon) & !isnan(macrostrat.age) & notsed

    s = vec(isnan.(macrostrat.age))
    ageuncert = macrostrat.age .* age_min_error
    for i in eachindex(ageuncert)
        calc_uncert = (macrostrat.agemax[i] .- macrostrat.agemin[i])/2
        ageuncert[i] = ifelse(calc_uncert > ageuncert[i], calc_uncert, ageuncert[i])
    end

    k = invweight_age(macrostrat.age[t])     # Already spatially corrected
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [mbulk.SiO2[t] macrostrat.age[t]]
    uncertainty = [fill(SiO2_error, count(t)) ageuncert[t]]
    sim_mbulk = bsresample(data, uncertainty, nsims, p)

    # Make a 2d-histogram / heatmap. Normalize each time step to between 0 and 1
    out_mbulk = zeros(ybins, xbins)                # Preallocate
    for i = 1:ybins
        # Filter for samples in this age bin
        t = @. yedges[i] <= sim_mbulk[:,2] < yedges[i+1]

        # Count and normalize distribution of silica
        c, n = bincounts(sim_mbulk[:,1][t], xmin, xmax, xbins)
        n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
        # n = float(n) ./ nansum(float(n) .* step(c))

        # Put the output in the output array
        out_mbulk[i,:] .= n
    end


## --- Build plot 
    # EarthChem
    h1 = Plots.plot(
        ylims=(0,380),
        yticks=(0:50:380, string.(0:500:3800)),
        yflip=true,
        xlabel="SiO₂ [wt.%]", 
        ylabel="Age [Ma]",
        size=(600,500),
    )
    Plots.heatmap!(h1, out_bulk, 
        # aspectratio=:equal,
        # colorbar_title="Relative Sample Density",
        colorbar=false,
        color=c_gradient,
        title="A. Resampled Geochemical Samples\n"
    )

    # Matched samples
    h2 = Plots.plot(
        ylims=(0,380),
        # yticks=(0:50:380, string.(0:500:3800)),
        yticks=false,
        yflip=true,
        xlabel="SiO₂ [wt.%]", 
        # ylabel="Age [Ma]",
        size=(600,500),
    )
    Plots.heatmap!(h2, out_mbulk, 
        # aspectratio=:equal,
        # colorbar_title="Relative Sample Density",
        colorbar=false,
        color=c_gradient,
        title="B. Resampled Matched Samples\n"
    )

    # Assemble all plots with a common color bar
    l = @layout [a{0.95w} b{0.05w}]
    a = Plots.plot(h1, h2, layout=(1,2),
        size=(1000, 400),
        framestyle=:box,
        grid=false,
        fontfamily=:Helvetica,
        xticks=(0:60:240, string.(40:10:80)),
        left_margin=(25,:px), bottom_margin=(15,:px),
        titleloc=:left, titlefont = font(12),
    )
    b = Plots.heatmap(rand(2,2), clims=(0,1), 
        framestyle=:none, color=c_gradient, colorbar_title="Relative Sample Density", 
        lims=(-1,0)
    )
    h = Plots.plot(a, b, layout=l)

    display(h)
    savefig("$filepath/silica_heatmap.pdf")


## --- Alternatively, do a series of curves every 100 million years
    # I'm particularly interested in the matched samples, and if there's points where
    # the Daly Gap isn't present
    using KernelDensity

    # Set up a cute little color pallette 
    p = Plots.palette(:glasgow, 50)[1:38]

    # Figure out how many rows to grab to get the bin size we want 
    target = collect(ymin:100:ymax)     # Age edges
    target[end] = 3800                  # Make sure we include all samples to 3800
    pᵢ = 1:round(Int, length(p)/length(target)):38

    # Preallocate for subplots 
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 2)
    fig_types = (:volc, :plut, :ign)
    fig_names = ("A. Matched Data (Temporal Resample)", "B. Whole Dataset (Spatiotemporal Resample)")
    sim = [sim_mbulk, simbulk]

    for f in eachindex(fig)
        h = Plots.plot(
            xlims=(40,80),
            yticks=false,
            grid=false,
            xlabel="SiO₂ [wt.%]", 
            title=fig_names[f],
            size=(600,400),
            framestyle=:box
        )

        # Use the raw simulation data for KDE reasons
        for i in 1:(length(target)-1)
            # Get the samples in these age bins
            t = @. (target[i] <= sim[f][:,2] < target[i+1]) & !isnan(sim[f][:,1])

            # Estimate KDE and plot
            u = kde(sim[f][:,1][t])
            plot!(h, u.x, u.density, label="", color=p[pᵢ[i]])
        end
        
        fig[f] = h
    end

    # Get subplots together
    ylabel!(fig[1], "Relative Abundance")
    a = Plots.plot(fig..., layout=(1, 2), size=(1200,500))

    # We love legends. Very hack-y colorbar
    b = Plots.plot(1:38, ones(38), label="",
        seriestype=:bar, barwidths=1.1,
        color=p, 
        linealpha=0,
        xlims=(0,38), xticks=([0.5, 19, 37], string.(0:1900:3800)),
        xlabel="Age [Ma]", 
        yticks=false, yaxis=false, 
        tick_direction=:none, bordercolor=:white,
        size=(600, 80),
    )

    # Assemble everything
    l = @layout [a{0.95h} 
                b{0.05h}]
    h = Plots.plot(a, b, layout=l, size=(1200,600), 
        left_margin=(40,:px), bottom_margin=(15,:px)
    )
    display(h)


## --- End of file 