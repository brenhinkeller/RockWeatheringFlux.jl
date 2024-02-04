## --- Set up 
    # 2D histograms of igneous and metamorphic samples as a function of age and silica 
    # content. After Figure 6.9 from Keller, 2016 (10.31237/osf.io/q7yra)

    # Load data and base packages
    include("Definitions.jl");

    # Preallocate / Local definitions
    nsims = Int(1e7)                         # 10 M simulations
    SiO2_error = 1.0                         # Assumed SiO₂ error
    age_error = 0.05                         # Minimum age error

    xmin, xmax, xbins = 40, 80, 240          # Silica
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380         # Age
    yedges = ymin:(ymax-ymin)/ybins:ymax

    # Rock classes to resample
    target = (:ign, :plut, :volc)


## --- Resample (spatiotemporal) bulk geochemical data 
    # Preallocate 
    simbulk = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

    # Compute age uncertainties 
    ageuncert = nanadd.(bulk.Age_Max, .- bulk.Age_Min) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = max(bulk.Age[i]*age_error, ageuncert[i])
    end

    # Restrict to samples with data and resample 
    t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude) & !isnan(bulk.Age) & bulk_cats.ign
    for key in target 
        s = t .& bulk_cats[key]
        k = invweight(bulk.Latitude[t], bulk.Longitude[t], bulk.Age[t])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        data = [bulk.SiO2[t] bulk.Age[t]]
        uncertainty = [fill(SiO2_error, count(t)) ageuncert[t]]
        simbulk[key] .= bsresample(data, uncertainty, nsims, p)
    end


## --- Resample (temporal) matched samples
    # Preallocate 
    sim_mbulk = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

    # Use the sample age, unless the sample doesn't have an age: then use map age
    t = @. !isnan(mbulk.Age);
    sampleage = copy(mbulk.Age);
    ageuncert = nanadd.(mbulk.Age_Max, .- mbulk.Age_Min) ./ 2;
    sampleage[t] .= macrostrat.age[t]
    ageuncert[t] .= nanadd.(macrostrat.agemax[t], .- macrostrat.agemin[t]) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = max(sampleage[i]*age_error, ageuncert[i])
    end

    # Restrict to only samples with data and resample 
    t = @. !isnan.(sampleage);
    for key in target 
        s = t .& match_cats[key]
        k = invweight_age(sampleage[s])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        data = [mbulk.SiO2[s] sampleage[s]]
        uncertainty = [fill(SiO2_error, count(s)) ageuncert[s]]
        sim_mbulk[key] .= bsresample(data, uncertainty, nsims, p)
    end


## --- Build plot (bulk igneous data)
    # Geochemical data
    out_bulk = zeros(ybins, xbins)
    for i = 1:ybins
        t = @. yedges[i] <= simbulk.ign[:,2] < yedges[i+1]
        c, n = bincounts(simbulk.ign[:,1][t], xmin, xmax, xbins) 
        out_bulk[i,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    end

    # Matched data
    out_mbulk = zeros(ybins, xbins)
    for i = 1:ybins
        t = @. yedges[i] <= sim_mbulk.ign[:,2] < yedges[i+1]
        c, n = bincounts(sim_mbulk.ign[:,1][t], xmin, xmax, xbins)
        out_mbulk[i,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    end

    # Bulk geochemical data
    h1 = Plots.plot(
        ylims=(0,ybins),
        yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
        yflip=true,
        xlabel="SiO2 [wt.%]", 
        ylabel="Age [Ma]",
        size=(600,500),
    )
    Plots.heatmap!(h1, out_bulk, 
        colorbar=false,
        color=c_gradient,
        title="A. Bulk Geochemical Data\n"
    )

    # Matched samples
    h2 = Plots.plot(
        ylims=(0,ybins),
        yticks=false,
        yflip=true,
        xlabel="SiO2 [wt.%]", 
        size=(600,500),
    )
    Plots.heatmap!(h2, out_mbulk, 
        colorbar=false,
        color=c_gradient,
        title="B. Matched Samples\n"
    )

    # Assemble all plots with a common color bar
    l = @layout [a{0.95w} b{0.05w}]
    a = Plots.plot(h1, h2, layout=(1,2),
        size=(1000, 400),
        framestyle=:box,
        grid=false,
        fontfamily=:Helvetica,
        xticks=(0:xbins/4:xbins, string.(40:10:80)),
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


## --- Build plot (subdivide igneous into volcanic and plutonic classes)
    # Preallocate 
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    subclass = (:ign, :volc, :plut)
    labels = ("Igneous", "Volcanic", "Plutonic")
    subfig = (("A", "B"), ("C", "D"), ("E", "F"))

    p = Progress(6*ybins, desc="Makin' plots:");
    for i in eachindex(fig)
        # Geochemical data
        out_bulk = zeros(ybins, xbins)
        for j = 1:ybins
            t = @. yedges[j] <= simbulk[subclass[i]][:,2] < yedges[j+1]
            c, n = bincounts(simbulk[subclass[i]][:,1][t], xmin, xmax, xbins) 
            out_bulk[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
            next!(p)
        end

        # Matched data
        out_mbulk = zeros(ybins, xbins)
        for j = 1:ybins
            t = @. yedges[j] <= sim_mbulk[subclass[i]][:,2] < yedges[j+1]
            c, n = bincounts(sim_mbulk[subclass[i]][:,1][t], xmin, xmax, xbins)
            out_mbulk[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
            next!(p)
        end

        # Bulk geochemical
        h1 = Plots.plot(
            ylims=(0,ybins),
            yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
            yflip=true,
            xlabel="SiO2 [wt.%]", 
            ylabel="Age [Ma]",
            size=(600,500),
        )
        Plots.heatmap!(h1, out_bulk, 
            colorbar=false,
            color=c_gradient,
            title="$(subfig[i][1]). Bulk $(labels[i]) Data\n"
        )

        # Matched samples
        h2 = Plots.plot(
            ylims=(0,ybins),
            yticks=false,
            yflip=true,
            xlabel="SiO2 [wt.%]", 
            size=(600,500),
        )
        Plots.heatmap!(h2, out_mbulk, 
            colorbar=false,
            color=c_gradient,
            title="$(subfig[i][2]). Matched $(labels[i]) Samples\n"
        )

        # Assemble all plots with a common color bar
        l = @layout [a{0.95w} b{0.05w}]
        a = Plots.plot(h1, h2, layout=(1,2),
            size=(1000, 400),
            framestyle=:box,
            grid=false,
            fontfamily=:Helvetica,
            xticks=(0:xbins/4:xbins, string.(40:10:80)),
            # xticks=(0:60:240, string.(40:10:80)),
            # xticks=false,
            left_margin=(25,:px), bottom_margin=(15,:px),
            titleloc=:left, titlefont = font(12),
        )
        b = Plots.heatmap(rand(2,2), clims=(0,1), 
            framestyle=:none, color=c_gradient, colorbar_title="Relative Sample Density", 
            lims=(-1,0)
        )

        h = Plots.plot(a, b, layout=l)
        fig[i] = h
        savefig(h, "$filepath/silica_heatmap_$(subclass[i]).pdf")
    end

    # Assemble The Big Plot™
    # Actually, probably do this in illustrator because this keeps moving stuff around
    h = Plots.plot(fig..., layout=(3, 1), size=(1000,1300))

    display(h)
    savefig("$filepath/silica_heatmap_subclasses.pdf")


## --- End of file 