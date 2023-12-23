## --- Set up 
    # 2D histograms of igneous and metamorphic samples as a function of age and silica 
    # content. After Figure 6.9 from Keller, 2016 (10.31237/osf.io/q7yra)

    # Load data and base packages
    if !@isdefined(filepath)
        include("Definitions.jl")
    end

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

        # Put the output in the output array
        out_bulk[i,:] .= n
    end


## --- Matched samples
    notsed = macro_cats.ign .| macro_cats.met
    t = @. !isnan(mbulk.Latitude) & !isnan(mbulk.Longitude) & !isnan(mbulk.Age) & notsed

    s = vec(isnan.(mbulk.Age))
    ageuncert = mbulk.Age .* age_min_error
    for i in eachindex(ageuncert)
        calc_uncert = (mbulk.Age_Max[i] .- mbulk.Age_Min[i])/2
        ageuncert[i] = ifelse(calc_uncert > ageuncert[i], calc_uncert, ageuncert[i])
    end

    k = invweight_age(mbulk.Age[t])     # Already spatially corrected
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [mbulk.SiO2[t] mbulk.Age[t]]
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
        title="A. Resampled EarthChem Samples\n"
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

    
## --- End of file 