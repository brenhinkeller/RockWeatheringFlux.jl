## --- Set up
    # Where did the Daly Gap go in matched sample set :(

    # Packages 
    using RockWeatheringFlux
    using HDF5, DelimitedFiles
    using Plots

    # Preallocate / Local definitions
    nsims = Int(1e7)                         # 10 M simulations
    SiO2_error = 5.0                         # Large SiO₂ error to smooth data
    age_min_error = 0.15                     # Minimum 15% age error

    xmin, xmax, xbins = 40, 80, 240          # Silica
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380         # Age
    yedges = ymin:(ymax-ymin)/ybins:ymax
    
    # Echo info 
    @info """ Files:
    Geochemical data: $geochem_fid
    Lithologic data:  $macrostrat_io
    """

## --- Load data
    # Indices of matched samples
    fid = readdlm(matchedbulk_io)
    matches = Int.(vec(fid[:,1]))
    t = @. matches != 0

    # Geochemical Data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][matches[t]] 
        for i in eachindex(header)])
    close(fid)

    # Macrostrat
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    # Unmatched geochemcial data 
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    # Major classes include minors, delete cover
    include_minor!(macro_cats)
    include_minor!(bulk_cats)
    macro_cats = delete_cover(macro_cats)
    bulk_cats = delete_cover(bulk_cats)
    

## --- 2D histograms of age and silica, but divided by volcanic / plutonic classes
    # Preallocate 
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    class = (:plut, :volc, :ign)

    # Construct plots
    for i in eachindex(fig)
        # Preallocate
        target = bulk_cats[class[i]]
        out_bulk = zeros(ybins, xbins)
        out_mbulk = zeros(ybins, xbins)

        # EarthChem 
        t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude) & !isnan(bulk.Age) & bulk_cats[class[i]]
        ageuncert = bulk.Age .* age_min_error
        for j in eachindex(ageuncert)
            calc_uncert = (bulk.Age_Max[j] .- bulk.Age_Min[j])/2
            ageuncert[j] = nanmax(calc_uncert, ageuncert[j])
        end

        k = invweight(bulk.Latitude[t], bulk.Longitude[t], bulk.Age[t])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        data = [bulk.SiO2[t] bulk.Age[t]]
        uncertainty = [fill(SiO2_error, count(t)) ageuncert[t]]
        simbulk = bsresample(data, uncertainty, nsims, p)

        for j = 1:ybins
            t = @. yedges[j] <= simbulk[:,2] < yedges[j+1]                # Get this age bin
            c, n = bincounts(simbulk[:,1][t], xmin, xmax, xbins)          # Count
            n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))   # Normalize
            out_bulk[j,:] .= n                                            # Put in output array
        end

        # Matched samples 
        t = @. !isnan(macrostrat.rocklat) & !isnan(macrostrat.rocklon) & macro_cats[class[i]];
        t .&= (.!isnan.(mbulk.Age) .| .!isnan.(macrostrat.age));

        age = copy(mbulk.Age)
        age[isnan.(age)] .= macrostrat.age[isnan.(age)]
        ageuncert = age .* age_min_error
        for i in eachindex(ageuncert)
            calcuncert = abs(nanadd(mbulk.Age_Max[i], -mbulk.Age_Min[i])/2)
            ageuncert[i] = nanmax(calcuncert, ageuncert[i])

        end

        k = invweight_age(age[t])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        data = [mbulk.SiO2[t] age[t]]
        uncertainty = [fill(SiO2_error, count(t)) ageuncert[t]]
        sim_mbulk = bsresample(data, uncertainty, nsims, p)

        for j = 1:ybins
            t = @. yedges[j] <= sim_mbulk[:,2] < yedges[j+1]
            c, n = bincounts(sim_mbulk[:,1][t], xmin, xmax, xbins)
            n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
            out_mbulk[j,:] .= n
        end

        # Plot 
        h1 = Plots.plot(
            ylims=(0,380),
            yticks=(0:50:380, string.(0:500:3800)),
            yflip=true,
            xlabel="SiO₂ [wt.%]", 
            ylabel="Age [Ma]",
            size=(600,500),
        )
        Plots.heatmap!(h1, out_bulk, 
            colorbar=false,
            color=c_gradient,
            title="A. Resampled $geochem_fid Samples\n"
        )

        # Matched samples
        h2 = Plots.plot(
            ylims=(0,380),
            yticks=false,
            yflip=true,
            xlabel="SiO₂ [wt.%]", 
            size=(600,500),
        )
        Plots.heatmap!(h2, out_mbulk, 
            colorbar=false,
            color=c_gradient,
            title="B. Resampled Matched $(class[i]) Samples\n"
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

        fig[i] = h
    end


## --- Construct a predicted histogram over 100 Ma year bins
    # We can get a completely spatially representative dataset of the rock types exposed
    # at Earth's surface from the lithologic map, although the dataset is subject to 
    # preservation bias. Temporally resample to correct for this.
    a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    for i in eachindex(keys(macro_cats))
        for j in eachindex(macro_cats[i])
            a[j,i] = ifelse(macro_cats[i][j], 1, 0)
        end
    end
    
    ageuncert = nanadd.(macrostrat.agemax, .- macrostrat.agemin) ./ 2;
    ageuncert[isnan.(macrostrat.age) .| (ageuncert .== 0)] .= NaN;
    for i in eachindex(ageuncert)
        # Default 5% age uncertainty if bounds do not exist (or error is 0)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), macrostrat.age[i]*0.05 , ageuncert[i])
    end

    k = invweight_age(macrostrat.age)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [a macrostrat.age]
    uncertainty = [zeros(size(a)) ageuncert]
    simout = bsresample(data, uncertainty, nsims, p)

    cats = simout[:,1:end-1] .> 0
    lith_cats = NamedTuple{keys(macro_cats)}([cats[:,i] for i in eachindex(keys(macro_cats))])
    sim_lith_age = simout[:,end]
    sim_lith_age[.!(0 .< sim_lith_age .< 4000)] .= NaN

    # Alternatively, don't correct for preservation bias--we're not technically interested
    # in the change over time, but in what's exposed **right now**. Resample to propogate
    # error, but don't weight that resample
    data = [a macrostrat.age]
    uncertainty = [zeros(size(a)) ageuncert]
    simout = bsresample(data, uncertainty, nsims, ones(length(ageuncert)))

    cats = simout[:,1:end-1] .> 0
    lith_cats_unweighted = NamedTuple{keys(macro_cats)}([cats[:,i] for i in eachindex(keys(macro_cats))])
    sim_lith_age_unweighted = simout[:,end]
    sim_lith_age_unweighted[.!(0 .< sim_lith_age_unweighted .< 4000)] .= NaN

    # We can get the distributions of silica in different rock types from the geochemical
    # dataset. This is subject to sampling bias as well as preservation bias. 
    # Spatiotemporally resample the dataset to remove sampling and preservation bias. 
    # Tack the rock types on too, so we can get them later when we correlate this with 
    # spatial abundance of each rock class
    a = Array{Int64}(undef, length(bulk_cats[1]), length(bulk_cats))
    for i in eachindex(keys(bulk_cats))
        for j in eachindex(bulk_cats[i])
            a[j,i] = ifelse(bulk_cats[i][j], 1, 0)
        end
    end
    
    ageuncert = nanadd.(bulk.Age_Max, .- bulk.Age_Min) ./ 2;
    ageuncert[isnan.(bulk.Age) .| (ageuncert .== 0)] .= NaN;
    for i in eachindex(ageuncert)
        # Default 5% age uncertainty if bounds do not exist (or error is 0)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), bulk.Age[i]*0.05 , ageuncert[i])
    end

    k = invweight(bulk.Latitude, bulk.Longitude, bulk.Age)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2 bulk.Age a]
    uncertainty = [fill(1.0, length(bulk.SiO2)) ageuncert zeros(size(a))]     # 1 wt.% SiO₂ error
    simout = bsresample(data, uncertainty, nsims, p)

    sim_silica = simout[:,1]
    sim_chem_age = simout[:,2]
    sim_chem_age[.!(0 .< sim_chem_age .< 4000)] .= NaN
    cats = simout[:,3:end] .> 0
    chem_cats = NamedTuple{keys(macro_cats)}([cats[:,i] for i in eachindex(keys(macro_cats))])

    # We should be able to construct an expected silica distribution and time histogram
    # from these two things... 
    xmin, xmax, xbins = 40, 80, 80          # Silica bins
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 38         # Age bins
    yedges = ymin:(ymax-ymin)/ybins:ymax

    # For each rock class (volcanic / plutonic / igneous), for each time bin, compute a 
    # silica distribution 
    minorvolc, minorplut, minorign = get_rock_class()[3:5];
    class = (:plut, :volc, :ign)
    minorclass = (minorplut, minorvolc, (minorplut..., minorvolc..., :carbonatite))
    for i in eachindex(class) 
        out = zeros(ybins, xbins)
        for j in 1:ybins
            # # Get relative abundances of rock types exposed on the surface for this age bin
            # t = @. yedges[j] <= sim_lith_age < yedges[j+1];
            # prop = float.([count(lith_cats[k] .& t) for k in minorclass[i]])
            # prop ./= nansum(prop)

            # Alternatively, don't correct for preservation bias in exposed rocks!
            t = @. yedges[j] <= sim_lith_age_unweighted < yedges[j+1];
            prop = float.([count(lith_cats_unweighted[k] .& t) for k in minorclass[i]])
            prop ./= nansum(prop)

            # Get the cumulative distribution of silica in this age bin, constructed of 
            # the distribution of each minor class weighted by it's surficial abundance 
            s = @. yedges[j] <= sim_chem_age < yedges[j+1];
            dist = zeros(xbins)
            for k in eachindex(prop)
                c, n = bincounts(sim_silica[chem_cats[minorclass[i]][k] .& s], xmin, xmax, xbins)
                n = float(n) ./ nansum(float(n) .* step(c))
                dist .+= (n*prop[k])
            end

            # Enter into the output array
            out[j,:] .= dist
        end

        # Plot 
        h = Plots.plot(
            ylims=(0,38),
            yticks=(0:5:38, string.(0:500:3800)),
            yflip=true,
            xlabel="SiO₂ [wt.%]", 
            ylabel="Age [Ma]",
            size=(600,500),
        )
        Plots.heatmap!(h, out, 
            # colorbar=false,
            color=c_gradient,
            title="$(class[i])"
        )
        display(h)
    end


## --- End of File 