## --- Set up
    # Where did the Daly Gap go in matched sample set :(

    # Packages 
    using RockWeatheringFlux
    using HDF5, DelimitedFiles, KernelDensity
    using Plots

    # Preallocate / Local definitions
    nsims = Int(1e7)                         # 10 M simulations
    SiO2_error = 1.0                         # SiO₂ error
    age_error = 0.05                         # Default age error

    # Echo info 
    @info """ Files:
    Geochemical data: $geochem_fid
    Lithologic data:  $macrostrat_io
    """

## --- Load data
    # Indices of matched samples
    fid = readdlm(matchedbulk_io)
    matches = Int.(vec(fid[:,1]))
    t = @. matches != 0;

    # Assigned rock class during sample matching 
    match_class = string.(vec(fid[:,2]))[t];
    match_cats = match_rocktype(match_class);

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
        rocktype = read(fid["vars"]["rocktype"])[t],
        rockname = read(fid["vars"]["rockname"])[t],
        rockdescrip = read(fid["vars"]["rockdescrip"])[t],
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
    include_minor!(macro_cats);
    include_minor!(bulk_cats);
    macro_cats = delete_cover(macro_cats);
    bulk_cats = delete_cover(bulk_cats);
    

## --- Resample bulk geochemical igneous / volcanic / plutonic 
    # Preallocate 
    simbulk = NamedTuple{(:ign, :plut, :volc)}(Array{Float64}(undef, nsims, 2) for _ in 1:3)

    # Compute age uncertainties
    ageuncert = nanadd.(bulk.Age_Max, .- bulk.Age_Min) ./ 2;
    ageuncert[isnan.(bulk.Age) .| (ageuncert .== 0)] .= NaN;
    for i in eachindex(ageuncert)
        # Default 5% age uncertainty if bounds do not exist (or error is 0)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), bulk.Age[i]*age_error, ageuncert[i])
    end

    # Restrict to only samples with data and reset any weirdness with types
    t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude) & !isnan(bulk.Age)
    include_minor!(bulk_cats)

    # Bulk igneous
    s = t .& bulk_cats.ign
    k = invweight(bulk.Latitude[s], bulk.Longitude[s], bulk.Age[s])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2[s] bulk.Age[s]]
    uncertainty = [fill(SiO2_error, count(s)) ageuncert[s]]
    simbulk.ign .= bsresample(data, uncertainty, nsims, p)

    # Bulk plutonic
    s = t .& bulk_cats.volc
    k = invweight(bulk.Latitude[s], bulk.Longitude[s], bulk.Age[s])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2[s] bulk.Age[s]]
    uncertainty = [fill(SiO2_error, count(s)) ageuncert[s]]
    simbulk.plut .= bsresample(data, uncertainty, nsims, p)

    # Bulk volcanic 
    s = t .& bulk_cats.plut
    k = invweight(bulk.Latitude[s], bulk.Longitude[s], bulk.Age[s])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2[s] bulk.Age[s]]
    uncertainty = [fill(SiO2_error, count(s)) ageuncert[s]]
    simbulk.volc .= bsresample(data, uncertainty, nsims, p)


## --- Resample matched dataset igneous / volcanic / plutonic 
    # Preallocate 
    target = (:ign, :plut, :volc, :basalt, :rhyolite)
    sim_mbulk = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

    # Compute age uncertainties
    ageuncert = nanadd.(macrostrat.agemax, - macrostrat.agemin) ./ 2;
    ageuncert[isnan.(macrostrat.age) .| (ageuncert .== 0)] .= NaN;
    for i in eachindex(ageuncert)
        # Default 5% age uncertainty if bounds do not exist (or error is 0)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), macrostrat.age[i]*age_error, ageuncert[i])
    end
    
    # Restrict to only samples with data (all samples have lat / lon) and reset any 
    # weirdness with types
    t = @. !isnan(macrostrat.age);
    include_minor!(macro_cats);

    for key in target 
        s = t .& macro_cats[key]
        k = invweight_age(macrostrat.age[s])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        data = [mbulk.SiO2[s] macrostrat.age[s]]
        uncertainty = [fill(SiO2_error, count(s)) ageuncert[s]]
        sim_mbulk[key] .= bsresample(data, uncertainty, nsims, p)
    end


## --- Resample observed rock types, correcting for preservation bias
    # We can get a completely spatially representative dataset of the rock types exposed
    # at Earth's surface from the lithologic map, although the dataset is subject to 
    # preservation bias.

    # Get rock classes
    include_minor!(macro_cats)
    a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    for i in eachindex(keys(macro_cats))
        for j in eachindex(macro_cats[i])
            a[j,i] = ifelse(macro_cats[i][j], 1, 0)
        end
    end

    # Compute age uncertainties
    ageuncert = nanadd.(macrostrat.agemax, - macrostrat.agemin) ./ 2;
    ageuncert[isnan.(macrostrat.age) .| (ageuncert .== 0)] .= NaN;
    for i in eachindex(ageuncert)
        # Default 5% age uncertainty if bounds do not exist (or error is 0)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), macrostrat.age[i]*age_error, ageuncert[i])
    end

    # Resample with temporal weights 
    k = invweight_age(macrostrat.age)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [a macrostrat.age]
    uncertainty = [zeros(size(a)) ageuncert]
    simout = bsresample(data, uncertainty, nsims, p)

    # Parse output
    cats = simout[:,1:end-1] .> 0
    lith_cats = NamedTuple{keys(macro_cats)}([cats[:,i] for i in eachindex(keys(macro_cats))])
    sim_lith_age = simout[:,end]
    sim_lith_age[.!(0 .< sim_lith_age .< 4000)] .= NaN


## --- Resample observed rock types, WITHOUT correcting for preservation bias
    # If we're just interested in what's on the surface **right now** we don't want to 
    # correct for preservation bias

    # Get rock classes
    include_minor!(macro_cats)
    a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    for i in eachindex(keys(macro_cats))
        for j in eachindex(macro_cats[i])
            a[j,i] = ifelse(macro_cats[i][j], 1, 0)
        end
    end

    # Compute age uncertainties
    ageuncert = nanadd.(macrostrat.agemax, - macrostrat.agemin) ./ 2;
    ageuncert[isnan.(macrostrat.age) .| (ageuncert .== 0)] .= NaN;
    for i in eachindex(ageuncert)
        # Default 5% age uncertainty if bounds do not exist (or error is 0)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), macrostrat.age[i]*age_error, ageuncert[i])
    end

    # Resample to propogate age uncertainty
    data = [a macrostrat.age]
    uncertainty = [zeros(size(a)) ageuncert]
    simout = bsresample(data, uncertainty, nsims, ones(length(ageuncert)))

    # Parse output
    cats = simout[:,1:end-1] .> 0
    lith_cats_unweighted = NamedTuple{keys(macro_cats)}([cats[:,i] for i in eachindex(keys(macro_cats))])
    sim_lith_age_unweighted = simout[:,end]
    sim_lith_age_unweighted[.!(0 .< sim_lith_age_unweighted .< 4000)] .= NaN


## --- Resample silica and age distributions in the bulk geochemical datset 
    # We can get the distributions of silica in different rock types from the geochemical
    # dataset. This is subject to sampling bias as well as preservation bias. 
    # Spatiotemporally resample the dataset to remove sampling and preservation bias. 
    # Tack the rock types on too, so we can get them later when we correlate this with 
    # spatial abundance of each rock class

    include_minor!(bulk_cats)
    a = Array{Int64}(undef, length(bulk_cats[1]), length(bulk_cats))
    for i in eachindex(keys(bulk_cats))
        for j in eachindex(bulk_cats[i])
            a[j,i] = ifelse(bulk_cats[i][j], 1, 0)
        end
    end

    # Compute age uncertainty 
    ageuncert = nanadd.(bulk.Age_Max, .- bulk.Age_Min) ./ 2;
    ageuncert[isnan.(bulk.Age) .| (ageuncert .== 0)] .= NaN;
    for i in eachindex(ageuncert)
        # Default 5% age uncertainty if bounds do not exist (or error is 0)
        ageuncert[i] = ifelse(isnan(ageuncert[i]), bulk.Age[i]*age_error, ageuncert[i])
    end

    # Resample with spatemporal weights 
    k = invweight(bulk.Latitude, bulk.Longitude, bulk.Age)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    data = [bulk.SiO2 bulk.Age a]
    uncertainty = [fill(SiO2_error, length(bulk.SiO2)) ageuncert zeros(size(a))]
    simout = bsresample(data, uncertainty, nsims, p)

    # Parse output
    sim_silica = simout[:,1]
    sim_chem_age = simout[:,2]
    sim_chem_age[.!(0 .< sim_chem_age .< 4000)] .= NaN
    cats = simout[:,3:end] .> 0
    chem_cats = NamedTuple{keys(bulk_cats)}([cats[:,i] for i in eachindex(keys(bulk_cats))])


## --- 2D histograms of age and silica, but divided by volcanic / plutonic classes
    # Preallocate 
    xmin, xmax, xbins = 40, 80, 380          # Silica
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380         # Age
    yedges = ymin:(ymax-ymin)/ybins:ymax

    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    class = (:plut, :volc, :ign)

    # Construct plots
    p = Progress(6*ybins, desc="Making plots:");
    for i in eachindex(fig)
        # Preallocate
        target = bulk_cats[class[i]]
        out_bulk = zeros(ybins, xbins)
        out_mbulk = zeros(ybins, xbins)

        # Bulk geochemcial data
        for j = 1:ybins
            t = @. yedges[j] <= simbulk[class[i]][:,2] < yedges[j+1]        # Get this age bin
            c, n = bincounts(simbulk[class[i]][:,1][t], xmin, xmax, xbins)  # Count
            n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))     # Normalize
            out_bulk[j,:] .= n                                              # Put in output array
            next!(p)
        end

        # Matched samples
        for j = 1:ybins
            t = @. yedges[j] <= sim_mbulk[class[i]][:,2] < yedges[j+1]
            c, n = bincounts(sim_mbulk[class[i]][:,1][t], xmin, xmax, xbins)
            n = (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
            out_mbulk[j,:] .= n
            next!(p)
        end

        # Plot 
            h1 = Plots.plot(
                ylims=(0,ybins),
                yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
                yflip=true,
                xlabel="SiO₂ [wt.%]", 
                ylabel="Age [Ma]",
                size=(600,500),
            )
            Plots.heatmap!(h1, out_bulk, 
                colorbar=false,
                color=:nipy_spectral,
                title="A. Resampled $geochem_fid Samples\n"
            )

            # Matched samples
            h2 = Plots.plot(
                ylims=(0,ybins),
                yticks=false,
                yflip=true,
                xlabel="SiO₂ [wt.%]", 
                size=(600,500),
            )
            Plots.heatmap!(h2, out_mbulk, 
                colorbar=false,
                color=:nipy_spectral,
                title="B. Resampled Matched $(class[i]) Samples\n"
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
            # display(h)

        fig[i] = h
    end

    # Display plots here so the error message doesn't mess up the progress bar
    display(fig[1])
    display(fig[2])
    display(fig[3])


## --- Predicted age / silica distribution with 100 Ma year bins
    # We should be able to construct an expected silica distribution and time histogram
    # from the relative proportion of rock types on Earth's surface and the distributions
    # of silica for each rock class, given we know the ages of both datasets 

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


## --- Visualize component distributions of volc predicted histogram 
    # Basically, instead of a composite histogram, I want to plot all the distributions of 
    # igneous rock types, weighted by their surficial abundance, on top of each other  
    minorvolc, minorplut, minorign = get_rock_class()[3:5];
    prop = float.([count(macro_cats[k]) for k in minorvolc])
    prop ./= nansum(prop)
    here = @. !isnan(mbulk.SiO2) .& (macrostrat.age .> 1500);
    macro_cats.basalt .&= .!(macro_cats.rhyolite .| macro_cats.granite .| macro_cats.siliciclast);
    macro_cats.rhyolite .&= .!(macro_cats.shale .| macro_cats.andesite);

    # Assemble plots
    h = plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica, 
        xlims=(40,80),
        xticks=(40:10:80, string.(40:10:80)),
        yticks=false, 
        xlabel="SiO2 [wt.%]", ylabel="Relative Abundance",
        legend=:outerright,
        size=(800, 500), 
        left_margin=(20,:px)
    );
    for i in eachindex(minorvolc)
        u = kde(mbulk.SiO2[macro_cats[minorvolc[i]] .& here])
        u.density .*= prop[i]
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[minorvolc[i]], linecolor=colors[minorvolc[i]], 
            label="$(minorvolc[i])"
        )
    end
    display(h)
    
    # # Do the same thing for the bulk geochemical dataset 
    # here = @. !isnan(bulk.SiO2) .& (bulk.Age .> 1500);
    # h = plot(
    #     framestyle=:box, 
    #     grid = false,
    #     fontfamily=:Helvetica, 
    #     xlims=(40,80),
    #     xticks=(40:10:80, string.(40:10:80)),
    #     yticks=false, 
    #     xlabel="SiO2 [wt.%]", ylabel="Relative Abundance",
    #     legend=:outerright,
    #     size=(800, 500), 
    #     left_margin=(20,:px)
    # );
    # for i in eachindex(minorvolc)
    #     u = kde(bulk.SiO2[bulk_cats[minorvolc[i]] .& here])
    #     u.density .*= prop[i]
    #     Plots.plot!(u.x, u.density, 
    #         seriestype=:path, linewidth=4,
    #         color=colors[minorvolc[i]], linecolor=colors[minorvolc[i]], 
    #         label="$(minorvolc[i])"
    #     )
    # end
    # display(h)


## --- Contamination of basalts and rhyolites 
    here = @. !isnan(mbulk.SiO2) .& (macrostrat.age .> 1500);

    # We know that if we use the assigned types from the sample matching file (one type
    # per sample, no major classes), our distributions look like we expect them to.
    h = plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica, 
        xlims=(40,80),
        xticks=(40:10:80, string.(40:10:80)),
        yticks=false, 
        xlabel="SiO2 [wt.%]", ylabel="Relative Abundance",
        legend=:outerright,
        size=(800, 500), 
        left_margin=(20,:px)
    );
    for i in eachindex(minorvolc)
        u = kde(mbulk.SiO2[match_cats[minorvolc[i]] .& here])
        u.density .*= prop[i]
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[minorvolc[i]], linecolor=colors[minorvolc[i]], 
            label="$(minorvolc[i])"
        )
    end
    display(h)

    # # Get the BitVectors to where they would have been before type assignment 
    # exclude_minor!(macro_cats)

    # # Make metamorphic rocks only those without a protolith 
    # undiff_met = copy(macro_cats.met)
    # for type in keys(macro_cats)
    #     type==:met && continue
    #     macro_cats.met .&= .!macro_cats[type]
    # end 

    # Problematic basalts
    b = (mbulk.SiO2 .> 60) .& (macrostrat.age .> 1500) .& macro_cats.basalt .& .!match_cats.basalt;
    count(b)        # 359
    unique(macrostrat.rocktype[b])
    unique([match_class[b] macrostrat.rocktype[b]], dims=1)

    # Problematic rhyolites
    r = (mbulk.SiO2 .< 68) .& (macrostrat.age .> 1500) .& macro_cats.rhyolite .& .!match_cats.rhyolite;
    count(r)        # 466
    unique([match_class[r] macrostrat.rocktype[r]], dims=1)

    # Hypothesis: this is from the matched rock types that are matched with basalt but 
    # were assigned to a different rock type during sample matching 

    # Plot matched sample classes that make up basalts, weighted by their presence in the 
    # match_type dataset
    c = countmap(match_class[b]);
    c = Dict([k => c[k]/sum(values(c)) for k in keys(c)])

    h = plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica, 
        xlims=(40,80),
        xticks=(40:10:80, string.(40:10:80)),
        yticks=false, 
        xlabel="SiO2 [wt.%]", ylabel="Relative Abundance",
        legend=:outerright,
        size=(800, 500), 
        left_margin=(20,:px)
    );
    for k in Symbol.(keys(c))
        s = match_cats[k] .& macro_cats.rhyolite .& here
        count(s)==0 && continue
        u = kde(mbulk.SiO2[s])
        u.density .*= c[string(k)]
        Plots.plot!(u.x, u.density, 
            seriestype=:path, linewidth=4,
            color=colors[k], linecolor=colors[k], 
            label="$k"
        )
    end
    display(h)


## --- Contamination of basalts and rhyolites in 2D histogram form 
    # Without resampling because I want things fast 
    xmin, xmax, xbins = 40, 80, 240           # Silica bins
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380          # Age bins
    yedges = ymin:(ymax-ymin)/ybins:ymax

    # Bulk geochemical data
    out = zeros(ybins, xbins)
    for j in 1:ybins
        t = @. yedges[j] <= macrostrat.age < yedges[j+1];
        c, n = bincounts(sim_silica[chem_cats.basalt .& t], xmin, xmax, xbins)
        out[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    end

    h = Plots.plot(
        ylims=(0,ybins),
        yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
        xticks=(0:xbins/4:xbins, string.(40:10:80)),
        yflip=true,
        xlabel="SiO₂ [wt.%]", 
        ylabel="Age [Ma]",
        size=(600,500),
    );
    Plots.heatmap!(h, out, color=c_gradient)
    display(h)


## --- Options to mitigate contamination and why I dislike them 
    # 1) only use the assigned rock classes 
        # IDK bad vibes. I feel like if everything works right we shouldn't need to use 
        # this??
        #
        # Maybe some kind of compromise where we exclude e.g. seds... like rather than doing
        # just igneous rocks we could do *only* igneous rocks for this kind of analysis
        # 
        # On the other hand, does the system really work if it fucks up the distributions?

    # 2) Screen rock classes for silica distributions that make sense (e.g. something with
    # very high silica is unlikely to be a basalt, so the likelihood of deleting it from the
    # basalts dataset is proportional to the distance from the expected basalt distribution)
        # This feels too much like forcing the data to be what I want it to be


## --- Is the proportion of sed / met / ign constant through time?
    # Let's do sed / volc / plut / met every 100 Ma.
    ymin, ymax, ybins = 0, 3800, 38         # Age
    yedges = ymin:(ymax-ymin)/ybins:ymax
    
    # Preallocate 
    p_sed = Array{Float64}(undef, ybins)
    p_volc = Array{Float64}(undef, ybins)
    p_plut = Array{Float64}(undef, ybins)
    p_met = Array{Float64}(undef, ybins)

    # Count the proportion of major classes that are undifferentiated
    lith_cats_incl = deepcopy(lith_cats) 
    include_minor!(lith_cats_incl) 
    exclude_minor!(lith_cats) 
    for i = 1:ybins
        t = @. yedges[i] <= sim_lith_age < yedges[i+1]      # Number of possible samples

        p_sed[i] = count(lith_cats.sed .& t) / count(lith_cats_incl.sed .& t)
        p_volc[i] = count(lith_cats.volc .& t) / count(lith_cats_incl.volc .& t)
        p_plut[i] = count(lith_cats.plut .& t) / count(lith_cats_incl.plut .& t)
        p_met[i] = count(lith_cats.met .& t) / count(t)
    end

    h = Plots.plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica, 
        xlims=(0,3800),
        legend=:outerright,
        size=(800, 500), 
        left_margin=(20,:px)
    );
    Plots.plot!(h, cntr(yedges), p_sed, markershape=:circle, color=colors.sed, label="sed")
    Plots.plot!(h, cntr(yedges), p_volc, markershape=:circle, color=colors.volc, label="volc")
    Plots.plot!(h, cntr(yedges), p_plut, markershape=:circle, color=colors.plut, label="plut")
    Plots.plot!(h, cntr(yedges), p_met, markershape=:circle, color=colors.met, label="met")
    vline!([1500], label="1500 Ma.")


## --- Proportion of volcanic rocks made of basalt / rhyolite over time, 100 Myr. bins
    ymin, ymax, ybins = 0, 3800, 38
    yedges = ymin:(ymax-ymin)/ybins:ymax

    p_basalt = Array{Float64}(undef, ybins)
    p_rhyolite = Array{Float64}(undef, ybins)

    include_minor!(lith_cats)
    for i = 1:ybins
        t = @. yedges[i] <= sim_lith_age < yedges[i+1]
        p_basalt[i] = count(lith_cats.basalt .& t) / count(lith_cats.volc .& t)
        p_rhyolite[i] = count(lith_cats.rhyolite .& t) / count(lith_cats.volc .& t)
    end

    h = Plots.plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica, 
        xlims=(0,3800),
        legend=:topright,
        left_margin=(20,:px),
        xlabel="Age [Ma.]", ylabel="Relative Proportion"
    );
    Plots.plot!(h, cntr(yedges), p_basalt, markershape=:circle, color=colors.basalt,
        label="basalt"
    )
    Plots.plot!(h, cntr(yedges), p_rhyolite, markershape=:circle, color=colors.rhyolite, 
        label="rhyolite"
    )
    display(h)


## --- 2D histogram of basaltic silica over time
    xmin, xmax, xbins = 40, 80, 240           # Silica bins
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380          # Age bins
    yedges = ymin:(ymax-ymin)/ybins:ymax

    # Bulk geochemical data
    out = zeros(ybins, xbins)
    for j in 1:ybins
        t = @. yedges[j] <= sim_chem_age < yedges[j+1];
        c, n = bincounts(sim_silica[chem_cats.basalt .& t], xmin, xmax, xbins)
        out[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    end

    h = Plots.plot(
        ylims=(0,ybins),
        yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
        xticks=(0:xbins/4:xbins, string.(40:10:80)),
        yflip=true,
        xlabel="SiO₂ [wt.%]", 
        ylabel="Age [Ma]",
        size=(600,500),
    );
    Plots.heatmap!(h, out, color=c_gradient)
    display(h)

    # Matched samples
    out = zeros(ybins, xbins)
    for j in 1:ybins
        t = @. yedges[j] <= sim_mbulk.basalt[:,2] < yedges[j+1];
        c, n = bincounts(sim_mbulk.basalt[:,1][t], xmin, xmax, xbins)
        out[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    end

    h = Plots.plot(
        ylims=(0,ybins),
        yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
        xticks=(0:xbins/4:xbins, string.(40:10:80)),
        yflip=true,
        xlabel="SiO₂ [wt.%]", 
        ylabel="Age [Ma]",
        size=(600,500),
    );
    Plots.heatmap!(h, out, color=c_gradient,)
    display(h)


## --- Same as above, with rhyolites
    xmin, xmax, xbins = 40, 80, 240           # Silica bins
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380          # Age bins
    yedges = ymin:(ymax-ymin)/ybins:ymax

    # Bulk geochemical data
    out = zeros(ybins, xbins)
    for j in 1:ybins
        t = @. yedges[j] <= sim_chem_age < yedges[j+1];
        c, n = bincounts(sim_silica[chem_cats.rhyolite .& t], xmin, xmax, xbins)
        out[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    end

    h = Plots.plot(
        ylims=(0,ybins),
        yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
        xticks=(0:xbins/4:xbins, string.(40:10:80)),
        yflip=true,
        xlabel="SiO₂ [wt.%]", 
        ylabel="Age [Ma]",
        size=(600,500),
    );
    Plots.heatmap!(h, out, color=c_gradient)
    display(h)

    # Matched samples
    out = zeros(ybins, xbins)
    for j in 1:ybins
        t = @. yedges[j] <= sim_mbulk.rhyolite[:,2] < yedges[j+1];
        c, n = bincounts(sim_mbulk.rhyolite[:,1][t], xmin, xmax, xbins)
        out[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))
    end

    h = Plots.plot(
        ylims=(0,ybins),
        yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
        xticks=(0:xbins/4:xbins, string.(40:10:80)),
        yflip=true,
        xlabel="SiO₂ [wt.%]", 
        ylabel="Age [Ma]",
        size=(600,500),
    );
    Plots.heatmap!(h, out, 
        color=:nipy_spectral,
    )
    display(h)


## --- End of File 