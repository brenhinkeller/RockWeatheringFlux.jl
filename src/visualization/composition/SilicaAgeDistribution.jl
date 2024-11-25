## --- Set up 
    # 2D histograms of igneous and metamorphic samples as a function of age and silica 
    # content. After Figure 6.9 from Keller, 2016 (10.31237/osf.io/q7yra)

    # Load data and base packages
    include("Definitions.jl");

    # Preallocate / Local definitions
    # elem = :MgO                              # Element of interest
    elem = :SiO2
    elem_error = 1.0                         # Assumed element of interest (SiO₂) error
    age_error = 5                            # Minimum age error (%)
    age_error_abs = 50                       # Minimum age error (Ma)
    nsims = Int(1e7)                         # 10 M simulations

    xmin, xmax, xbins = 40, 80, 240          # Silica
    # xmin, xmax, xbins  = 0, 50, 300          # MgO
    xedges = xmin:(xmax-xmin)/xbins:xmax
    ymin, ymax, ybins = 0, 3800, 380         # Age
    yedges = ymin:(ymax-ymin)/ybins:ymax

    # Rock classes to resample
    target = (:ign, :plut, :volc)    

    # Set file 
    suffix = RockWeatheringFlux.version * "_" * RockWeatheringFlux.tag
    fpath = "src/visualization/composition/SilicaAgeDistribution_$(elem)_" * suffix * ".h5"

    # Preallocate 2D histogram data storage
    out_bulk = NamedTuple{target}(zeros(ybins, xbins) for _ in target)
    out_mbulk = NamedTuple{target}(zeros(ybins, xbins) for _ in target)


## --- Set switches 
    # [CHANGE ME] 
    # Set to TRUE to resample bulk geochemical or matched data 
    redo_bulk = true
    redo_matched = false

    # If the file doesn't exist, we gotta create it 
    !isfile(fpath) && (redo_bulk = redo_matched = true)

    # Typically, re-doing the 2D histogram is gonna depend on if we resampled or not 
    # If we resampled, we should redo the histogram. If we didn't resample, we don't have 
    # the data to redo the histogram so we gotta load from the file 
    redo_2D_bulk = redo_bulk 
    redo_2D_matched = redo_matched 

    # [CHANGE ME]
    # The exception to this is debugging. Say we resampled something, and we want to recalculate 
    # the 2D histogram for some reason using that data.
    redo_2D_bulk = true
    redo_2D_matched = false

    if (redo_bulk != redo_2D_bulk) || (redo_matched != redo_2D_matched)
        @warn "Your switches are in debugging mode"
    end


## --- Resample data 
    # Bulk geochemical data: spatiotemporal resample 
    if redo_bulk==true
        simbulk = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

        # Compute age uncertainties 
        ageuncert = nanadd.(bulk.Age_Max, .- bulk.Age_Min) ./ 2;
        for i in eachindex(ageuncert)
            ageuncert[i] = max(bulk.Age[i]*age_error, ageuncert[i], age_error_abs)
        end

        # Restrict to samples with data and resample 
        t = @. !isnan(bulk.Latitude) & !isnan(bulk.Longitude) & !isnan(bulk.Age);
        for key in target 
            s = t .& bulk_cats[key]
            k = invweight(bulk.Latitude[s], bulk.Longitude[s], bulk.Age[s])
            p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
            data = [bulk[elem][s] bulk.Age[s]]
            uncertainty = [fill(elem_error, count(s)) ageuncert[s]]
            simbulk[key] .= bsresample(data, uncertainty, nsims, p)
        end
    end

    # Matched data: temporal 
    if redo_matched==true
        sim_mbulk = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

        # Calculate sample age and uncertainty
        sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
            macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
            uncert_rel=age_error, uncert_abs=age_error_abs
        )

        # Restrict to only samples with data and resample 
        t = @. !isnan.(sampleage);
        for key in target 
            s = t .& match_cats[key]
            k = invweight_age(sampleage[s])
            p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
            data = [mbulk[elem][s] sampleage[s]]
            uncertainty = [fill(elem_error, count(s)) ageuncert[s]]
            sim_mbulk[key] .= bsresample(data, uncertainty, nsims, p)
        end
    end


## --- Calculate data for 2D histogram: sort matched data into time bins and normalize
    # All geochemical data 
    if redo_2D_bulk
        p = Progress((ybins % 10) * length(target), desc="Bulk 2D histogram...")

        for i in eachindex(target)
            out = zeros(ybins, xbins)
            for j = 1:ybins
                t = @. yedges[j] <= simbulk[target[i]][:,2] < yedges[j+1]
                c, n = bincounts(simbulk[target[i]][:,1][t], xmin, xmax, xbins)
                out[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))

                (i % 10 == 0) && next!(p)
            end
            out_bulk[target[i]] .= out
        end
    else
        fid = h5open(fpath, "r")
        out_bulk = NamedTuple{target}(read(fid["vars"]["resampled"]["$key"]) for key in target)
        close(fid)
    end

    # Matched geochemical data 
    if redo_2D_matched
        p = Progress((ybins % 10) * length(target), desc="Matched 2D histogram...")

        for i in eachindex(target)
            out = zeros(ybins, xbins)
            for j = 1:ybins
                t = @. yedges[j] <= sim_mbulk[target[i]][:,2] < yedges[j+1]
                c, n = bincounts(sim_mbulk[target[i]][:,1][t], xmin, xmax, xbins)
                out[j,:] .= (n .- nanminimum(n)) ./ (nanmaximum(n) - nanminimum(n))

                (i % 10 == 0) && next!(p)
            end
            out_mbulk[target[i]] .= out
            next!(p)
        end
    else
        fid = h5open(fpath, "r")
        out_mbulk = NamedTuple{target}(read(fid["vars"]["matched"]["$key"]) for key in target)
        close(fid)
    end

    # If we changed anything, just re-do the file
    if redo_2D_bulk || redo_2D_matched
        fid = h5open(fpath, "w")
        g = create_group(fid, "vars")

        g_resam = create_group(g, "resampled")
        for i in eachindex(target)
            write(g_resam, "$(target[i])", out_bulk[target[i]])
        end

        g_match = create_group(g, "matched")
        for i in eachindex(target)
            write(g_match, "$(target[i])", out_mbulk[target[i]])
        end

        close(fid)
    end


## --- Build plot (subdivide igneous into volcanic and plutonic classes)
    # Preallocate 
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, 3)
    subclass = (:ign, :volc, :plut)
    labels = ("Igneous", "Volcanic", "Plutonic")
    subfig = (("A", "B"), ("C", "D"), ("E", "F"))
    colorgrad = :inferno

    for i in eachindex(fig)
        # Base plot 
        h = Plots.plot(
            ylims=(0,ybins), yflip=true, 
            xticks=(0:xbins/4:xbins, string.(40:10:80)),
            yticks=(0:ybins/7.6:ybins, string.(0:500:3800)),
            xlabel="SiO2 [wt.%]", ylabel="Age [Ma]",
            size=(600,500),
            titleloc=:left, titlefont = font(12),
            tickfontsize=12, labelfontsize=14,
            fontfamily=:Helvetica,
        )

        # Resampled
        h1 = deepcopy(h)
        Plots.heatmap!(h1, out_bulk[subclass[i]],
            colorbar=false,
            color=colorgrad,
            title="$(subfig[i][1]). Resampled $(labels[i])\n",
        )
        savefig(h1, "$filepath/heatmap_$(subclass[i])_resampled.pdf")

        # Matched samples
        h2 = deepcopy(h)
        Plots.heatmap!(h2, out_mbulk[subclass[i]], 
            colorbar=false,
            color=colorgrad,
            title="$(subfig[i][2]). Matched $(labels[i])\n",
        )
        savefig(h2, "$filepath/heatmap_$(subclass[i])_matched.pdf")

        # Temporary display 
        h = plot(h1, h2, size=(1200,500))
        display(h)
        
        # # Assemble all plots with a common color bar
        # l = @layout [a{0.95w} b{0.05w}]
        # a = Plots.plot(h1, h2, layout=(1,2),
        #     size=(1000, 400),
        #     framestyle=:box,
        #     grid=false,
        #     fontfamily=:Helvetica,
        #     xticks=(0:xbins/4:xbins, string.(40:10:80)),
        #     # xticks=(0:60:240, string.(40:10:80)),
        #     # xticks=false,
        #     left_margin=(25,:px), bottom_margin=(15,:px),
        #     titleloc=:left, titlefont = font(12),
        # )
        # b = Plots.heatmap(rand(2,2), clims=(0,1), 
        #     framestyle=:none, color=colorgrad, 
        #     colorbar_title="Relative Sample Density", 
        #     lims=(-1,0)
        # )

        # h = Plots.plot(a, b, layout=l)
        # fig[i] = h
        # savefig(h, "$filepath/heatmap_$(subclass[i]).pdf")
    end

    # Save a colorbar 
    b = Plots.heatmap(rand(2,2), clims=(0,1), 
        framestyle=:none, color=colorgrad, 
        lims=(-1,0),
    )
    savefig(b, "$filepath/heatmap_colorbar.pdf")

    # # Assemble The Big Plot™
    # # But just for looks. Do this in illustrator because this keeps moving stuff around
    # h = Plots.plot(fig..., layout=(3, 1), size=(1000,1300))
    # savefig(h, "$filepath/heatmap_all_classes.pdf")
    # display(h)


## --- End of file 