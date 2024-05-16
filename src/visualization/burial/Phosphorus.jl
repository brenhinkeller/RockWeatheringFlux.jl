## --- Set up 
    # Potential correlation between the ratio of phosphorus to alkalininty and 
    # carbon burial 
    
    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots, StatsPlots, Colors

    # Save figures to: 
    filepath = "results/figures/burial"

    # Conversion factors from wt.% element oxide to moles per 100g sample
    # Note that these are molar masses and must be **divided** from the wt.% [g/g] value, 
    # hence the inversion (1 / value)
    CaO_to_Ca =   1 / (molarmass["Ca"]   + molarmass["O"]  )
    MgO_to_Mg =   1 / (molarmass["Mg"]   + molarmass["O"]  )
    K2O_to_K =    2 / (molarmass["K"] *2 + molarmass["O"]  )    # 2 mol K / 1 mol K₂O
    Na2O_to_Na =  2 / (molarmass["Na"]*2 + molarmass["O"]  )    # 2 mol Na / 1 mol Na₂O
    FeO_to_Fe =   1 / (molarmass["Fe"]   + molarmass["O"]  )
    P2O5_to_mol = 2 / (molarmass["P"] *2 + molarmass["O"]*5)    # 2 mol P / 1 mol P₂O₅


## --- Load data 
    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];
    classes = keys(match_cats)
    minor_classes = (minorsed..., minorvolc..., minorplut...,)

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)

    # Macrostrat and mapped rock classes
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        rocktype = read(fid["vars"]["rocktype"])[t],
        rockdescrip = read(fid["vars"]["rockdescrip"])[t],
        rockname = read(fid["vars"]["rockname"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    include_minor!(macro_cats)
    close(fid)


## --- Prepare data for resampling and calculate temporal weights 
    # Definitions 
    nsims = Int(1e6)
    xmin, xmax, nbins = 0, 3800, 38
    age_error = 0.05                   # Minimum age error (%)
    age_error_abs = 50                 # Minimum age error (Ma)
    data_uncert = 0.05                 # 5% uncertainty to all geochemical data

    # Elements in geochemical dataset 
    majors, minors = get_elements()
    allelements = Tuple([majors; minors])

    # Ages and uncertainties (prefer geochemical age unless value is missing)
    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
        uncert_rel=age_error, uncert_abs=age_error_abs
    )

    # Mole abundance of **CHARGE** for all major alkalinity-forming cations
    Ca²⁺ = @. mbulk.CaO * CaO_to_Ca .* 2
    Mg²⁺ = @. mbulk.MgO * MgO_to_Mg .* 2
    K⁺ = @. mbulk.K2O * K2O_to_K
    Na⁺ = @. mbulk.Na2O * Na2O_to_Na
    Fe²⁺ = @. mbulk.FeOT * FeO_to_Fe .* 2
    alkalinity = nansum([Ca²⁺ Mg²⁺ K⁺ Na⁺], dims=2)
    alkalinity_fe = nansum([Ca²⁺ Mg²⁺ K⁺ Na⁺ Fe²⁺], dims=2)

    # Mole abundance of phosphorus 
    phosphorus = mbulk.P2O5 .* P2O5_to_mol

    # Get all geochemical data into an array 
    geochem_in = Array{Float64}(undef, length(mbulk.SiO2), length(allelements))
    for i in eachindex(allelements)
        geochem_in[:,i] .= mbulk[allelements[i]]
    end

    # Rock class assigned during matching 
    matched_in = Array{Int64}(undef, length(match_cats.sed), length(classes))
    for i in eachindex(classes)
        for j in eachindex(match_cats[classes[i]])
            matched_in[j,i] = ifelse(match_cats[classes[i]][j], 1, 0)
        end
    end

    # Rock class as mapped
    mapped_in = Array{Int64}(undef, length(macro_cats.sed), length(classes))
    for i in eachindex(classes)
        for j in eachindex(macro_cats[classes[i]])
            mapped_in[j,i] = ifelse(macro_cats[classes[i]][j], 1, 0)
        end
    end

    # Create a filtered / groundtruthed rock class: all rock classes must be BOTH matched 
    # and mapped as that rock class, and metamorphic rocks are only undifferentiated 
    filter_cats = get_cats(false, length(match_cats.sed))[2]
    for k in classes
        filter_cats[k] .= match_cats[k] .& macro_cats[k]
    end
    filter_cats.met .&= megaclass.met_undiff

    filtered_in = Array{Int64}(undef, length(filter_cats.sed), length(classes))
    for i in eachindex(classes)
        for j in eachindex(filter_cats[classes[i]])
            filtered_in[j,i] = ifelse(filter_cats[classes[i]][j], 1, 0)
        end
    end


## --- Create a pseudo-filtered Tuple for P/Alkalinity calculations 
    # Metamorphic rocks are still only undifferentiated metamorphic samples 
    pseudo_cats = deepcopy(filter_cats)

    # Sedimentary rocks exclude undifferentiated sedimentary rocks 
    exclude_minor!(pseudo_cats)
    sed_undiff = copy(pseudo_cats.sed)
    include_minor!(pseudo_cats)
    pseudo_cats.sed .&= .!sed_undiff

    # We allow igneous rocks that were re-assigned to a rock type, because those samples 
    # weren't causing problems (and restricting samples means larger uncertainties)
    for k in (minorvolc..., minorplut..., minorign...,)
        pseudo_cats[k] .= match_cats[k]
    end


## --- Calculate resampling weights 
    k = invweight_age(sampleage)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)


## --- Resample data 
    # Remove any negative or 0 values from data 
    phosphorus[phosphorus .< 0] .= NaN;
    alkalinity[alkalinity .< 0] .= NaN;
    geochem_in[geochem_in .< 0] .= NaN;

    # Collect all class data and uncertainties
    class_all = hcat(matched_in, mapped_in, filtered_in)
    class_uncert = zeros(size(class_all));
    cats = get_cats(false, nsims)[2];
    cats_blank = (;
        match_cats = deepcopy(cats),
        macro_cats = deepcopy(cats),
        filter_cats = deepcopy(cats),
    )

    # Phosphorus / alkalinity ratio 
    target = (:sed, :carb, :shale, :volc, :plut)
    out = (:c,:m,:el,:eu)
    sim_ratio = NamedTuple{target}(
        NamedTuple{out}(Array{Float64}(undef, nbins) for _ in out) for _ in target
    )
    for k in target
        t = pseudo_cats[k]

        c,m,el,eu = bin_bsr_ratios(sampleage[t], vec(phosphorus)[t], vec(alkalinity)[t], 
            xmin, xmax, nbins,
            x_sigma = ageuncert[t],
            num_sigma = vec(phosphorus .* data_uncert)[t],
            denom_sigma = vec(alkalinity .* data_uncert)[t],
            p = p[t]
        )
        sim_ratio[k].c .= c 
        sim_ratio[k].m .= m 
        sim_ratio[k].eu .= eu 
        sim_ratio[k].el .= el
    end
    
    # Phosphorus and alkalinity data including cation mass balance 
    data_mol = hcat(Ca²⁺, Mg²⁺, K⁺, Na⁺, Fe²⁺, alkalinity, alkalinity_fe, phosphorus);
    resampled = bsresample(
        [data_mol sampleage class_all], 
        [data_mol .* data_uncert ageuncert class_uncert], 
        nsims, p
    )
    target = (:Ca, :Mg, :K, :Na, :Fe, :Alk, :Alk_Fe, :P, :Age)
    sim_mol = (;
        data = NamedTuple{target}(resampled[:,i] for i in eachindex(target)),
        cats = deepcopy(cats_blank),
    )
    i = length(target) + 1
    for k in keys(cats_blank)
        for c in classes
            sim_mol.cats[k][c] .= resampled[:,i] .> 0
            global i += 1
        end
    end
    
    # All element and rock class data 
    resampled = bsresample(
        [geochem_in sampleage class_all], 
        [geochem_in .* data_uncert ageuncert class_uncert], 
        nsims, p
    )
    target = (allelements..., :Age)
    sim_wt = (;
        data = NamedTuple{target}(resampled[:,i] for i in eachindex(target)),
        cats = deepcopy(cats_blank),
    )
    i = length(target) + 1
    for k in keys(cats_blank)
        for c in classes
            sim_wt.cats[k][c] .= resampled[:,i] .> 0
            global i += 1
        end
    end

    # Major rock classes inclusive of minors for all resampled data
    for k in keys(cats_blank)
        include_minor!(sim_mol.cats[k])
        include_minor!(sim_wt.cats[k])
    end


## --- Save resampled data to a file 
    fid = h5open("src/visualization/burial/resampled_geochem.h5", "w")
    g = create_group(fid, "vars")

    # Ratio data
    g_prime = create_group(g, "ratio")
        out = (:c,:m,:el,:eu)
        g_prime["ratio_head"] = string.(collect(out))
        for k in keys(sim_ratio)
            g_prime["$k"] = [sim_ratio[k].c sim_ratio[k].m sim_ratio[k].el sim_ratio[k].eu]
        end

    # Mol. data 
    g_prime = create_group(g, "mole")
    g_data = create_group(g_prime, "data")
        out = keys(sim_mol.data)
        g_data["data_head"] = string.(collect(out))
        a = Array{Float64}(undef, nsims, length(out))
        for i in eachindex(out)
            a[:,i] .= sim_mol.data[out[i]]
        end
        g_data["data"] = a

    g_cats = create_group(g_prime, "cats")
        g_cats["cats_head"] = string.(collect(keys(sim_mol.cats.match_cats)))
        for k in keys(sim_mol.cats)
            a = Array{Int64}(undef, nsims, length(sim_mol.cats[k]))
            for i in eachindex(keys(sim_mol.cats[k]))
                for j in eachindex(sim_mol.cats[k][i])
                    a[j,i] = ifelse(sim_mol.cats[k][i][j], 1, 0)
                end
            end
            g_cats["$k"] = a
        end

    # Wt.% data 
    g_prime = create_group(g, "wt")
    g_data = create_group(g_prime, "data")
        out = keys(sim_wt.data)
        g_data["data_head"] = string.(collect(out))
        a = Array{Float64}(undef, nsims, length(out))
        for i in eachindex(out)
            a[:,i] .= sim_wt.data[out[i]]
        end
        g_data["data"] = a
    g_cats = create_group(g_prime, "cats")
        g_cats["cats_head"] = string.(collect(keys(sim_wt.cats.match_cats)))
        for k in keys(sim_wt.cats)
            a = Array{Int64}(undef, nsims, length(sim_mol.cats[k]))
            for i in eachindex(keys(sim_wt.cats[k]))
                for j in eachindex(sim_wt.cats[k][i])
                    a[j,i] = ifelse(sim_wt.cats[k][i][j], 1, 0)
                end
            end
            g_cats["$k"] = a
        end

    close(fid)


## --- [PLOT] P/Alk ratios over time
    h = plot(
        xlabel="Age [Ma.]", ylabel="P / Alk [mol. ratio]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );

    target = keys(sim_ratio)
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(figs)
        hᵢ = deepcopy(h)
        plot!(hᵢ, sim_ratio[target[i]].c, sim_ratio[target[i]].m,
            yerror=(2*sim_ratio[target[i]].el, 2*sim_ratio[target[i]].eu),
            ribbon=(2*sim_ratio[target[i]].el, 2*sim_ratio[target[i]].eu),
            fillalpha=0.25,
            label="",
            title="$(target[i])",
            color=colors[target[i]], lcolor=colors[target[i]], msc=:auto,
            markershape=:circle,
            linewidth=2,
        )
        figs[i] = hᵢ
    end

    h = plot(figs..., layout=(2,3), size=(1800,800))
    display(h)
    savefig(h, "$filepath/p_alk_ratio.pdf")


## --- [PLOT] Phosphorus / alkalinity timeseries
    h = plot(
        xlabel="Age [Ma.]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );

    target = (:sed, :carb, :shale, :volc, :plut,)
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(target)
        f = sim_mol.cats.filter_cats[target[i]]

        hᵢ = deepcopy(h)
        c,m,e = binmeans(sim_mol.data.Age[f], sim_mol.data.Alk[f], xmin, xmax, nbins)
        plot!(hᵢ, c, m,
            yerror=2e, 
            label="", 
            title="$(target[i])",
            ylabel="Alkalinity [mol.]", 
            color=colors[target[i]], lcolor=colors[target[i]], msc=:auto,
                y_foreground_color_border=colors[target[i]],
                y_foreground_color_text=colors[target[i]],
                y_foreground_color_axis=colors[target[i]],
                y_guidefontcolor=colors[target[i]],
            markershape=:circle,
            linewidth=2,
        )

        c,m,e = binmeans(sim_mol.data.Age[f], sim_mol.data.P[f], xmin, xmax, nbins)
        plot!(twinx(), c, m,
            yerror=2e, 
            label="", 
            ylabel="Phosphorus [mol.]", 
            color=:black, lcolor=:black, msc=:auto,
            markershape=:circle,
            linewidth=2,
        )
        figs[i] = hᵢ
    end

    nplots = length(target)
    h = plot(figs..., layout=(2,3), size=(3*600,800))
    display(h)
    savefig(h, "$filepath/p_alk_abs.pdf")


## --- [PLOT] Fraction of matched samples mapped as the assigned class
    h = plot(
        xlabel="Age [Ma.]", ylabel="Fraction",
        # yaxis=:log10,
        # ylims=(0,1),
        xlims=(xmin, xmax),
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );

    target = (:sed, :carb, :shale, :volc, :plut, :met)
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(target)
        hᵢ = deepcopy(h)

        f = sim_wt.cats.filter_cats[target[i]]                  # Matched and mapped
        c,n₁ = bincounts(sim_wt.data.Age[f], xmin, xmax, nbins)
        f = sim_wt.cats.match_cats[target[i]]                   # Matched
        c,n₂ = bincounts(sim_wt.data.Age[f], xmin, xmax, nbins)

        n = n₁ ./ n₂ 
        plot!(hᵢ, c, n,
            label="", 
            title="$(target[i])",
            color=colors[target[i]], lcolor=colors[target[i]], msc=:auto,
            # markershape=:circle,
            # linewidth=2,
            seriestype=:bar,
            barwidths=((xmax-xmin)/nbins),
        )
        ylims!(0, min(1, nanmaximum(n)*1.1))
        figs[i] = hᵢ
    end

    nplots = length(target)
    h = plot(figs..., layout=(2,3), size=(3*600,800))
    display(h)
    savefig(h, "$filepath/data_density.pdf")


## --- End of file  