## --- Set up 
    # Potential correlation between the ratio of phosphorus to alkalininty and 
    # carbon burial 
    
    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots, StatsPlots, Colors

    # Save figures to: 
    filepath = "results/figures/burial"

    # Script-wide definitions 
    nsims = Int(1e6)
    xmin, xmax, nbins = 0, 3800, 38
    age_error = 0.05                   # Minimum age error (%)
    age_error_abs = 50                 # Minimum age error (Ma)

    # Conversion factors from wt.% element oxide to moles per 100g sample
    # Note that these are molar masses and must be **divided** from the wt.% [g/g] value
    CaO_to_Ca =   1 / (molarmass["Ca"]   + molarmass["O"]  )
    MgO_to_Mg =   1 / (molarmass["Mg"]   + molarmass["O"]  )
    K2O_to_K =    2 / (molarmass["K"] *2 + molarmass["O"]  )    # 2 mol K / 1 mol K₂O
    Na2O_to_Na =  2 / (molarmass["Na"]*2 + molarmass["O"]  )    # 2 mol Na / 1 mol Na₂O
    FeO_to_Fe =   1 / (molarmass["Fe"]   + molarmass["O"]  )
    P2O5_to_mol = 2 / (molarmass["P"] *2 + molarmass["O"]*5)    # 2 mol P / 1 mol P₂O₅


## --- Load data 
    # # Carbon isotope data
    # carbon = importdataset("data/carbonisotope/compilation.csv", ',', importas=:Tuple)

    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
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
    close(fid)


## --- Temporal resample P/Alk ratio 
    # L + ratio + we can't resample everything and then take a ratio, because that's bad 
    # for some reason. So calculate the ratio [actually, calculate X/(X+Y)], and then 
    # resample that before re-converting back into ratio. It's true! Google Keller and Schoene 
    # 2012 extended methods for more info

    # Definitions
    xmin, xmax, nbins = 0, 3800, 38
    uncert = 0.05                       # Percent error 

    # Get sample ages and uncertainties
    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
        uncert_rel=age_error, uncert_abs=age_error_abs
    )

    # Preallocate 
    alkalinity = Array{Float64}(undef, length(mbulk.SiO2), 1)
    phosphorus = Array{Float64}(undef, length(mbulk.SiO2), 1)

    # target = (:sed, :volc, :plut)
    target = (:plut, :volc, :sed)
    simout_ratio = NamedTuple{target}(Array{Float64}(undef, nbins, 4) for _ in target)
    simout_bulk = NamedTuple{target}(Array{Float64}(undef, nsims, 3) for _ in target)

    # Moles of alkalinity (charge), phosphorus, 
    for i in eachindex(alkalinity)
        Ca²⁺ = mbulk.CaO[i] * CaO_to_Ca * 2                 # +2 change
        Mg²⁺ = mbulk.MgO[i] * MgO_to_Mg * 2                 # +2 change
        K⁺ = mbulk.K2O[i] * K2O_to_K
        Na⁺ = mbulk.Na2O[i] * Na2O_to_Na
        alkalinity[i] = nansum([Ca²⁺, Mg²⁺, K⁺, Na⁺])

        phosphorus[i] = mbulk.P2O5[i] * P2O5_to_mol
    end
    alkalinity = vec(alkalinity)
    phosphorus = vec(phosphorus)

    # Calculate temporal weights and resample ratio
    for key in target 
        t = @. match_cats[key] & !isnan(phosphorus) & !isnan(alkalinity)
        t .&= .!(match_cats.phosphorite) .| (mbulk.P2O5 .< 4); 

        k = invweight_age(sampleage[t])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

        # Resample ratios 
        c,m,el,eu = bin_bsr_ratios(sampleage[t], phosphorus[t], alkalinity[t], 
            xmin, xmax, nbins,
            x_sigma = ageuncert[t],
            num_sigma = phosphorus[t]*uncert,
            denom_sigma = alkalinity[t]*uncert,
            p = p
        )
        simout_ratio[key] .= [c m el eu]

        # Resample values
        simout_bulk[key] .= bsresample([sampleage[t] phosphorus[t] alkalinity[t]], 
            [ageuncert[t] phosphorus[t]*uncert alkalinity[t]*uncert], 
            nsims, p
        )
    end

    # Sanity check plot 
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(simout_ratio))
    p = palette([:red, :hotpink, :seagreen], 3)
    c,m,el,eu = 1,2,3,4;
    h = plot(
        # xlabel="Age [Ma.]", ylabel="P / Alk [mol. ratio]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
    );
    for i in eachindex(target)
        key = target[i]
        h1 = deepcopy(h)
        plot!(h1, simout_ratio[key][:,c], simout_ratio[key][:,m], 
            yerror=(2*simout_ratio[key][:,el], 2*simout_ratio[key][:,eu]), 
            label="$key", 
            color=p[i], lcolor=p[i], msc=:auto, 
            markershape=:circle,
        )
        figs[i] = h1
    end
    
    h2 = plot(figs..., layout=(1,3), size=(1800,400))
    display(h2)
    savefig(h2, "$filepath/p_alk.pdf")


## --- Plot resampled alkalinity / phosphorus over time 
    # Sedimentary samples (likely problems)
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        title="Resampled P and Alk"
    );
    c,m,e = binmeans(simout_bulk.sed[:,1], simout_bulk.sed[:,2], xmin, xmax, nbins)
    plot!(c,m,yerror=2e, label="", ylabel="P [mol.]", 
        color=:red, lcolor=:red, msc=:auto,
        y_foreground_color_border=:red,
        y_foreground_color_text=:red,
        y_foreground_color_axis=:red,
        y_guidefontcolor=:red,
    )
    
    c,m,e = binmeans(simout_bulk.sed[:,1], simout_bulk.sed[:,3], xmin, xmax, nbins)
    plot!(twinx(), c,m,yerror=2e, label="", ylabel="Alk [mol.]", 
        color=:blue, lcolor=:blue, msc=:auto,
        y_foreground_color_border=:blue,
        y_foreground_color_text=:blue,
        y_foreground_color_axis=:blue,
        y_guidefontcolor=:blue,
    )


## --- ??? What if we calculate the ratio after resampling 
    # ratio = simout_bulk.sed[:,2] ./ simout_bulk.sed[:,3]
    # c,m,e = binmeans(simout_bulk.sed[:,1], ratio, xmin, xmax, nbins)
    # h = plot(c,m,yerror=2e, label="", ylabel="P / Alk [mol.]", 
    #     color=:blue, lcolor=:blue, msc=:auto,
    # )
    # display(h)

    # # That's fucked; try a new way
    # c,m₁,e₁ = binmeans(simout_bulk.sed[:,1], simout_bulk.sed[:,2], xmin, xmax, nbins)
    # c,m₂,e₂ = binmeans(simout_bulk.sed[:,1], simout_bulk.sed[:,3], xmin, xmax, nbins)
    # r = (m₁ .± e₁) ./ (m₂ .± e₂)
    # m = Measurements.value.(r)
    # e = Measurements.uncertainty.(r)
    # h = plot(c,m,yerror=2e, label="", ylabel="P / Alk [mol.]", 
    #     color=:blue, lcolor=:blue, msc=:auto,
    # )
    # display(h)


## --- Resample moles of each alkalinity cation in seds, preserving class data 
    # Definitions 
    uncert = 0.05                       # Percent error 
    xmin, xmax = 0, 3800
    nbins = Int((xmax-xmin)/100)

    # Get sample ages and uncertainties
    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
        uncert_rel=age_error, uncert_abs=age_error_abs
    )

    # Get rock classes
    include_minor!(match_cats)
    seds = (get_rock_class()[2]..., :sed)
    a = Array{Int64}(undef, length(match_cats.sed), length(seds))
    for i in eachindex(seds)
        for j in eachindex(match_cats[seds[i]])
            a[j,i] = ifelse(match_cats[seds[i]][j], 1, 0)
        end
    end

    # Get cation mole abundance. Calculate Fe²⁺, but don't include it in calculations 
    # since we're concerned with the all geologic time. Potential to include this for 
    # Archean sediments in a later iteration....
    Ca²⁺ = @. mbulk.CaO * CaO_to_Ca     # x2 for charge!!
    Mg²⁺ = @. mbulk.MgO * MgO_to_Mg     # x2 for charge!!
    K⁺ = @. mbulk.K2O * K2O_to_K
    Na⁺ = @. mbulk.Na2O * Na2O_to_Na
    Fe²⁺ = @. mbulk.FeOT * FeO_to_Fe    # x2 for charge!!
    alkalinity = nansum([2*Ca²⁺ 2*Mg²⁺ K⁺ Na⁺], dims=2)
    # alkalinity = nansum([2*Ca²⁺ 2*Mg²⁺ K⁺ Na⁺ 2*Fe²⁺], dims=2)

    # Calculate resampling weights
    k = invweight_age(sampleage[match_cats.sed])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    # Restrict data and remove any zero values                                                       
    t = match_cats.sed
    data = [Ca²⁺[t] Mg²⁺[t] K⁺[t] Na⁺[t] Fe²⁺[t] alkalinity[t]]
    nanzero!(data)
    resampled = bsresample(
        [data sampleage[t] a[t,:]], 
        [data.*uncert ageuncert[t] zeros(size(a[t,:]))], 
        nsims, p
    )

    # Indexing lookup to avoid indexing errors
    cations_ferrous = (:Ca, :Mg, :K, :Na, :Fe,)
    cations = (:Ca, :Mg, :K, :Na,)
    target = (cations_ferrous..., :Alk, :Age, seds...,)
    r_index = NamedTuple{target}(1:length(target))

    # Create a new resampled array for charges to avoid confusion between cation 
    # abundance and alkalinity 
    resampled_charge = copy(resampled)
    resampled_charge[:,r_index.Ca] .*= 2;
    resampled_charge[:,r_index.Mg] .*= 2;
    resampled_charge[:,r_index.Fe] .*= 2;

    # Re-parse rock class data
    sim_cats = NamedTuple{seds}(resampled[:,r_index[k]] .> 0 for k in seds)
    for k in seds
        sim_cats.sed .|= sim_cats[k]
    end

    # Get age data since we'll be accessing it a lot 
    sim_age = resampled[:,r_index.Age]
    arc = @. 2500 < sim_age < 3800;       # Filter for Archean samples 


## --- [PLOT] Agreement in alkalinity vs. sum of charges in resampled data
    sim_alk = Array{Float64}(undef, size(resampled)[1], length(cations))
    for i in eachindex(cations)
        sim_alk[:,i] .= resampled_charge[:,r_index[cations[i]]]
    end
    sim_alk = nansum(sim_alk, dims=2)
    alk_difference = sim_alk .- resampled[:,r_index.Alk];
    h = histogram(alk_difference, label="",
        xlabel="Difference [mol.] -- Negative if resampled larger than calculated", 
        ylabel="Abundance",
        framestyle=:box,
        fontfamily=:Helvetica,
        color=:black,
        title="Resampled Alkalinity vs. Sum of Charges", titleloc=:left,
    )
    display(h)


## --- [PLOT] Bulk sedimentary cation abundance and alkalinity
    p = palette(:rainbow, 5)
    labels = ["Ca²⁺", "Mg²⁺", "K⁺", "Na⁺", "Fe²⁺"]
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        ylabel="Cation [mol.]",
        legend=:top
    ); 
    for i in eachindex(labels)
        c,m,e = binmeans(sim_age[arc], resampled[:,i][arc], 2500, 3800, 13)
        plot!(h, c,m,yerror=2e, 
            color=p[i], lcolor=p[i], msc=:auto, 
            markershape=:circle, label=labels[i]
        )
    end
    c,m,e = binmeans(sim_age[arc], resampled[:,r_index.Alk][arc], 2500, 3800, 13)
    plot!(twinx(), c,m,yerror=2e, 
        label="Alkalinity", legend=:topright, fg_color_legend=:white,
        ylabel="Alk [mol.]", color=:black,
        title="Resampled Alkalinity, Cation Abundance", titleloc=:left,
    )
    display(h)


## --- [PLOT] Alkalinity by rock class of interest
    p = palette(:berlin, (length(seds)))
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        xlabel="Age [Ma.]", ylabel="Alkalinity [mol.]",
        legend=:outerright, size=(800,600),
        title="Alkalinity by Rock Class", titleloc=:left
    )
    for i in eachindex(seds)
        s = arc .& sim_cats[seds[i]]
        count(s) == 0 && continue
        c,m,e = binmeans(sim_age[s], resampled[:,r_index.Alk][s], 2500, 3800, 13)
        s = .!isnan.(m)
        plot!(h, c[s], m[s], yerror=2e, label="$(seds[i])",
            linewidth=2, markershape=:circle,
            color=p[i], lcolor=p[i], msc=:auto
        )
    end
    display(h)


## --- [PLOT] Do low values correlate with missing data? 
    p = palette(:rainbow, length(cations))
    labels = ["Ca²⁺", "Mg²⁺", "K⁺", "Na⁺", "Fe²⁺"]

    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlabel="Age [Ma.]", ylabel="Sample Count",
        legend=:top, fg_color_legend=:white,
        y_foreground_color_border=:red,
        y_foreground_color_text=:red,
        y_foreground_color_axis=:red,
        y_guidefontcolor=:red,
    )
    for i in eachindex(cations)
        s = arc .& .!isnan.(resampled[:,i])
        c,n = bincounts(sim_age[s], 2500, 3800, 13)
        plot!(c, n, markershape=:circle, label="$(labels[i])",
            color=p[i], lcolor=p[i], msc=:auto, 
        )
    end
    c,m,e = binmeans(sim_age[arc], resampled[:,r_index.Alk][arc], 2500, 3800, 13)
    plot!(twinx(), c,m,yerror=2e, label="", ylabel="Alk [mol.]", color=:black,)
    display(h)


## --- [PLOT] Are low values caused by a spike in chert abundance? 
    # Still restricting to the Archean...
    c,n₁ = bincounts(resampled[:,r_index.Age][arc], 2500,3800,13)
    c,n₂ = bincounts(resampled[:,r_index.Age][arc .& sim_cats.chert], 2500,3800,13)
    n = float.(n₂) ./ float.(n₁) .* 100
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlabel="Age [Ma.]", ylabel="Fractional Abundance",
        fg_color_legend=:white, legend=:top,
    )
    c,m,e = binmeans(sim_age[arc], resampled[:,r_index.Alk][arc], 2500,3800,13)
    plot!(c,m,yerror=2e, label="All Seds", ylabel="Alk [mol.]", 
        color=:seagreen, lcolor=:seagreen, msc=:auto
    )
    c,m,e = binmeans(sim_age[arc .& sim_cats.chert], resampled[:,r_index.Alk][arc .& sim_cats.chert], 2500,3800,13)
    plot!(c,m,yerror=2e, label="Chert", ylabel="Alk [mol.]", 
        color=:black, lcolor=:black, msc=:auto
    )
    plot!(twinx(), c, n, label="",
        ylabel="Relative Chert Abundance [%]",  # % of rocks that are chert 
        color=colors.chert,
        linewidth=2,
        markershape=:circle, msc=:auto,
        y_foreground_color_border=colors.chert,
        y_foreground_color_text=colors.chert,
        y_foreground_color_axis=colors.chert,
        y_guidefontcolor=colors.chert,
    )
    display(h)

    # I wonder if there's something about age uncertainty that could explain the 
    # decrease in alkalinity... the big spike in chert abundance at ~2850 does 
    # line up with a decrease in the alkalinity recorded in chert, and a bit with a 
    # decrease in alkalinity for all seds... but it doesn't line up with the big drop. 
    # Could that just be because we don't have good ages for these rocks? Do note that 
    # we see drops in alkalinity in most rock types though...


## --- [PLOT] Marine sediment alkalinity over time?
    p = (;
        carb=:teal,
        shale=:darkorange,
    )
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlabel="Age [Ma.]", 
        fg_color_legend=:white,
        title="Marine Alkalinity [GOE Marked]", titleloc=:left,
    )
    c,m,e = binmeans(sim_age[sim_cats.shale], resampled[:,r_index.Alk][sim_cats.shale], 
        xmin, xmax, nbins)
    plot!(c,m,yerror=2e, label="", markershape=:circle, 
        color=p.shale, lcolor=p.shale, msc=:auto,
        legend=:topleft, fg_color_legend=:white,
        ylabel="Shale Alkalinity [mol.]",
        y_foreground_color_border=p.shale,
        y_foreground_color_text=p.shale,
        y_foreground_color_axis=p.shale,
        y_guidefontcolor=p.shale,
    )
    c,m,e = binmeans(sim_age[sim_cats.carb], resampled[:,r_index.Alk][sim_cats.carb], 
        xmin, xmax, nbins)
    plot!(twinx(), c,m,yerror=2e, label="", markershape=:circle, 
        color=p.carb, lcolor=p.carb, msc=:auto,
        legend=:bottomleft, fg_color_legend=:white,
        ylabel="Carbonate Alkalinity [mol.]",
        y_foreground_color_border=p.carb,
        y_foreground_color_text=p.carb,
        y_foreground_color_axis=p.carb,
        y_guidefontcolor=p.carb,
    )
    vline!([2500], label="", color=:black, linestyle=:dash)
    display(h)


## --- End of file  