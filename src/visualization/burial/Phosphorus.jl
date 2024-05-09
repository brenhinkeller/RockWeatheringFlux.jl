## --- Set up 
    # Potential correlation between the ratio of phosphorus to alkalininty and 
    # carbon burial 
    
    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots, Colors

    # Save figures to: 
    filepath = "results/figures/burial"

    # Script-wide definitions 
    nsims = Int(1e6)
    xmin, xmax, nbins = 0, 3800, 38
    age_error = 0.05                   # Minimum age error (%)
    age_error_abs = 50                 # Minimum age error (Ma)


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

    # Conversion factors from wt.% element oxide to moles per 100g sample
    # Note that these are molar masses and must be **divided** from the wt.% [g/g] value
    CaO_to_Ca =   1 / (molarmass["Ca"]   + molarmass["O"]  )
    MgO_to_Mg =   1 / (molarmass["Mg"]   + molarmass["O"]  )
    K2O_to_K =    2 / (molarmass["K"] *2 + molarmass["O"]  )    # 2 mol K / 1 mol K₂O
    Na2O_to_Na =  2 / (molarmass["Na"]*2 + molarmass["O"]  )    # 2 mol Na / 1 mol Na₂O
    FeO_to_Fe =   1 / (molarmass["Fe"]   + molarmass["O"]  )
    P2O5_to_mol = 2 / (molarmass["P"] *2 + molarmass["O"]*5)    # 2 mol P / 1 mol P₂O₅

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


## --- Plot moles of each alkalinity cation in Archean seds     
    # Get individual cation abundances and resample.
    Ca²⁺ = @. mbulk.CaO * CaO_to_Ca     # x2 for charge!!
    Mg²⁺ = @. mbulk.MgO * MgO_to_Mg     # x2 for charge!!
    K⁺ = @. mbulk.K2O * K2O_to_K
    Na⁺ = @. mbulk.Na2O * Na2O_to_Na
    Fe²⁺ = @. mbulk.FeOT * FeO_to_Fe    # x2 for charge!!
    # alkalinity2 = nansum([2*Ca²⁺ 2*Mg²⁺ K⁺ Na⁺], dims=2)
    alkalinity2 = nansum([2*Ca²⁺ 2*Mg²⁺ K⁺ Na⁺ 2*Fe²⁺], dims=2)

    k = invweight_age(sampleage[match_cats.sed])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    t = @. match_cats.sed & (2500 < sampleage < 3800)
    # t .&= .!isnan.(Ca²⁺) .& .!isnan.(Mg²⁺) .& .!isnan.(K⁺) .& .!isnan.(Na⁺)

    data = [Ca²⁺[t] Mg²⁺[t] K⁺[t] Na⁺[t] Fe²⁺[t] alkalinity2[t]]
    nanzero!(data)
    uncert = ones(size(data)) ./ 100
    resampled = bsresample([data sampleage[t]], [uncert ageuncert[t]], nsims, p)


## --- Plot data 
    xmin, xmax = 2500, 3800
    nbins = Int((xmax-xmin)/100)

    p = palette(:rainbow, 5)
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        ylabel="Cation [mol.]",
        legend=:top
    );

    # # Real data 
    # h1 = deepcopy(h)
    # t = @. match_cats.sed .& (xmin .< sampleage .< xmax)
    # # t = trues(length(match_cats.sed))
    # c,m,e = binmeans(sampleage[t], Ca²⁺[t], xmin, xmax, nbins)
    # plot!(h1, c,m,yerror=2e, color=p[1], lcolor=p[1], msc=:auto, markershape=:circle, label="Ca²⁺")
    # c,m,e = binmeans(sampleage[t], Mg²⁺[t], xmin, xmax, nbins)
    # plot!(h1, c,m,yerror=2e, color=p[2], lcolor=p[2], msc=:auto, markershape=:circle, label="Mg²⁺")
    # c,m,e = binmeans(sampleage[t], K⁺[t], xmin, xmax, nbins)
    # plot!(h1, c,m,yerror=2e, color=p[3], lcolor=p[3], msc=:auto, markershape=:circle, label="K⁺")
    # c,m,e = binmeans(sampleage[t], Na⁺[t], xmin, xmax, nbins)
    # plot!(h1, c,m,yerror=2e, color=p[4], lcolor=p[4], msc=:auto, markershape=:circle, label="Na⁺")
    # c,m,e = binmeans(sampleage[t], Na⁺[t], xmin, xmax, nbins)
    # plot!(h1, c,m,yerror=2e, color=p[4], lcolor=p[4], msc=:auto, markershape=:circle, label="Na⁺")
    # c,m,e = binmeans(sampleage[t], Fe²⁺[t], xmin, xmax, nbins)
    # plot!(h1, c,m,yerror=2e, color=p[5], lcolor=p[5], msc=:auto, markershape=:circle, label="Fe²⁺")
    # c,m,e = binmeans(sampleage[t], alkalinity2[t], xmin, xmax, nbins)
    # plot!(twinx(), c,m,yerror=2e, label="", ylabel="Alk [mol.]", color=:black,
    #     title="Observed", titleloc=:left
    # )
    # # xlims!(2500,3800)
    # display(h1)

    # Shouldn't adding Fe change the alkalinity??
    # t = @. match_cats.sed .& (xmin .< sampleage .< xmax)
    # c,m,e = binmeans(sampleage[t], alkalinity2[t], xmin, xmax, nbins)
    # h1 = deepcopy(h)
    # plot!(h1,c,m,yerror=2e, label="", ylabel="Alk [mol.]", color=:black,)

    # Resampled data 
    h1 = deepcopy(h)
    age = resampled[:,end]
    t = @. 2500 < age < 3800;
    labels = ["Ca²⁺", "Mg²⁺", "K⁺", "Na⁺", "Fe²⁺"]
    for i in eachindex(labels)
        c,m,e = binmeans(age[t], resampled[:,i][t], xmin, xmax, nbins)
        plot!(h1, c,m,yerror=2e, 
            color=p[i], lcolor=p[i], msc=:auto, 
            markershape=:circle, label=labels[i]
        )
    end

    c,m,e = binmeans(age[t], resampled[:,end-1][t], xmin, xmax, nbins)
    plot!(twinx(), c,m,yerror=2e, 
        label="Alkalinity", legend=:topright, fg_color_legend=:white,
        ylabel="Alk [mol.]", color=:black,
        title="Resampled", titleloc=:left,
    )
    xlims!(2500, 3900)
    display(h1)


## --- Missing data or low values?
    # Concerned with K, Na, Mg, and Ca in that order 
    # Look at resampled data 
    xmin, xmax = 2500, 3800
    nbins = Int((xmax-xmin)/100)

    age = resampled[:,end]
    t = @. 2500 < age < 3800;
    labels = ["Ca²⁺", "Mg²⁺", "K⁺", "Na⁺", "Fe²+"]

    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlabel="Age [Ma.]", ylabel="Sample Count",
        legend=:top
    )
    for i in eachindex(labels)
        s = t .& .!isnan.(resampled[:,i])
        c,n = bincounts(age[s], xmin, xmax, nbins)
        plot!(c, n, markershape=:circle, label="$(labels[i])",
            color=p[i], lcolor=p[i], msc=:auto, 
        )
    end
    c,m,e = binmeans(age[t], resampled[:,end-1][t], xmin, xmax, nbins)
    plot!(twinx(), c,m,yerror=2e, label="", ylabel="Alk [mol.]", color=:black,
    )


## --- Check chert abundance over time 
    # Definitions
    xmin, xmax, nbins = 0, 3800, 38
    nsims = Int(1e6)
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];

    # Try geochemical age, but it may make more sense to use map age
    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
        uncert_rel=age_error, uncert_abs=age_error_abs
    )

    # Resample sedimentary rock types over geologic time 
    a = Array{Int64}(undef, count(match_cats.sed), length(minorsed))
    for i in eachindex(minorsed)
        target = match_cats[minorsed[i]][match_cats.sed]
        for j in eachindex(target)
            a[j,i] = ifelse(target[j], 1, 0)
        end
    end

    k = invweight_age(sampleage[match_cats.sed])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    resampled = bsresample([sampleage[match_cats.sed] a], 
        [ageuncert[match_cats.sed] zeros(size(a))], 
        nsims, p
    )

    # Count abundance of cherts in each bin over time 
    sedindex = NamedTuple{minorsed}(collect(1:length(minorsed)) .+ 1)
    t = resampled[:,sedindex.chert] .== 1;                   # Find cherts 
    # s = resampled[]

    c,n₁ = bincounts(resampled[:,1], xmin, xmax, nbins)      # Absolute counts of all seds
    c,n₂ = bincounts(resampled[:,1][t], xmin, xmax, nbins)   # Absolute counts of cherts 
    n = float.(n₂) ./ float.(n₁)                             # Relative abundance

    plot(c, n, label="",
        color=colors.chert,
        markershape=:circle,
        framestyle=:box,
        fontfamily=:Helvetica,
        xlabel="Age [Ma.]", ylabel="Fractional Abundance"    
    )
    # plot!(twinx())


    # cherts = Array{Float64}(undef, length(c))
    # for i in eachindex(cherts)
    #     t = 
    # end


## --- End of file  