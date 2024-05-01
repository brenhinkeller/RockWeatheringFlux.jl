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
    P2O5_err = 2.0      # Error set as 10x the error given in 
    alk_err = 1.0       # volcanic.mat (Keller et al., 2015)   

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
    CaO_to_Ca =   (molarmass["Ca"]   + molarmass["O"]  )
    MgO_to_Mg =   (molarmass["Mg"]   + molarmass["O"]  )
    K2O_to_K =    (molarmass["K"] *2 + molarmass["O"]  ) * 2   # 2 mol K / 1 mol K₂O
    Na2O_to_Na =  (molarmass["Na"]*2 + molarmass["O"]  ) * 2   # 2 mol Na / 1 mol Na₂O
    P2O5_to_mol = (molarmass["P"] *2 + molarmass["O"]*5)

    # Moles of alkalinity (charge), phosphorus, 
    for i in eachindex(alkalinity)
        Ca²⁺ = mbulk.CaO[i] * CaO_to_Ca * 2     # +2 change
        Mg²⁺ = mbulk.MgO[i] * MgO_to_Mg * 2     # +2 change
        K⁺ = mbulk.K2O[i] * K2O_to_K
        Na⁺ = mbulk.Na2O[i] * Na2O_to_Na
        alkalinity[i] = Ca²⁺ + Mg²⁺ + K⁺ + Na⁺

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
            num_sigma = fill(P2O5_err, count(t)),
            denom_sigma = fill(alk_err, count(t)),
            p = p
        )
        simout_ratio[key] .= [c m el eu]

        # Resample values
        simout_bulk[key] .= bsresample([sampleage[t] phosphorus[t] alkalinity[t]], 
            [ageuncert[t] fill(P2O5_err, count(t)) fill(alk_err, count(t))], 
            nsims, p
        )
    end

    # Sanity check plot 
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(simout))
    p = palette([:red, :hotpink, :seagreen], 3)
    c,m,el,eu = 1,2,3,4;
    h = plot(
        # xlabel="Age [Ma.]", ylabel="P / Alk [mol. ratio]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topleft,
    );
    for i in eachindex(target)
        key = target[i]
        h1 = deepcopy(h)
        plot!(h1, simout[key][:,c], simout[key][:,m], 
            yerror=(2*simout[key][:,el], 2*simout[key][:,eu]), 
            label="$key", 
            color=p[i], lcolor=p[i], msc=:auto, 
            markershape=:circle,
        )
        figs[i] = h1
    end
    
    h2 = plot(figs..., layout=(1,3), size=(1800,400))
    display(h2)
    # savefig(h2, "$filepath/p_alk.pdf")


## --- Plot resampled alkalinity / phosphorus over time 
    # Sedimentary samples (likely problems)
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
    )
    c,m,e = binmeans(simout_bulk.sed[:,1], simout_bulk.sed[:,2], xmin, xmax, nbins)
    plot!(c,m,yerror=2e, label="", ylabel="P [mol.]", color=:red)
    
    c,m,e = binmeans(simout_bulk.sed[:,1], simout_bulk.sed[:,3], xmin, xmax, nbins)
    plot!(twinx(), c,m,yerror=2e, label="", ylabel="Alk [mol.]", color=:blue)


## --- Plot moles of each alkalinity cation over the Archean
    Ca²⁺ = Array{Float64}(undef, length(mbulk.CaO))
    Mg²⁺ = Array{Float64}(undef, length(mbulk.CaO))
    K⁺ = Array{Float64}(undef, length(mbulk.CaO))
    Na⁺ = Array{Float64}(undef, length(mbulk.CaO))
    for i in eachindex(alkalinity)
        Ca²⁺[i] = mbulk.CaO[i] * CaO_to_Ca
        Mg²⁺[i] = mbulk.MgO[i] * MgO_to_Mg
        K⁺[i] = mbulk.K2O[i] * K2O_to_K
        Na⁺[i] = mbulk.Na2O[i] * Na2O_to_Na
        # alkalinity[i] = Ca²⁺ + Mg²⁺ + K⁺ + Na⁺
    end

    t = match_cats.ign
    # t = trues(length(match_cats.sed))
    p = palette(:rainbow, 4)
    h = plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        ylabel="Cation [mol.]",
        legend=:top
    );
    c,m,e = binmeans(sampleage[t], Ca²⁺[t], xmin, xmax, nbins)
    plot!(c,m,yerror=2e, color=p[1], lcolor=p[1], msc=:auto, markershape=:circle, label="Ca²⁺")
    c,m,e = binmeans(sampleage[t], Mg²⁺[t], xmin, xmax, nbins)
    plot!(c,m,yerror=2e, color=p[2], lcolor=p[2], msc=:auto, markershape=:circle, label="Mg²⁺")
    c,m,e = binmeans(sampleage[t], K⁺[t], xmin, xmax, nbins)
    plot!(c,m,yerror=2e, color=p[3], lcolor=p[3], msc=:auto, markershape=:circle, label="K⁺")
    c,m,e = binmeans(sampleage[t], Na⁺[t], xmin, xmax, nbins)
    plot!(c,m,yerror=2e, color=p[4], lcolor=p[4], msc=:auto, markershape=:circle, label="Na⁺")

    c,m,e = binmeans(sampleage[t], alkalinity[t], xmin, xmax, nbins)
    plot!(twinx(), c,m,yerror=2e, label="", ylabel="Alk [mol.]", color=:black)

    xlims!(2500, 3900)
    display(h)

## --- End of file 