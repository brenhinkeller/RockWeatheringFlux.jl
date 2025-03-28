## --- Set up 
    # Parse GloRiSe (suspended sediment) and GLORICH (dissolved river chemistry) data 

    # Packages 
    using RockWeatheringFlux

    # Load data 
    ID = importdataset("data/GloRiSe/SedimentDatabase_ID.csv", ',', importas=:Tuple);
    TE = importdataset("data/GloRiSe/SedimentDatabase_TE.csv", ',', importas=:Tuple);
    ME = importdataset("data/GloRiSe/SedimentDatabase_ME_Nut.csv", ',', importas=:Tuple);
    hydro = importdataset("data/GLORICH/hydrochemistry.csv", ',', importas=:Tuple);

    # Required elements
    majors, minors = get_elements();
    allelements = Tuple([majors; minors])
        

## --- [GloRiSE] Remove outliers 
    # Outlier removal can't deal with any missing keys, and it needs all data to be in 
    # wt.%. We can't combine the GloRiSe data into a single dataset because not all 
    # samples in the ME and TE datasets have sample IDs, and the datasets are different 
    # lengths. Therefore, make datasets with the required elements, but have most of them
    # blank. We can ignore any location / sample ID data because we can't use it anyway

    # Outlier screening differentiates by rock type, just set everything to true here 
    ME_cats = get_cats(false, length(ME[1]))[2];
    TE_cats = get_cats(false, length(TE[1]))[2];
    [ME_cats[k] .= true for k in keys(ME_cats)];
    [TE_cats[k] .= true for k in keys(TE_cats)];

    # Major elements
    out_ME = NamedTuple{allelements}(Array{Float64}(undef, length(ME[1])) for _ in allelements);
    for e in allelements 
        if haskey(ME, Symbol(string(e)*"_wt"))
            out_ME[e] .= ME[Symbol(string(e)*"_wt")]
        elseif e==:P2O5
        # P2O5 
            out_ME.P2O5 .= ME.P2O5_tot_wt
        elseif e==:FeOT
        # Fe2O3 to FeO
            out_ME.FeOT .= feoconversion.(NaN, NaN, NaN, ME.Fe2O3T_wt)
        else 
        # Otherwise, no data and set to NaN 
            out_ME[e] .= NaN
        end
    end
    screen_outliers!(out_ME, ME_cats)

    # Trace elements 
    out_TE = NamedTuple{allelements}(Array{Float64}(undef, length(TE[1])) for _ in allelements);
    for e in allelements 
        if e==:I
        # This data is messed up (no data except for a "Zr" string)
            out_TE[e] .= NaN
        elseif haskey(TE, Symbol(string(e)*"_ppm"))
            out_TE[e] .= TE[Symbol(string(e)*"_ppm")] ./ 10_000
        else
            out_TE[e] .= NaN 
        end
    end
    screen_outliers!(out_TE, TE_cats)

    # How far off is bulk averaging major elements from the Muller et al., estimate?
    # Make sure to anhydrous normalize our things, since they're anhydrous
    majoravg = normalize!([nanmean(out_ME[k]) for k in majors])


## --- [GloRiSE] Save to a file 
    # Majors 
    anhydrous_majors = majors[1:end .!= findfirst(x->x==:Volatiles, majors)]
    result = Array{Float64}(undef, length(out_ME[1]), length(anhydrous_majors))
    rows = reshape(string.(anhydrous_majors), 1, :)
    for i in eachindex(anhydrous_majors)
        result[:,i] .= out_ME[anhydrous_majors[i]]
    end
    writedlm("output/GlobalRivers/GloRiSe_major_screened.tsv", vcat(rows, result))

    # Minors 
    result = Array{Float64}(undef, length(out_TE[1]), length(minors))
    rows = reshape(string.(minors), 1, :)
    for i in eachindex(minors)
        result[:,i] .= out_TE[minors[i]]
    end
    writedlm("output/GlobalRivers/GloRiSe_minor_screened.tsv", vcat(rows, result))


## --- [GLORICH]
    # Preallocate 
    nhydro = length(hydro[1])
    hydro_out = NamedTuple{allelements}(Array{Float64}(undef, nhydro) for _ in allelements);  

    # The river chemistry is... tricky. I really want to just look at the major elements
    # (REEs are generally insoluable so like... not reported) but I also do want to screen 
    # for outliers. So What I'm gonna do is as above, but I'm gonna make a NamedTuple that's 
    # mostly blank to send it through screen_outliers!(). Then I'll just take the stuff I 
    # want and save it to a file.

    # Convert umol / L to wt.% element oxide
        # umol -> mol (1e-6); g/kg -> wt.% (1e-1), 
        # mol -> g (molar mass Ca), Ca -> CaO (molar mass CaO / molar mass Ca)
        # divide Na2O and K2O by 2 because 2 cations / oxide
    hydro_out.CaO .= hydro.Ca * 1e-7 * molarmasspercation["CaO"]
    hydro_out.MgO .= hydro.Mg * 1e-7 * molarmasspercation["MgO"]
    hydro_out.Na2O .= hydro.Na * 1e-7 * molarmasspercation["Na2O"] / 2
    hydro_out.K2O .= hydro.K * 1e-7 * molarmasspercation["K2O"] / 2
    hydro_out.SiO2 .= hydro.SiO2 * 1e-7 * molarmasspercation["SiO2"]
    # No TiO2, Al2O3, Fe2O3 or FeO (these are insoluable so like... of course)

    # Screen outliers as described above
    hydro_cats = get_cats(false, nhydro)[2];
    [hydro_cats[k] .= true for k in keys(hydro_cats)];
    screen_outliers!(hydro_out, hydro_cats)

    # What's the total concentration of all dissolved load? 
    otherload = nansum(hcat(
        # Add only non-carbon elements from carbonate ions so we don't double-count carbon
        hydro.HCO3 * 1e-7 * (molarmass["H"] + molarmass["O"]*3),
        hydro.CO3 * 1e-7 * (molarmass["O"]*3),

        # Everything else as normal 
        hydro.Cl * 1e-7 * (molarmass["Cl"]),
        hydro.F * 1e-7 * (molarmass["F"]),
        hydro.DSr * 1e-7 * (molarmass["Sr"]),   # Dissolved Sr 
        hydro.TC * 1e-7 * (molarmass["C"]),     # Total carbon 
        hydro.TN * 1e-7 * (molarmass["N"]),     # Total nitrogen 
        hydro.TP * 1e-7 * (molarmass["P"]),     # Total phosphorus 
    ), dims=2)
    totalload = nansum(hcat(otherload, hydro_out.CaO, hydro_out.MgO, hydro_out.Na2O, 
        hydro_out.K2O, hydro_out.SiO2), dims=2
    )

    # Normalize to 100% 

    # Save just a file of major elements (just so things don't totally go off the rails) 
    out = vcat(
        ["SiO2" "Al2O3" "FeOT" "TiO2" "MgO" "CaO" "Na2O" "K2O" "Volatiles"],
        hcat(
            hydro_out.SiO2,     
            fill(NaN, nhydro),  # Al2O3
            fill(NaN, nhydro),  # FeOT
            fill(NaN, nhydro),  # TiO2 
            hydro_out.MgO,
            hydro_out.CaO,
            hydro_out.Na2O,
            hydro_out.K2O,
            fill(NaN, nhydro),  # Volatiles
        )
    )
    writedlm("output/GlobalRivers/GLORICH_major_screened.tsv", out)

    
## --- End of File 