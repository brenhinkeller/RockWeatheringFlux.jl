## --- Set up 
    # Parse GloRiSe database data 

    # Packages 
    using RockWeatheringFlux

    # Load data 
    ID = importdataset("data/GloRiSe/SedimentDatabase_ID.csv", ',', importas=:Tuple);
    TE = importdataset("data/GloRiSe/SedimentDatabase_TE.csv", ',', importas=:Tuple);
    ME = importdataset("data/GloRiSe/SedimentDatabase_ME_Nut.csv", ',', importas=:Tuple);

    
## --- Remove outliers 
    # Outlier removal can't deal with any missing keys, and it needs all data to be in 
    # wt.%. We can't combine the GloRiSe data into a single dataset because not all 
    # samples in the ME and TE datasets have sample IDs, and the datasets are different 
    # lengths. Therefore, make datasets with the required elements, but have most of them
    # blank. We can ignore any location / sample ID data because we can't use it anyway

    # Required elements
    majors, minors = get_elements();
    allelements = Tuple([majors; minors])
    
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


## --- How far off is bulk averaging major elements from the Muller et al., estimate?
    # Answer: very :/
    majoravg = [nanmean(out_ME[k]) for k in majors]

## --- Save to a file 
    # Majors 
    anhydrous_majors = majors[1:end .!= findfirst(x->x==:Volatiles, majors)]
    result = Array{Float64}(undef, length(out_ME[1]), length(anhydrous_majors))
    rows = reshape(string.(anhydrous_majors), 1, :)
    for i in eachindex(anhydrous_majors)
        result[:,i] .= out_ME[anhydrous_majors[i]]
    end
    writedlm("output/GloRiSe_major_screened.tsv", vcat(rows, result))

    # Minors 
    result = Array{Float64}(undef, length(out_TE[1]), length(minors))
    rows = reshape(string.(minors), 1, :)
    for i in eachindex(minors)
        result[:,i] .= out_TE[minors[i]]
    end
    writedlm("output/GloRiSe_minor_screened.tsv", vcat(rows, result))


## --- End of File 