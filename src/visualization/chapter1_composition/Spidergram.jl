## --- Set up 
    # Spider... gram. Spidergram.
    # Check out the REE patterns on this guy

    # Load data and base packages
    include("Definitions.jl")

    # Get a list of REEs
    REEs = get_REEs()
    spider_REEs = REEs[1:end .!= findfirst(x->x==:Pm, REEs)]       # Pm isn't in datasets

    # Load data
    @info "Upper crust data: $ucc_out"
    ucc = importdataset(ucc_out, '\t', importas=:Tuple);
    ucc_err = importdataset(ucc_out_err, '\t', importas=:Tuple);
    rudnick_gao = importdataset("data/rudnick_gao_2014_table1-2.csv",  ',', importas=:Tuple);
    GloRiSe = importdataset("output/GloRiSe_minor_screened.tsv", '\t', importas=:Tuple);

    # Until StatGeochem PR gets updated on this machine 
    # include("../../../src/dev.jl")


## --- Load data 
    # Matched data, converting to ppm
    elementindex = NamedTuple{Tuple(Symbol.(ucc.element))}(i for i in eachindex(ucc.element))
    ucc = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ucc[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
    );
    ucc_err = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ucc_err[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
    );

    # Previous estimates, units are already ppm
    elementindex = NamedTuple{Tuple(Symbol.(rudnick_gao.Element))}(
        i for i in eachindex(rudnick_gao.Element)
    )
    condie = NamedTuple{Tuple(spider_REEs)}(rudnick_gao.Condie_1993[elementindex[k]] 
        for k in spider_REEs)
    taylor_mclennan = NamedTuple{Tuple(spider_REEs)}(rudnick_gao.Taylor_and_McLennan_1985[elementindex[k]] 
        for k in spider_REEs)
    rudnick_gao = NamedTuple{Tuple(spider_REEs)}(rudnick_gao.This_Study[elementindex[k]] 
        for k in spider_REEs
    );

    # Gao et al., 1998 carbonate free estimate 
    gao = (La=31.0,Ce=58.7,Nd=26.1,Sm=4.45,Eu=1.08,Tb=0.69,Yb=1.91,Lu=0.30)

    # Condie by lithology and time (Late Archean, Middle Proterozoic, and Meso-Cenozoic)
    condiearcheanshale=(La=30.7,Ce=60.9,Nd=27.7,Sm=4.85,Eu=1.12,Gd=4.55,Tb=0.71,Yb=2.43,Lu=0.39)
    condiearcheansed = (La=26,Ce=52,Nd=22,Sm=3.9,Eu=1.1,Gd=3.69,Tb=0.58,Yb=1.4,Lu=0.25,)
    condieproterozoicsed = (La=28,Ce=60,Nd=26,Sm=4.9,Eu=0.93,Gd=4.34,Tb=0.66,Yb=2.2,Lu=0.38,)
    condiemidphansed = (La=28,Ce=61,Nd=26,Sm=4.9,Eu=0.9,Gd=4.34,Tb=0.66,Yb=2.2,Lu=0.38,)

    # Suspended load, convert to ppm
    GloRiSe = NamedTuple{Tuple(spider_REEs)}(nanmean(GloRiSe[k].*10_000) for k in spider_REEs)


## --- Resample (temporal) Archean vs. post-Archean shales 
    # Resample 
    nsims = Int(1e7)
    err = 0.01
    age_error = 0.05
    simout = Array{Float64}(undef, nsims, length(spider_REEs)+1)

    t = @. !isnan(mbulk.Age);
    sampleage = copy(mbulk.Age);
    ageuncert = nanadd.(mbulk.Age_Max, .- mbulk.Age_Min) ./ 2;
    sampleage[t] .= macrostrat.age[t]
    ageuncert[t] .= nanadd.(macrostrat.agemax[t], .- macrostrat.agemin[t]) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = max(sampleage[i]*age_error, ageuncert[i])
    end

    t = @. !isnan.(sampleage) .& match_cats.shale;
    data = Array{Float64}(undef, count(t), length(spider_REEs)+1)
    uncertainty = similar(data)

    k = invweight_age(sampleage[t])
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    [data[:,i] .= mbulk[spider_REEs[i]][t].*10_000 for i in eachindex(spider_REEs)]
    data[:,end] .= sampleage[t]
    uncertainty[:,1:end-1] .= err
    uncertainty[:,end] .= ageuncert[t]

    simout .= bsresample(data, uncertainty, nsims, p)

    # Create NamedTuples of Archean vs. post-Archean shales 
    archean = 3800 .>= simout[:,end] .> 2500;
    proterozoic = 2500 .>= simout[:,end] .> 541;
    phanerozoic = 541 .>= simout[:,end] .> 0;

    archeanshale = NamedTuple{Tuple(spider_REEs)}(
        reshape(nanmean(simout[:,1:end-1][archean,:], dims=1), :, 1)
    )
    proterozoicshale = NamedTuple{Tuple(spider_REEs)}(
        reshape(nanmean(simout[:,1:end-1][proterozoic,:], dims=1), :, 1)
    )
    phanerozoicshale = NamedTuple{Tuple(spider_REEs)}(
        reshape(nanmean(simout[:,1:end-1][phanerozoic,:], dims=1), :, 1)
    )
    allshale = NamedTuple{Tuple(spider_REEs)}(
        reshape(nanmean(simout[:,1:end-1], dims=1), :, 1)
    )


## --- Terminal printout: what's the ratio of my estimate to Rudnick and Gao?
    me = [ucc.bulk[k] for k in spider_REEs]
    rg = [rudnick_gao[k] for k in spider_REEs]
    ratio = round.(me ./ rg, digits=2)

    @info """ REEs are higher are higher than Rudnick and Gao by a factor of:
    $(join(rpad.(spider_REEs, 8), " "))
    $(join(rpad.(ratio, 8), " "))

    Average: $(round(nanmean(ratio), digits=2))
    2σ s.d:  $(round(nanstd(ratio)*2, digits=2))
    """
    

## --- Assemble plots
    p = Plots.palette(colorpalette, 5)

    # Bulk Earth
    h1 = spidergram(rudnick_gao, label="Rudnick and Gao (Whole Earth)", 
        markershape=:diamond, seriescolor=p[1], msc=:auto, markersize=6,
        legend=:topright, legendfont=10, titlefont=16,
        size=(700,400), title="A. Whole Earth Averages", titleloc=:left,
        left_margin=(15,:px),
        fontfamily=:Helvetica, 
    )
    spidergram!(h1, condie, label="Condie (Whole Earth)",
        markershape=:diamond, seriescolor=p[2], msc=:auto, markersize=6,)
    spidergram!(h1, gao, label="Gao et al. (Central East China)",
        markershape=:diamond, seriescolor=p[3], msc=:auto, markersize=6,)
    spidergram!(h1, GloRiSe, label="Muller et al. (Suspended Sediment)",
        markershape=:diamond, seriescolor=p[4], msc=:auto, markersize=6,
        linestyle=:dot)
    spidergram!(h1, ucc.bulk, label="This Study",
        markershape=:circle, seriescolor=p[5], msc=:auto, markersize=5)
    ylims!(4,200)
    savefig("$filepath/spidergram_bulk.pdf")
    # display(h1)

    # Major lithologies 
    h4 = spidergram(ucc.bulk, label="Whole Earth Average",
    markershape=:circle, seriescolor=p[5], msc=:auto, markersize=5,
        legend=:topright, legendfont=10, titlefont=16,
        size=(700,400), title="B. Major Lithology Averages", titleloc=:left,
        left_margin=(15,:px),
        fontfamily=:Helvetica, 
    )
    spidergram!(h4, ucc.sed, label="This Study (Sedimentary)",
        markershape=:circle, seriescolor=p[1], msc=:auto, markersize=5,)
    spidergram!(h4, ucc.ign, label="This Study (Igneous)",
        markershape=:circle, seriescolor=p[2], msc=:auto, markersize=5,)
    ylims!(4,200)
    savefig("$filepath/spidergram_lith.pdf")

    # Igneous rocks
    h2 = spidergram(ucc.bulk, label="Whole Earth Average",
        markershape=:circle, seriescolor=p[5], msc=:auto, markersize=5,
        legend=:topright, legendfont=10, titlefont=16,
        size=(700,400), title="C. Igneous", titleloc=:left,
        left_margin=(15,:px),
        fontfamily=:Helvetica, 
    )
    spidergram!(h2, ucc.ign, label="All Igneous",
        markershape=:circle, seriescolor=p[1], msc=:auto, markersize=5)
    spidergram!(h2, ucc.granite, label="Granite",
        markershape=:circle, seriescolor=p[2], msc=:auto, markersize=5)
    spidergram!(h2, ucc.basalt, label="Basalt",
        markershape=:circle, seriescolor=p[3], msc=:auto, markersize=5)
    ylims!(4,200)
    savefig("$filepath/spidergram_igneous.pdf")
    # display(h2)

    # Shales
    h3 = spidergram(ucc.bulk, label="Whole Earth Average (This Study)",
        markershape=:circle, seriescolor=p[5], msc=:auto, markersize=5,
        legend=:topright, legendfont=10, titlefont=16,
        size=(700,400), title="D. Shales and Greywacke", titleloc=:left,
        left_margin=(15,:px),
        fontfamily=:Helvetica, 
    )
    spidergram!(h3, condiearcheanshale, label="Archean Shale (Condie)",
        markershape=:diamond, seriescolor=p[1], msc=:auto, markersize=6)
    spidergram!(h3, condiearcheansed, label="L. Archean Graywacke (Condie)",
        markershape=:diamond, seriescolor=p[2], msc=:auto, markersize=6)
    spidergram!(h3, archeanshale, label="Archean Shale (This Study)",
        markershape=:circle, seriescolor=p[3], msc=:auto, markersize=5)
    spidergram!(h3, phanerozoicshale, label="Phanerozoic Shale  (This Study)",
        markershape=:circle, seriescolor=p[4], msc=:auto, markersize=5)
    ylims!(4,200)
    savefig("$filepath/spidergram_shale.pdf")
    # display(h3)

    # Assemble plots, but this is a placeholder because the y axis gets all messed up :(
    h = Plots.plot(h1, h4, h2, h3, layout=(2, 2), size=(1200,800))
    display(h)
    savefig(h, "$filepath/spidergram.pdf")


## --- End of file 