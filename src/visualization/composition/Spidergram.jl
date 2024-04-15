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
    ersn = importdataset(erodedcomp_out, '\t', importas=:Tuple);
    rudnick_gao = importdataset("data/rudnick_gao_2014_table1-2.csv",  ',', importas=:Tuple);
    GloRiSe = importdataset("output/GloRiSe_minor_screened.tsv", '\t', importas=:Tuple);


## --- Load data 
    # Matched data, converting to ppm
    elementindex = NamedTuple{Tuple(Symbol.(ucc.element))}(i for i in eachindex(ucc.element))
    ucc = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ucc[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
    );
    ersn = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ersn[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
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

    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, age_error, age_error_abs  
    )

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
    

## --- Assemble plots
    # Base plot (surface earth)
    h = spidergram(ucc.bulk, label="Surface Earth (This Study)",
        seriescolor=:black, msc=:auto,
        markershape=:circle, markersize=5,
        fontfamily=:Helvetica, 
        legend=:topright, titleloc=:left,
        legendfont=10, titlefont=16,
        ylims=(3,175),
        size=(700,400), 
        left_margin=(15,:px),
    );

    # Surface Earth
    h1 = deepcopy(h)
    spidergram!(h1, rudnick_gao, label="Rudnick and Gao", 
        seriescolor=colors_source.rudnick, msc=:auto,
        markershape=:diamond, markersize=6,
        title="A. Surface Earth Average",
    )
    spidergram!(h1, condie, label="Condie", 
        seriescolor=colors_source.condie, msc=:auto,
        markershape=:diamond, markersize=6,
    )
    spidergram!(h1, gao, label="Gao et al.", 
        seriescolor=colors_source.gao, msc=:auto,
        markershape=:diamond, markersize=6,
    )
    spidergram!(h1, GloRiSe, label="Suspended Sediment", 
        seriescolor=colors_source.muller, msc=:auto,
        markershape=:diamond, markersize=6,
    )
    spidergram!(h1, ersn.bulk, label="Eroded Sediment (This Study)", 
        seriescolor=colors_source.this_study, msc=:auto,
        markershape=:circle, markersize=5,
    )
    savefig("$filepath/spidergram_bulk.pdf")

    # Igneous 
    h2 = deepcopy(h)
    spidergram!(h2, ucc.ign, label="All Igneous",
        seriescolor=colors.ign, msc=:auto, 
        markershape=:circle, markersize=5,
        title="B. Igneous Rocks",
    )
    spidergram!(h2, ucc.granite, label="Granite",
        seriescolor=colors.plut, msc=:auto, 
        markershape=:circle, markersize=5
    )
    spidergram!(h2, ucc.basalt, label="Basalt",
        seriescolor=colors.basalt, msc=:auto, 
        markershape=:circle, markersize=5
    )
    savefig("$filepath/spidergram_ign.pdf")

    # Sedimentary
    h3 = deepcopy(h)
    spidergram!(h3, ersn.bulk, label="Eroded Sediment", 
        seriescolor=colors_source.this_study, msc=:auto,
        markershape=:circle, markersize=5,
        title="C. Sedimentary Rocks",
    )
    spidergram!(h3, ucc.sed, label="Sedimentary", 
        seriescolor=colors_dark.sed, msc=:auto,
        markershape=:circle, markersize=5,
    )
    spidergram!(h3, ucc.shale, label="Shale", 
        seriescolor=colors_dark.shale, msc=:auto,
        markershape=:circle, markersize=5,
    )
    spidergram!(h3, ucc.carb, label="Carbonate", 
        seriescolor=colors.carb, msc=:auto,
        markershape=:circle, markersize=5,
    )
    savefig("$filepath/spidergram_sed.pdf")

    # Shales 
    h4 = deepcopy(h)
    spidergram!(h4, archeanshale, label="Archean Shale", 
        seriescolor=parse(Colorant, "#007a6c"), msc=:auto,
        markershape=:circle, markersize=5,
        title="D. Fine Grained Sedimentary Rocks",
    )
    spidergram!(h4, condiearcheansed, label="Archean Graywacke (Condie)",
        seriescolor=colors_source.condie, msc=:auto,
        markershape=:diamond, markersize=6
    )
    spidergram!(h4, phanerozoicshale, label="Phanerozoic Shale", 
        seriescolor=parse(Colorant, "#f58220"), msc=:auto,
        markershape=:circle, markersize=5,
    )
    savefig("$filepath/spidergram_shale.pdf")
    # display(h3)

    # Assemble plots, but this is a placeholder because the y axis gets all messed up :(
    h = Plots.plot(h1, h2, h3, h4, layout=(2, 2), size=(1200,800))
    display(h)
    savefig(h, "$filepath/spidergram.pdf")


## --- End of file 