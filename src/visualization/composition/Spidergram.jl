## --- Set up 
    # Spider... gram. Spidergram.
    # Check out the REE patterns on this guy

    # Load data and base packages
    include("Definitions.jl")
    using Measurements

    # Get a list of REEs
    REEs = get_REEs()
    spider_REEs = REEs[1:end .!= findfirst(x->x==:Pm, REEs)]       # Pm isn't in datasets

    # Load data
    @info "Upper crust data: $ucc_out"
    ucc = importdataset(ucc_out, '\t', importas=:Tuple)
    ucc_errs = importdataset(ucc_out_err, '\t', importas=:Tuple)

    ersn = importdataset(comp_eroded, '\t', importas=:Tuple)
    ersn_errs = importdataset(comp_eroded_err, '\t', importas=:Tuple)

    rudnick_gao = importdataset("data/rudnickgao2014.csv",  ',', importas=:Tuple);
    GloRiSe = importdataset("output/GlobalRivers/GloRiSe_minor_screened.tsv", '\t', importas=:Tuple);


## --- Load data 
    # Matched data, converting to ppm
    elementindex = NamedTuple{Tuple(Symbol.(ucc.element))}(i for i in eachindex(ucc.element))
    ucc = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ucc[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
    );
    ucc_err = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ucc_errs[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
    );

    ersn = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ersn[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
    );
    ersn_err = NamedTuple{keys(class)}([NamedTuple{Tuple(spider_REEs)}(
        [ersn_errs[f][elementindex[k]].*10_000 for k in spider_REEs]) for f in keys(class)]
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

    # Gao et al., 1998 carbonate inclusive estimate 
    gao = (La=31.0,Ce=58.7,Nd=26.1,Sm=4.45,Eu=1.08,Tb=0.69,Yb=1.91,Lu=0.30)
    gao_carb = (La=13.8,Ce=25.9,Nd=11.5,Sm=1.50,Eu=NaN,Tb=0.22,Yb=0.58,Lu=0.08)

    # Condie by lithology and time (Late Archean, Middle Proterozoic, and Meso-Cenozoic)
    condiearcheanshale=(La=30.7,Ce=60.9,Nd=27.7,Sm=4.85,Eu=1.12,Gd=4.55,Tb=0.71,Yb=2.43,Lu=0.39)
    condiearcheansed = (La=26,Ce=52,Nd=22,Sm=3.9,Eu=1.1,Gd=3.69,Tb=0.58,Yb=1.4,Lu=0.25,)
    condieproterozoicsed = (La=28,Ce=60,Nd=26,Sm=4.9,Eu=0.93,Gd=4.34,Tb=0.66,Yb=2.2,Lu=0.38,)
    condiemidphansed = (La=28,Ce=61,Nd=26,Sm=4.9,Eu=0.9,Gd=4.34,Tb=0.66,Yb=2.2,Lu=0.38,)

    # Suspended load, convert to ppm
    GloRiSe = NamedTuple{Tuple(spider_REEs)}(nanmean(GloRiSe[k].*10_000) for k in spider_REEs)


## --- Calculate REE abundance for TTGs
    ttg = match_cats.trondhjemite .| match_cats.tonalite .| match_cats.granodiorite;
    npoints = NamedTuple{Tuple(spider_REEs)}(
        [unique_sample(mbulk.Sample_ID[.!isnan.(mbulk[k])], 90) for k in spider_REEs]
    );

    ucc = merge(ucc, (;
        ttg = NamedTuple{Tuple(spider_REEs)}(
            [nanmean(mbulk[j][ttg])*10_000 for j in spider_REEs])
    ));
    ucc_err = merge(ucc_err, (;
        ttg = NamedTuple{Tuple(spider_REEs)}(
            [nanstd(mbulk[k][ttg])*10_000 ./sqrt(npoints[k]).*2 for k in spider_REEs])
    ));
        

## --- Assemble plots
    # Base plot (surface earth)
    # Error bars smaller than markers
    h = spidergram(ucc.bulk, 
        # collect(values(ucc.bulk) .± values(ucc_err.bulk)),
        label="Exposed Continental Crust",
        seriescolor=:black, msc=:auto,
        markershape=:circle, markersize=5,
        fontfamily=:Helvetica, 
        legend=:topright,
        titleloc=:left,
        legendfont=10, titlefont=16,
        ylims=(3,175),
        size=(550,400), 
        left_margin=(15,:px),
    );

    # Surface Earth
    h_surf = deepcopy(h)
    spidergram!(h_surf, rudnick_gao, label="Rudnick and Gao, 2014", 
        seriescolor=colors_source.rudnick, msc=:auto,
        markershape=:diamond, markersize=6,
        title="A. Exposed / Upper Continental Crust",
    )
    spidergram!(h_surf, gao, label="Gao et al., 1998", 
        seriescolor=colors_source.gao, msc=:auto,
        markershape=:diamond, markersize=6,
    )
    
    # Selected lithologies
    h_lith = deepcopy(h)
    spidergram!(h_lith, ucc.granite, 
        # collect(values(ucc.granite) .± values(ucc_err.granite)),
        label="Granite",
        seriescolor=colors.plut, msc=:auto, 
        markershape=:circle, markersize=5,
        title="B. Selected Lithologies",
    )
    spidergram!(h_lith, ucc.basalt, 
        # collect(values(ucc.basalt) .± values(ucc_err.basalt)),
        label="Basalt",
        seriescolor=:sandybrown, msc=:auto, 
        markershape=:circle, markersize=5
    )
    spidergram!(h_lith, ucc.ttg, 
        # collect(values(ucc.ttg) .± values(ucc_err.ttg)),
        label="TTG",
        seriescolor=:purple, msc=:auto, 
        markershape=:circle, markersize=5
    )
    spidergram!(h_lith, ucc.shale, 
        # collect(values(ucc.shale) .± values(ucc_err.shale)),
        label="Shale", 
        seriescolor=colors_dark.shale, msc=:auto,
        markershape=:circle, markersize=5,
    )

    # Eroded material 
    h_ersn = deepcopy(h)
    spidergram!(h_ersn, GloRiSe, label="Suspended Sediment (Mueller et al.)", 
        seriescolor=colors_source.muller, msc=:auto,
        markershape=:diamond, markersize=6,
        title="C. Eroded Material"
    )
    spidergram!(h_ersn, ersn.bulk, 
        # collect(values(ersn.bulk) .± values(ersn_err.bulk)),
        label="Eroded Sediment", 
        seriescolor=colors_source.this_study, msc=:auto,
        markershape=:circle, markersize=5,
    )
    
    display(h_surf)
    display(h_lith)
    display(h_ersn)


## --- Save files 
    savefig(h_surf, "$filepath/spidergram_surface.pdf")
    savefig(h_lith, "$filepath/spidergram_lithologies.pdf")
    savefig(h_ersn, "$filepath/spidergram_eroded.pdf")

    # Assemble plots, but this is a placeholder because the y axis gets all messed up :(
    bigH = Plots.plot(h_surf, h_lith, h_ersn, layout=(1, 3), size=(1800, 400))
    display(bigH)
    savefig(bigH, "$filepath/spidergram.pdf")


## --- Look at all lithologies
    include("spidergramtest.jl")
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];
    cats = (
        sed = minorsed, 
        volc = (minorvolc[collect(minorvolc) .!= :volcaniclast]..., :carbonatite),
        plut = minorplut,
    )

    m = maximum(length.([cats[k] for k in keys(cats)])) + 1
    p = palette(:jet1, m)

    n = chondrite
    h = spidergram(ucc.bulk, normalizer=n,
        # collect(values(ucc.bulk) .± values(ucc_err.bulk)),
        label="Surface Earth",
        seriescolor=:black, msc=:auto,
        markershape=:circle, markersize=5,
        fontfamily=:Helvetica, 
        legend=:outerright, titleloc=:left,
        legendfont=8, titlefont=16,
        ylims=(1,500),
        # ylims=(10^-2, 10^1),
        # yticks=(10.0 .^(-2:1), ["0.01", "0.1", "1", "10"]),
        size=(700,400), 
        left_margin=(15,:px),
    )

    for k in keys(cats)
        hprime = deepcopy(h)
        p = palette(:jet1, length(cats[k])+1)

        spidergram!(hprime, ucc[k], normalizer=n,
            label="$k",
            color=p[1], msc=:auto,
        )
        for i in eachindex(cats[k])
            kprime = cats[k][i] 
            spidergram!(hprime, ucc[kprime], normalizer=n,
                label="$kprime",
                color=p[i+1], msc=:auto,
            )
        end
        display(hprime)
    end


## --- End of file 