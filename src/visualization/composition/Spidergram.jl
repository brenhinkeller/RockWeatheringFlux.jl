## --- Set up 
    # Spider... gram. Spidergram.
    # Check out the REE patterns on this guy

    # Load data and base packages
    include("Definitions.jl")
    using Measurements

    # Get a list of REEs
    REEs = get_REEs()
    spider_REEs = REEs[1:end .!= findfirst(x->x==:Pm, REEs)]       # Pm isn't in datasets


## --- Load data
    @info "Upper crust data: $ucc_out"
    ucc = importdataset(ucc_out, '\t', importas=:Tuple)
    ucc_errs = importdataset(ucc_out_err, '\t', importas=:Tuple)

    ersn = importdataset(comp_eroded, '\t', importas=:Tuple)
    ersn_errs = importdataset(comp_eroded_err, '\t', importas=:Tuple)

    rudnick_gao = importdataset("data/rudnickgao14_table3_processed.csv",  ',', importas=:Tuple);
    glorise = importdataset("output/GlobalRivers/GloRiSe_minor_screened.tsv", '\t', importas=:Tuple);


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

    # Rudnick and Gao, 2014 
    elementindex = NamedTuple{Tuple(Symbol.(rudnick_gao.Element))}(i for i in eachindex(rudnick_gao.Element))
    rudnick_gao_err = NamedTuple{Tuple(spider_REEs)}(
        rudnick_gao.SE_1_Sigma[elementindex[k]]*2 for k in spider_REEs
    )
    rudnick_gao = NamedTuple{Tuple(spider_REEs)}(
        rudnick_gao.Upper_crust[elementindex[k]] for k in spider_REEs
    )

    # Other estimates 
    gao = (La=31.0,Ce=58.7,Nd=26.1,Sm=4.45,Eu=1.08,Tb=0.69,Yb=1.91,Lu=0.30)         # Carbonate inclusive estimate 
    shaw = (La=32.3,Ce=65.6,Nd=25.9,Sm=4.61,Eu=0.937,Tb=0.481,Dy=2.9,Ho=0.62,Yb=1.47, Lu=0.233)     # Canadian Shield 
    condie = (La=28.4, Ce=57.5, Nd=25.6, Sm=4.59, Eu=1.05, Gd=4.21, Tb=0.66, Yb=1.91, Lu=0.32)      # Map model 

    
    # Just carbonates 
    gao_carb = (La=13.8,Ce=25.9,Nd=11.5,Sm=1.50,Eu=NaN,Tb=0.22,Yb=0.58,Lu=0.08)
    

    # # Martin et al. TTG 
    # martin_old_ttg = (La=35.3,Ce=61.7,Nd=25.8,Sm=4.2,Eu=1,Gd=3.2,Dy=1.8,Er=0.77,Yb=0.78,Lu=0.2)
    # margin_young_ttg = (La=30.8,Ce=58.5,Nd=23.2,Sm=3.5,Eu=0.6,Gd=2.3,Dy=1.6,Er=0.75,Yb=0.63,Lu=0.12)

    # # Condie by lithology and time (Late Archean, Middle Proterozoic, and Meso-Cenozoic)
    # condiearcheanshale=(La=30.7,Ce=60.9,Nd=27.7,Sm=4.85,Eu=1.12,Gd=4.55,Tb=0.71,Yb=2.43,Lu=0.39)
    # condiearcheansed = (La=26,Ce=52,Nd=22,Sm=3.9,Eu=1.1,Gd=3.69,Tb=0.58,Yb=1.4,Lu=0.25,)
    # condieproterozoicsed = (La=28,Ce=60,Nd=26,Sm=4.9,Eu=0.93,Gd=4.34,Tb=0.66,Yb=2.2,Lu=0.38,)
    # condiemidphansed = (La=28,Ce=61,Nd=26,Sm=4.9,Eu=0.9,Gd=4.34,Tb=0.66,Yb=2.2,Lu=0.38,)

    # Suspended load, convert to ppm
    glorise_err = NamedTuple{Tuple(spider_REEs)}(2 * nanstd(glorise[k].*10_000) / sqrt(count(.!isnan.(glorise[k]))) for k in spider_REEs)
    glorise = NamedTuple{Tuple(spider_REEs)}(nanmean(glorise[k].*10_000) for k in spider_REEs)
        

## --- Base plot (surface earth)
    h = spidergram(ucc.bulk, 
        label="Exposed Continental Crust",
        seriescolor=:black, msc=:auto,
        markershape=:circle, markersize=5,
        fontfamily=:Helvetica, 
        legend=false,
        titleloc=:left,
        legendfont=10, titlefont=16,
        ylims=(2, 200),
        yticks=([10, 100], ["10", "100"]),
        size=(550,400), 
        # left_margin=(15,:px),
    );

    # Error bars are smaller than data points for all values
    # spidergram!(h, 
    #     collect(values(ucc.bulk) .± values(ucc_err.bulk)),
    #     label="", 
    #     seriestype=:scatter, markersize=0,
    #     msw=2, msc=:black, 
    #     mcolor=:black, markershape=:none,
    # )


## --- Surface Earth
    h_surf = deepcopy(h)
    title!(h_surf, "A. Average Continental Crust")

    # Error bars
    spidergram!(h_surf, 
        collect(values(rudnick_gao) .± values(rudnick_gao_err)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors_source.rudnick, 
        mcolor=colors_source.rudnick, markershape=:none,
    )

    # Values 
    spidergram!(h_surf, shaw, label="Canadian Shield (Shaw et al., 1967, 1976)", 
        seriescolor=colors_source.shaw, msc=colors_source.shaw, markercolor=:white,
        markershape=:diamond, markersize=5,
        mswidth=1.5, linewidth=1.5,
    )
    spidergram!(h_surf, condie, label="Exposed Crust (Condie, 1993)", 
        seriescolor=colors_source.condie, msc=colors_source.condie, markercolor=:white,
        markershape=:diamond, markersize=5,
        mswidth=1.5, linewidth=1.5,
    )
    spidergram!(h_surf, gao, label="Exposed crust (Gao et al., 1998)", 
        seriescolor=colors_source.gao, msc=colors_source.gao, markercolor=:white,
        markershape=:diamond, markersize=5,
        mswidth=1.5, linewidth=1.5,
    )
    spidergram!(h_surf, rudnick_gao, 
        label="Upper crust (Rudnick and Gao, 2014)", 
        seriescolor=colors_source.rudnick, msc=colors_source.rudnick, markercolor=:white,
        markershape=:diamond, markersize=5,
        mswidth=1.5, linewidth=1.5,        
    )
    

## --- Selected lithologies
    h_lith = deepcopy(h)
    title!(h_lith, "B. Selected Lithologies")

    # Error bars 
    spidergram!(h_lith, 
        collect(values(ucc.carb) .± values(ucc_err.carb)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors.carb, 
        mcolor=colors.carb, markershape=:none,
    )
    spidergram!(h_lith, 
        collect(values(ucc.shale) .± values(ucc_err.shale)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors_dark.shale, 
        mcolor=colors_dark.shale, markershape=:none,
    )
    spidergram!(h_lith, 
        collect(values(ucc.trondhjemite) .± values(ucc_err.trondhjemite)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors.trondhjemite, 
        mcolor=colors.trondhjemite, markershape=:none,
    )
    spidergram!(h_lith, 
        collect(values(ucc.tonalite) .± values(ucc_err.tonalite)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors.tonalite, 
        mcolor=colors.tonalite, markershape=:none,
    )
    spidergram!(h_lith, 
        collect(values(ucc.granodiorite) .± values(ucc_err.granodiorite)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors.granodiorite, 
        mcolor=colors.granodiorite, markershape=:none,
    )

    # Values 
    # spidergram!(h_lith, ucc.granite, 
    #     label="Granite",
    #     seriescolor=colors.plut, msc=:auto, 
    #     markershape=:circle, markersize=5,
    # )
    spidergram!(h_lith, ucc.shale, 
        label="Shale", 
        seriescolor=colors_dark.shale, msc=:auto,
        markershape=:circle, markersize=5,
    )
    spidergram!(h_lith, ucc.carb, 
        label="Carbonate", 
        seriescolor=colors.carb, msc=:auto,
        markershape=:circle, markersize=5,
    )
    spidergram!(h_lith, ucc.trondhjemite, 
        label="Trondhjemite",
        seriescolor=colors.trondhjemite, msc=:auto,
        markershape=:circle, markersize=5,
    )
    spidergram!(h_lith, ucc.tonalite, 
        label="Tonalite",
        seriescolor=colors.tonalite, msc=:auto,
        markershape=:circle, markersize=5,
    )
    spidergram!(h_lith, ucc.granodiorite, 
        label="Granodiorite",
        seriescolor=colors.granodiorite, msc=:auto,
        markershape=:circle, markersize=5,
    )


## --- Eroded material 
    h_ersn = deepcopy(h)
    title!(h_ersn, "C. Weathered Material")

    # Error bars (smaller than point size for all values)
    spidergram!(h_ersn, 
        collect(values(glorise) .± values(glorise_err)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors_source.glorise, 
        mcolor=colors_source.glorise, markershape=:none,
    )
    spidergram!(h_ersn, 
        collect(values(ersn.bulk) .± values(ersn_err.bulk)),
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors.sed, 
        mcolor=colors.sed, markershape=:none,
    )

    # Values 
    spidergram!(h_ersn, glorise, 
        label="Suspended Sediment (Muller et al., 2021)", 
        seriescolor=colors_source.glorise, msc=colors_source.glorise, markercolor=:white,
        markershape=:diamond, markersize=5,
        mswidth=1.5, linewidth=1.5,
    )
    spidergram!(h_ersn, 
        ersn.bulk, 
        label="Eroded Sediment", 
        seriescolor=colors.sed, msc=:auto,
        markershape=:circle, markersize=5,
    )
    
    
## --- Communal legend?? 
    h_leg = plot(
        fg_color_legend=:transparent, 
        bg_color_legend=:transparent,
        legendfont=12, 
        legend=:inside,
        grid=false, axis=false,
        xlims=(2,3), ylims=(2,3),
        size=(550,400), 
    )

    plot!(h_leg, [1],[1], label="A. Average Continental Crust", color=:white)
    plot!(h_leg, [1],[1], label="Exposed Continental Crust", 
        seriescolor=:black, msc=:auto, 
        markershape=:circle, markersize=5,)
    plot!(h_leg, [1],[1], label="Rudnick and Gao, 2014 (Upper Continental Crust)", 
        seriescolor=colors_source.rudnick, msc=colors_source.rudnick, markercolor=:white,
        msw=1, markershape=:diamond, markersize=5,)
    plot!(h_leg, [1],[1], label="Gao et al., 1998 (Exposed Continental Crust)", 
        seriescolor=colors_source.gao, msc=colors_source.gao, markercolor=:white,
        msw=1, markershape=:diamond, markersize=5,)
    plot!(h_leg, [1],[1], label="Condie, 1993 (Exposed Continental Crust)", 
        seriescolor=colors_source.condie, msc=colors_source.condie, markercolor=:white,
        msw=1, markershape=:diamond, markersize=5,)
    plot!(h_leg, [1],[1], label="Shaw et al., 1967, 1976 (Canadian Shield)", 
        seriescolor=colors_source.shaw, msc=colors_source.shaw, markercolor=:white,
        msw=1, markershape=:diamond, markersize=5,)

    plot!(h_leg, [1],[1], label="\t", color=:white)
    plot!(h_leg, [1],[1], label="B. Selected Lithologies", color=:white)
    plot!(h_leg, [1],[1], label="Shale", 
        seriescolor=colors_dark.shale, msc=:auto, 
        markershape=:circle, markersize=5,) 
    plot!(h_leg, [1],[1], label="Carbonate", 
        seriescolor=colors.carb, msc=:auto, 
        markershape=:circle, markersize=5,) 
    plot!(h_leg, [1],[1], label="Trondhjemite", 
        seriescolor=colors.trondhjemite, msc=:auto, 
        markershape=:circle, markersize=5,) 
    plot!(h_leg, [1],[1], label="Tonalite", 
        seriescolor=colors.tonalite, msc=:auto, 
        markershape=:circle, markersize=5,)
    plot!(h_leg, [1],[1], label="Granodiorite", 
        seriescolor=colors.granodiorite, msc=:auto, 
        markershape=:circle, markersize=5,) 

    plot!(h_leg, [1],[1], label="\t", color=:white)
    plot!(h_leg, [1],[1], label="C. Weathered Material", color=:white)
    plot!(h_leg, [1],[1], label="Weathered Material", 
        seriescolor=colors.sed, msc=:auto, 
        markershape=:circle, markersize=5,) 
    plot!(h_leg, [1],[1], label="Muller et al., 2021 (Suspended Sediment)", 
        seriescolor=colors_source.glorise, msc=colors_source.glorise, markercolor=:white,
        msw=1, markershape=:diamond, markersize=5,)


## --- Save files 
    display(h_surf)
    display(h_lith)
    display(h_ersn)

    savefig(h_surf, "$filepath/spidergram_surface.pdf")
    savefig(h_lith, "$filepath/spidergram_lithologies.pdf")
    savefig(h_ersn, "$filepath/spidergram_eroded.pdf")
    savefig(h_leg, "$filepath/spidergram_legend.pdf")

    # Assemble plots, but this is a placeholder because the y axis gets all messed up :(
    bigH = Plots.plot(h_surf, h_lith, h_ersn, h_leg, layout=(2, 2), size=(1100, 800))
    display(bigH)
    savefig(bigH, "$filepath/spidergram.pdf")


## --- End of file 