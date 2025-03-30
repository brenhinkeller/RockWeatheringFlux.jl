## --- Set up 
    # After Figures 2 and 3 from Rudnick and Gao, 2014. Comparison of major element 
    # estimates for different models 

    # Load data and base packages
    include("Definitions.jl")

    # List of (major) elements 
    majors, minors = get_elements()
    anhydrous_majors = majors[1:end .!= findfirst(x->x==:Volatiles, majors)]


## --- Load data 
    # This study
    ucc = importdataset(ucc_out, '\t', importas=:Tuple);
    ucc_err = importdataset(ucc_out_err, '\t', importas=:Tuple);
    elementindex = NamedTuple{Tuple(Symbol.(ucc.element))}(i for i in eachindex(ucc.element))
    ucc = NamedTuple{keys(class)}([NamedTuple{Tuple(majors)}(
        [ucc[f][elementindex[k]] ± ucc_err[f][elementindex[k]] for k in majors]) 
        for f in keys(class)]
    );

    # Eroded material 
    ersn = importdataset(comp_eroded, '\t', importas=:Tuple)
    elementindex = NamedTuple{Tuple(Symbol.(ersn.elem))}(i for i in eachindex(ersn.elem))
    ersn = NamedTuple{keys(class)}([NamedTuple{Tuple(majors)}(
        [ersn[f][elementindex[k]] for k in majors]) for f in keys(class)]
    );

    # Rudnick and Gao 2014; Condie 1993
    rudnick_gao = importdataset("data/rudnickgao14_table3_processed.csv",  ',', importas=:Tuple);
    elementindex = NamedTuple{Tuple(Symbol.(rudnick_gao.Element))}(i for i in eachindex(rudnick_gao.Element));
    rudnick_gao = NamedTuple{Tuple(majors)}(k in keys(elementindex) ?
        rudnick_gao.Upper_crust[elementindex[k]] ± rudnick_gao.SE_1_Sigma[elementindex[k]]*2 : NaN ± NaN
        for k in majors
    )

    # Gao et al., 1998 (10.1016/S0016-7037(98)00121-5)
    gao = (SiO2=58.48,Al2O3=12.12,FeOT=4.6,TiO2=0.57,MgO=3.73,CaO=7.41,Na2O=2.57,K2O=2.27,
        Volatiles=(1.88+5.48),
    )

    # Condie, 1993 [map model] (10.1016/0009-2541(93)90140-E)
    condie = (SiO2=66.21,Al2O3=14.96,FeOT=4.70,TiO2=0.55,MgO=2.42,CaO=3.60,Na2O=3.51,K2O=2.73,
        Volatiles=0,
    )

    # Pease et al., 2023 (10.1029/2023JB026353)
    pease = (SiO2=59.5,Al2O3=17.3,FeOT=6.53,TiO2=0.091,MgO=2.88,CaO=5.39,Na2O=4.14,K2O=3.24,
        Volatiles=(0.0673 + 0.0615),
    )

    # GloRiSe ± 2 s.e.
    glorise = importdataset("output/GlobalRivers/GloRiSe_major_screened.tsv", '\t', importas=:Tuple);
    glorise = NamedTuple{Tuple(majors)}(k in keys(glorise) ? 
        (nanmean(glorise[k]) ± 2*nansem(glorise[k])) : NaN ± NaN for k in majors
    );
    
    # # GLORICH ± 2 s.e.
    # glorich = importdataset("output/GlobalRivers/GLORICH_major_screened.tsv", '\t', importas=:Tuple);
    # glorich = NamedTuple{Tuple(anhydrous_majors)}(
    #     (nanmean(glorich[k]) ± 2*nansem(glorich[k])) for k in anhydrous_majors
    # );


## --- Normalize all values to my whole earth UCC values
    ratio = (
        ucc = [ucc.bulk[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        ersn = [ersn.bulk[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        plut = [ucc.plut[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        volc = [ucc.volc[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        sed = [ucc.sed[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        rudnick_gao = [rudnick_gao[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        condie = [condie[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        gao = [gao[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        pease = [pease[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        glorise = [glorise[k] ./ ucc.bulk[k] for k in anhydrous_majors],
        # glorich = [glorich[k] ./ ucc.bulk[k] for k in anhydrous_majors],
    )


## --- Base plot
    h_base = plot(
        framestyle=:box,
        fontfamily=:Helvetica, 
        xlims=(0.5, length(anhydrous_majors)+0.5),
        xticks=(1:length(anhydrous_majors), string.(anhydrous_majors)), 
        ylims=(0.3, 2.5),
        yticks=(0.5:0.5:2, string.(((0.5:0.5:2)))),
        fg_color_legend=:transparent, bg_color_legend=:transparent, 
        legendfontsize=10, legend=:topleft,
        labelfontsize=12, tickfontsize=10,
        ylabel="Normalized to This Study",
        titleloc=:left,
        grid=false,
        size=(600,400),
    )

    # Plot data for this study and ± 20%, 10%
    x = 1:length(anhydrous_majors)
    plot!(h_base, 
        [xlims(h_base)[1], xlims(h_base)[2], xlims(h_base)[2], xlims(h_base)[1]], 
        [0.75, 0.75, 1.25, 1.25],
        seriestype=:shape, color=:lightgrey, lcolor=:match, label=""
    )
    plot!(h_base, 
        [xlims(h_base)[1], xlims(h_base)[2], xlims(h_base)[2], xlims(h_base)[1]], 
        [0.85, 0.9, 1.1, 1.1],
        seriestype=:shape, color=:darkgrey, lcolor=:match, label=""
    )
    plot!(h_base, collect(xlims(h_base)), [1,1], label="", color=:black, linewidth=3)

    plot!(h_base, x, Measurements.value.(ratio.ucc), 
        markershape=:circle, 
        label="Exposed Continental Crust",
        color=:black, msc=:auto, 
        markersize=6, linewidth=2,
    )


## --- Continental crust 
    h_rox = deepcopy(h_base)
    title!(h_rox, "A. Exposed Continental Crust")

    # Errors from all plotted series,
    plot!(h_rox, ratio.plut,
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors.plut, 
        mcolor=colors.plut, markershape=:none,
    )
    plot!(h_rox, ratio.rudnick_gao,
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors_source.rudnick, 
        mcolor=colors_source.rudnick, markershape=:none,
    )
    # plot!(h_rox, ratio.pease,
    #     label="", 
    #     seriestype=:scatter, markersize=0,
    #     msw=2, msc=colors_source.pease, 
    #     mcolor=colors_source.pease, markershape=:none,
    # )
    plot!(h_rox, ratio.gao,
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors_source.gao, 
        mcolor=colors_source.gao, markershape=:none,
    )
    plot!(h_rox, ratio.condie,
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors_source.condie, 
        mcolor=colors_source.condie, markershape=:none,
    )

    # Values from this study 
    plot!(h_rox, x, Measurements.value.(ratio.plut), 
        markershape=:circle, 
        label="Plutonic Rocks",
        color=colors.plut, msc=:auto, 
        markersize=6, linewidth=2
    )


    # Plot previous estimates
    # plot!(h_rox, x, Measurements.value.(ratio.pease), 
    #     markershape=:diamond, 
    #     label="Pease et al., 2023",
    #     msc=colors_source.pease, 
    #     linecolor=colors_source.pease,
    #     markercolor=:white, markersize=7, 
    #     mswidth=2, linewidth=2
    # )
    plot!(h_rox, x, Measurements.value.(ratio.rudnick_gao), 
        markershape=:diamond, 
        label="Rudnick and Gao, 2014",
        msc=colors_source.rudnick,
        linecolor=colors_source.rudnick,
        markercolor=:white, markersize=7, 
        mswidth=2, linewidth=2
    )
    plot!(h_rox, x, Measurements.value.(ratio.gao), 
        markershape=:diamond, 
        label="Gao et al., 1998",
        msc=colors_source.gao, 
        linecolor=colors_source.gao,
        markercolor=:white, markersize=7, 
        mswidth=2, linewidth=2
    )
    plot!(h_rox, x, Measurements.value.(ratio.condie), 
        markershape=:diamond, 
        label="Condie, 1993",
        msc=colors_source.condie, 
        linecolor=colors_source.condie,
        markercolor=:white, markersize=7, 
        mswidth=2, linewidth=2
    )
    
    # # Plot the Pease TiO2 value that goes off the screen 
    # i = 4 
    # TiO2 = round(Measurements.value.(ratio.pease[i]), sigdigits=2)
    # annotate!(h_rox, x[i], 0.5, ("$TiO2", colors_source.pease, 10, :bottom))
    # annotate!(h_rox, x[i], 0.5, ("l", colors_source.pease, :top))
    
    display(h_rox)
    savefig(h_rox, "$filepath/modelcomparison_surface.pdf")


## --- Eroded sediment 
    h_ersn = deepcopy(h_base)
    title!(h_ersn, "B. Weathered Material")

    # Errors
    plot!(h_ersn, ratio.ersn,
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors.sed, 
        mcolor=colors.sed, markershape=:none,
    ) 
    plot!(h_ersn, ratio.glorise,
        label="", 
        seriestype=:scatter, markersize=0,
        msw=2, msc=colors_source.glorise, 
        mcolor=colors_source.glorise, markershape=:none,
    )

    # Values 
    plot!(h_ersn, x, Measurements.value.(ratio.ersn), 
        markershape=:circle, 
        label="Weathered Material",
        color=colors.sed, msc=:auto, 
        markersize=6, linewidth=2,
    )
    plot!(h_ersn, x, Measurements.value.(ratio.glorise), 
        markershape=:dtriangle, 
        label="Muller et al., 2021 (Suspended Sediment)",
        color=:cadetblue, msc=:cadetblue, 
        markercolor=:white, markersize=5, 
        mswidth=2, linewidth=2
    )

    display(h_ersn)
    savefig(h_ersn, "$filepath/modelcomparison_erosion.pdf")


## --- Slam that bad boy together 
    bigH = plot(h_rox, h_ersn, layout=(2,1), size=(600, 800))
    savefig(bigH, "$filepath/modelcomparison.pdf")
    display(bigH)

    
## --- End of file 