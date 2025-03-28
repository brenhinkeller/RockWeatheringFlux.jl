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
    ucc = NamedTuple{keys(class)}([NamedTuple{Tuple(anhydrous_majors)}(
        [ucc[f][elementindex[k]] ± ucc_err[f][elementindex[k]] for k in anhydrous_majors]) 
        for f in keys(class)]
    );

    # Eroded material 
    ersn = importdataset(comp_eroded, '\t', importas=:Tuple)
    elementindex = NamedTuple{Tuple(Symbol.(ersn.elem))}(i for i in eachindex(ersn.elem))
    ersn = NamedTuple{keys(class)}([NamedTuple{Tuple(anhydrous_majors)}(
        [ersn[f][elementindex[k]] for k in anhydrous_majors]) for f in keys(class)]
    );

    # Rudnick and Gao 2014; Condie 1993
    # I don't know where the numbers for Condie in Rudnick and Gao come from, but Condie 
    # doesn't estimate volatiles so I guess that's fine to use the Rudnick and Gao version
    rudnick_gao = importdataset("data/rudnickgao2014.csv",  ',', importas=:Tuple);
    elementindex = NamedTuple{Tuple(Symbol.(rudnick_gao.Element))}(
        i for i in eachindex(rudnick_gao.Element))
    condie = NamedTuple{Tuple(anhydrous_majors)}((haskey(elementindex, k) ? 
        rudnick_gao.Condie_1993[elementindex[k]] : NaN) for k in anhydrous_majors
    );
    rudnick_gao = NamedTuple{Tuple(anhydrous_majors)}((haskey(elementindex, k) ? 
        rudnick_gao.This_Study[elementindex[k]] : NaN) for k in anhydrous_majors
    );

    # Gao et al., 1998 (10.1016/S0016-7037(98)00121-5)
    gao = (SiO2=58.48,Al2O3=12.12,FeOT=4.6,TiO2=0.57,MgO=3.73,CaO=7.41,Na2O=2.57,K2O=2.27,
        # Volatiles=7.36,
    )

    # Pease et al.
    pease = (SiO2=59.5,Al2O3=17.3,FeOT=6.53,TiO2=0.091,MgO=2.88,CaO=5.39,Na2O=4.14,K2O=3.24,
        # Volatiles=(0.0673 + 0.0615),
    )

    # GloRiSe ± 2 s.e.
    glorise = importdataset("output/GlobalRivers/GloRiSe_major_screened.tsv", '\t', importas=:Tuple);
    glorise = NamedTuple{Tuple(anhydrous_majors)}(
        (nanmean(glorise[k]) ± 2*nansem(glorise[k])) for k in anhydrous_majors
    );
    
    # GLORICH ± 2 s.e.
    glorich = importdataset("output/GlobalRivers/GLORICH_major_screened.tsv", '\t', importas=:Tuple);
    glorich = NamedTuple{Tuple(anhydrous_majors)}(
        (nanmean(glorich[k]) ± 2*nansem(glorich[k])) for k in anhydrous_majors
    );


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
        glorich = [glorich[k] ./ ucc.bulk[k] for k in anhydrous_majors],
    )


## --- Base plot
    h_base = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica, 
        xlims=(0.5, length(anhydrous_majors)+0.5),
        xticks=(1:length(anhydrous_majors), string.(anhydrous_majors)), 
        ylims=(0.3, 2.5),
        yticks=(0.5:0.5:2, string.(((0.5:0.5:2)))),
        fg_color_legend=:transparent, bg_color_legend=:transparent, 
        legendfontsize=12, legend=:topleft,
        labelfontsize=12, tickfontsize=10,
        ylabel="Normalized to This Study",
        grid=false,
    )

    # Plot data for this study and ±10%
    x = 1:length(anhydrous_majors)
    Plots.plot!(h_base, 
        [xlims(h_base)[1], xlims(h_base)[2], xlims(h_base)[2], xlims(h_base)[1]], 
        [0.9, 0.9, 1.1, 1.1],
        seriestype=:shape, color=:grey, alpha=0.15, lcolor=:match, label=""
    )
    Plots.plot!(h_base, collect(xlims(h_base)), [1,1], label="", color=:black, linewidth=3)

    Plots.plot!(h_base, x, Measurements.value.(ratio.ucc), markershape=:circle, 
        label="Exposed Continental Crust",
        color=:black, msc=:auto, 
        markersize=6, linewidth=2
    )


## --- Continental crust 
    h = deepcopy(h_base)
    Plots.plot!(h, x, Measurements.value.(ratio.plut), markershape=:circle, 
        label="Plutonic Rocks",
        color=colors.plut, msc=:auto, 
        markersize=6, linewidth=2
    )
    Plots.plot!(h, x, Measurements.value.(ratio.ersn), markershape=:circle, 
        label="Eroded Sediment",
        color=:darkgoldenrod, msc=:auto, 
        markersize=6, linewidth=2,
    )

    # Plot previous estimates
    Plots.plot!(h, x, Measurements.value.(ratio.rudnick_gao), markershape=:diamond, 
        label="Rudnick and Gao, 2014",
        msc=colors_source.rudnick,
        linecolor=colors_source.rudnick,
        markercolor=:white, markersize=7, 
        mswidth=2, linewidth=2
    )
    Plots.plot!(h, x, Measurements.value.(ratio.gao), markershape=:diamond, 
        label="Gao et al., 1998",
        msc=colors_source.gao, 
        linecolor=colors_source.gao,
        markercolor=:white, markersize=7, 
        mswidth=2, linewidth=2
    )
    Plots.plot!(h, x, Measurements.value.(ratio.glorise), markershape=:dtriangle, 
        label="Muller et al., 2021",
        color=:cadetblue, msc=:cadetblue, 
        markercolor=:white, markersize=5, 
        mswidth=2, linewidth=2
    )
    display(h)
    savefig(h, "$filepath/modelcomparison.pdf")


## --- End of file 