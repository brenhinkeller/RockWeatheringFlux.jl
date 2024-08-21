## --- Set up 
    # Visualization of the composition of continental crust and eroded material 

    # Load data and base packages
    include("Definitions.jl")

    # List of (major) elements 
    majors, minors = get_elements()
    anhyd_majors = Tuple(majors[1:end .!= findfirst(x->x==:Volatiles, majors)])
    x = 1:length(anhyd_majors)


## --- Load uppermost continental crust [This study, Rudnick and Gao, Gao et al.]
    # This study
    ucc = importdataset(ucc_out, '\t', importas=:Tuple);
    elementindex = NamedTuple{Tuple(Symbol.(ucc.element))}(i for i in eachindex(ucc.element))
    ucc_bulk = NamedTuple{anhyd_majors}([ucc.bulk[elementindex[k]] for k in anhyd_majors]) 
    ucc_plut =NamedTuple{anhyd_majors}([ucc.plut[elementindex[k]] for k in anhyd_majors]) 

    # Rudnick and Gao 
    rudnick_gao = importdataset("data/rudnick_gao_2014_table1-2.csv",  ',', importas=:Tuple);
    elementindex = NamedTuple{Tuple(Symbol.(rudnick_gao.Element))}(
        i for i in eachindex(rudnick_gao.Element))
    rudnick_gao = NamedTuple{Tuple(anhyd_majors)}((haskey(elementindex, k) ? 
        rudnick_gao.This_Study[elementindex[k]] : NaN) for k in anhyd_majors
    );

    # Gao et al. 
    gao = (SiO2=58.48,Al2O3=12.12,FeOT=4.6,TiO2=0.57,MgO=3.73,CaO=7.41,Na2O=2.57,K2O=2.27);


## --- Load eroded material [This study, this study anhydrous, Muller et al.] 
    # This study 
    ersn = importdataset(erodedcomp_out, '\t', importas=:Tuple);
    vᵢ = findfirst(x->x==:Volatiles, majors)
    bulk_volatiles = ersn.bulk[vᵢ]

    elementindex = NamedTuple{Tuple(Symbol.(ersn.elem))}(i for i in eachindex(ersn.elem))
    ersn = NamedTuple{anhyd_majors}([ersn.bulk[elementindex[k]] for k in anhyd_majors])

    # This study, normalized to 100% anhydrous 
    sum_anhydrous = 100 - bulk_volatiles
    ersn_anhydrous = NamedTuple{keys(ersn)}(values(ersn) ./ sum_anhydrous .* 100)

    # Muller et al.
    muller = (SiO2=65.1,Al2O3=18.7,FeOT=5.67,TiO2=0.8,MgO=2.1,CaO=2.9,Na2O=1.1,K2O=2.7,);


## --- Normalize all values to my uppermost continental crust values 
    ratio = (
        ucc = [ucc_bulk[k] ./ ucc_bulk[k] for k in anhyd_majors],
        plut = [ucc_plut[k] ./ ucc_bulk[k] for k in anhyd_majors],
        rudnick_gao = [rudnick_gao[k] ./ ucc_bulk[k] for k in anhyd_majors],
        gao = [gao[k] ./ ucc_bulk[k] for k in anhyd_majors],
        ersn = [ersn[k] ./ ucc_bulk[k] for k in anhyd_majors],
        ersn_anhyd = [ersn_anhydrous[k] ./ ucc_bulk[k] for k in anhyd_majors],
        muller = [muller[k] ./ ucc_bulk[k] for k in anhyd_majors],
    )


## --- Create base plot
    h = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica, 
        xticks=(1:length(anhyd_majors), string.(anhyd_majors)), 
        xlims=(0.5, length(anhyd_majors)+0.5),
        ylims=(0.25, 2.25),
        fg_color_legend=:white, legend=:topleft,
        # ylabel="Normalized to Surface Earth",
        grid=false,
        legendfontsize=12,
        tickfontsize=14,
    );

    Plots.plot!(h,
        [xlims(h)[1], xlims(h)[2], xlims(h)[2], xlims(h)[1]], [0.9, 0.9, 1.1, 1.1],
        seriestype=:shape, color=:grey, alpha=0.15, lcolor=:match, label=""
    )
    Plots.plot!(h, collect(xlims(h)), [1,1], label="", color=:black,)
    Plots.plot!(h, x, ratio.ucc, markershape=:circle, label="Surface Earth (This Study)",
        color=:black, msc=:auto, markersize=5)
    


## --- Plot upper continental crust composition
    h1 = deepcopy(h)
    Plots.plot!(h1, x, ratio.plut, markershape=:circle, 
        label="Plutonic (This Study)",
        color=pal[2], msc=:auto, markersize=5,
    )
    Plots.plot!(h1, x, ratio.gao, markershape=:diamond, 
        label="Gao et al., 1998",
        markercolor=:white, msc=pal[7], color=pal[7], markersize=6,
    )
    Plots.plot!(h1, x, ratio.rudnick_gao, markershape=:diamond, 
        label="Rudnick and Gao, 2014", markerstrokewidth=2,
        markercolor=:white, msc=pal[4], color=pal[4], markersize=6,
    )

    display(h1)
    savefig(h1, "$filepath/composition_ucc.pdf")


## --- Plot eroded material composition
    h2 = deepcopy(h)

    Plots.plot!(h2, x, ratio.ersn, markershape=:circle, 
        label="Eroded Material (This Study)",
        color=pal[2], msc=:auto, markersize=6,
    )
    # Plots.plot!(h2, x, ratio.ersn_anhyd, markershape=:circle, 
    #     label="Anhydrous Eroded Material (This Study)",
    #     color=pal[2], msc=:auto, markersize=6,
    # )    
    Plots.plot!(h2, x, ratio.muller, markershape=:diamond, 
        label="Suspended Sediment (Muller et al., 2021)",
        markercolor=:white, msc=pal[7], color=pal[7], markersize=6,
        markerstrokewidth=2,
    )

    display(h2)
    savefig(h2, "$filepath/composition_ersn.pdf")

    
## --- End of File 