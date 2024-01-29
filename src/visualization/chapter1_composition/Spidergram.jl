## --- Set up 
    # Spider... gram. Spidergram.
    # Check out the REE patterns on this guy

    # Load data and base packages
    include("Definitions.jl")

    # Get a list of REEs
    REEs = get_REEs()

    # Load data
    @info "Upper crust data: $ucc_out"
    ucc = importdataset(ucc_out, '\t', importas=:Tuple);
    rudnick_gao = importdataset("data/rudnick_gao_2014_table1-2.csv",  ',', importas=:Tuple);


## --- Load my data, normalized to 100% anhydrous
    # Elements (anhydrous)
    i = findfirst(x->x=="Volatiles", ucc.element)
    elemkeys = Tuple(Symbol.(ucc.element[1:end .!= i]))

    # Rock types, including bulk (all samples). Get all units as ppm
    class = merge(macro_cats, (bulk=trues(length(macro_cats[1])),))
    ucc = NamedTuple{keys(class)}(
        [NamedTuple{elemkeys}(normalize!(ucc[k][1:end .!= i]).*10_000) for k in keys(class)]
    );


## --- Load undifferentiated EarthChem data
    # Also convert units to ppm
    spider_REEs = REEs[1:end .!= findfirst(x->x==:Pm, REEs)]
    class = merge(bulk_cats, (bulk=trues(length(bulk_cats[1])),))
    bulk_REE = NamedTuple{keys(class)}(
        [NamedTuple{Tuple(spider_REEs)}([nanmean(bulk[k][f]).*10_000 for k in spider_REEs]) for f in class]
    );


## --- Load Rudnick and Gao, 2014; Taylor and McLennan, 1985
    # Convert data to dictionaries
    units = Dict(zip(rudnick_gao.Element, rudnick_gao.Units))

    taylor_mclennan = Dict(zip(rudnick_gao.Element, rudnick_gao.Taylor_and_McLennan_1985))
    rudnick_gao = Dict(zip(rudnick_gao.Element, rudnick_gao.This_Study))

    # Convert units to wt.% and normalize to 100%
    for d in (rudnick_gao, taylor_mclennan)
        for k in keys(rudnick_gao)
            if units[k] == "percent"
                continue
            elseif units[k] == "ppm"
                d[k] = d[k] / 10_000
            elseif units[k] == "ppb"
                d[k] = d[k] / 10_000_000
            end
        end
    end
    rudnick_gao = Dict(zip(keys(rudnick_gao), normalize!(collect(values(rudnick_gao)))))
    taylor_mclennan = Dict(zip(keys(taylor_mclennan), normalize!(collect(values(taylor_mclennan)))))

    # Convert normalized REE values to ppm for spidergram
    for d in (rudnick_gao, taylor_mclennan)
        for k in REEs
            haskey(d, string(k)) && (d[string(k)] *= 10_000)
        end
    end

    
## --- Spider... gram. Spidergram.
    h = spidergram(taylor_mclennan, label="Taylor and McLennan, 1985 / 1995", 
        markershape=:diamond, seriescolor=:lightgrey, legend=:bottomleft,
        title="Matched Samples [$geochem_fid]")
    spidergram!(h, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:grey)

    spidergram!(h, ucc.ign, label="This Study (Igneous)",
        markershape=:utriangle, seriescolor=colors.ign)
    spidergram!(h, ucc.granite, label="This Study (Granite)",
        markershape=:utriangle, seriescolor=colors.granite)
    spidergram!(h, ucc.sed, label="This Study (Sedimentary)",
        markershape=:dtriangle, seriescolor=colors.sed)
    spidergram!(h, ucc.shale, label="This Study (Shales)",
        markershape=:dtriangle, seriescolor=colors.shale)

    spidergram!(h, ucc.bulk, label="This Study (Bulk Earth)",
        markershape=:circle, seriescolor=:black)

    display(h)
    savefig("$filepath/spidergram.pdf")


## --- Spidergram comparing bulk geochemical REEs and Rudnick and Gao estimation
    h = spidergram(rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:grey, 
        title="Bulk [$geochem_fid] Geochemical Dataset")

    spidergram!(h, bulk_REE.sed, label="All Sedimentary", 
        markershape=:utriangle, seriescolor=colors.sed)
    spidergram!(h, bulk_REE.shale, label="All Shales", 
        markershape=:utriangle, seriescolor=colors.shale)
    spidergram!(h, bulk_REE.ign, label="All Igneous",
        markershape=:utriangle, seriescolor=colors.ign)

    spidergram!(h, bulk_REE.granite, label="All Granites",
        markershape=:utriangle, seriescolor=colors.granite)
    spidergram!(h, bulk_REE.basalt, label="All Basalts",
        markershape=:utriangle, seriescolor=colors.basalt)
    
    display(h)


## --- Break down igneous rocks a little more
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];

    # Volcanic
    p = Plots.palette(:managua, length(minorvolc))
    h1 = spidergram(bulk_REE.volc, label="Volcanic", markershape=:utriangle,
        seriescolor=colors.volc, legend=:outerright, size=(800, 400),
        title="Volcanic (All Geochemical Samples) [$geochem_fid]")
    h2 = spidergram(ucc.volc, label="Volcanic", markershape=:utriangle,
        seriescolor=colors.volc, legend=:outerright, size=(800, 400),
        title="Volcanic (Matched Geochemical Samples) [$geochem_fid]")
    for i in eachindex(minorvolc)
        spidergram!(h1, bulk_REE[minorvolc[i]], label="$(minorvolc[i])", markershape=:utriangle, color=p[i])
        spidergram!(h2, ucc[minorvolc[i]], label="$(minorvolc[i])", markershape=:utriangle, color=p[i])
    end
    spidergram!(h1, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:black)
    spidergram!(h2, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:black)
    display(h1)
    display(h2)


    # Plutonic
    p = Plots.palette(:managua, length(minorplut))
    h3 = spidergram(bulk_REE.plut, label="Plutonic", markershape=:utriangle,
        seriescolor=colors.plut, legend=:outerright, size=(800, 400),
        title="Plutonic (All Geochemical Samples) [$geochem_fid]")
    h4 = spidergram(ucc.plut, label="Plutonic", markershape=:utriangle,
        seriescolor=colors.plut, legend=:outerright, size=(800, 400),
        title="Plutonic (Matched Geochemical Samples) [$geochem_fid]")
    for i in eachindex(minorplut)
        spidergram!(h3, bulk_REE[minorplut[i]], label="$(minorplut[i])", markershape=:utriangle, color=p[i])
        spidergram!(h4, ucc[minorplut[i]], label="$(minorplut[i])", markershape=:utriangle, color=p[i])
    end
    spidergram!(h3, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:black)
    spidergram!(h4, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:black)
    display(h3)
    display(h4)


    # Bulk igneous
    p = Plots.palette(:managua, length(minorign))
    h5 = spidergram(bulk_REE.ign, label="Igneous", markershape=:utriangle,
        seriescolor=colors.ign, legend=:outerright, size=(800, 400),
        title="Igneous (All Geochemical Samples) [$geochem_fid]")
    h6 = spidergram(ucc.ign, label="Igneous", markershape=:utriangle,
        seriescolor=colors.ign, legend=:outerright, size=(800, 400),
        title="Igneous (Matched Geochemical Samples) [$geochem_fid]")
    for i in eachindex(minorign)
        spidergram!(h5, bulk_REE[minorign[i]], label="$(minorign[i])", markershape=:utriangle, color=p[i])
        spidergram!(h6, ucc[minorign[i]], label="$(minorign[i])", markershape=:utriangle, color=p[i])
    end
    spidergram!(h5, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:black)
    spidergram!(h6, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:black)
    display(h5)
    display(h6)


## --- End of file 