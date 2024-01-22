## --- Set up 
    # REE patterns in my data compared to Rudnick and Gao, 2014 (10.1016/B978-0-08-095975-7.00301-6)
    # REE patterns in the bulk EarthChem data

    # Load data and base packages
    include("Definitions.jl")

    # Get a list of REEs
    REEs = get_REEs()

    # Load data
    @info "Upper crust data: $ucc_out"
    ucc = importdataset(ucc_out, '\t', importas=:Tuple);
    rudnick_gao = importdataset("data/rudnick_gao_2014_table1-2.csv",  ',', importas=:Tuple);


## --- Load my data
    # Convert data to dictionaries
    ucc_ign = Dict(zip(ucc.element, ucc.ign))
    ucc_sed = Dict(zip(ucc.element, ucc.sed))
    ucc = Dict(zip(ucc.element, ucc.bulk))

    # Normalize to 100% anhydrous
    for d in (ucc, ucc_ign, ucc_sed)
        delete!(d, "Volatiles")
    end
    ucc = Dict(zip(keys(ucc), normalize!(collect(values(ucc)))))
    ucc_ign = Dict(zip(keys(ucc_ign), normalize!(collect(values(ucc_ign)))))
    ucc_sed = Dict(zip(keys(ucc_sed), normalize!(collect(values(ucc_sed)))))

    # Convert normalized REE values to ppm for spidergram
    for d in (ucc, ucc_ign, ucc_sed)
        for k in REEs
            haskey(d, string(k)) && (d[string(k)] *= 10_000)
        end
    end


## --- Load undifferentiated EarthChem data
    # Load data into dictionaries and convert to ppm
    absent = findfirst(x->x==:Pm, REEs)
    spider_REEs = REEs[1:end .!= absent]
    bulk_earth = Dict(zip(spider_REEs, 
        [nanmean(bulk[k]) for k in spider_REEs] .*= 10_000)
    )
    bulk_sed = Dict(zip(spider_REEs, 
        [nanmean(bulk[k][bulk_cats.sed]) for k in spider_REEs] .*= 10_000)
    )
    bulk_ign = Dict(zip(spider_REEs, 
        [nanmean(bulk[k][bulk_cats.ign]) for k in spider_REEs] .*= 10_000)
    )
    bulk_granite = Dict(zip(spider_REEs, 
        [nanmean(bulk[k][bulk_cats.granite]) for k in spider_REEs] .*= 10_000)
    )
    bulk_basalt = Dict(zip(spider_REEs, 
        [nanmean(bulk[k][bulk_cats.basalt]) for k in spider_REEs] .*= 10_000)
    )


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
        markershape=:diamond, seriescolor=:lightgrey)
    spidergram!(h, rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:grey)

    spidergram!(h, ucc_ign, label="This Study (Bulk Igneous)",
        markershape=:utriangle, seriescolor=colors.ign)
    spidergram!(h, ucc_sed, label="This Study (Bulk Sedimentary)",
        markershape=:dtriangle, seriescolor=colors.sed)

    spidergram!(h, ucc, label="This Study (Bulk Earth)",
        markershape=:circle, seriescolor=:black)

    display(h)
    savefig("$filepath/spidergram.pdf")


## --- Spidergram comparing bulk EarthChem REEs and Rudnick and Gao estimation
    h = spidergram(rudnick_gao, label="Rudnick and Gao, 2014", 
        markershape=:diamond, seriescolor=:grey)

    spidergram!(h, bulk_sed, label="Bulk Sedimentary", 
        markershape=:utriangle, seriescolor=colors.sed)
    spidergram!(h, bulk_ign, label="Bulk Igneous",
        markershape=:utriangle, seriescolor=colors.ign)

    spidergram!(h, bulk_granite, label="Bulk Granite",
        markershape=:utriangle, seriescolor=colors.granite)
    spidergram!(h, bulk_basalt, label="Bulk Basalt",
        markershape=:utriangle, seriescolor=colors.basalt)

    # spidergram!(h, bulk_earth, label="Bulk EarthChem",
    #     markershape=:utriangle, seriescolor=:brown)

    # spidergram!(h, ucc_ign, label="This Study (Bulk Igneous)",
    #     markershape=:circle, seriescolor=colors.ign)
    # spidergram!(h, ucc_sed, label="This Study (Bulk Sedimentary)",
    #     markershape=:circle, seriescolor=colors.sed)
    # spidergram!(h, ucc, label="This Study (Bulk Earth)",
    #     markershape=:circle, seriescolor=:brown)


## --- End of file 