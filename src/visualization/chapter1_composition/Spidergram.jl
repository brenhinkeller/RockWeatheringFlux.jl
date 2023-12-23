## --- Set up 
    # REE patterns in my data compared to Rudnick and Gao, 2014 (10.1016/B978-0-08-095975-7.00301-6)

    # Load data and base packages
    if !@isdefined(filepath)
        include("Definitions.jl")
    end


## --- Load and parse data
    # Load data
    rg = importdataset("data/rudnick_gao_2014_table1-2.csv",  ',', importas=:Tuple)
    ucc = importdataset(ucc_out, '\t', importas=:Tuple)

    units = Dict(zip(rg.Element, rg.Units))                     # Units for Rudnick and Gao

    # Get dictionaried
    ucc_ign = Dict(zip(ucc.element, ucc.ign))
    ucc_sed = Dict(zip(ucc.element, ucc.sed))
    ucc_met = Dict(zip(ucc.element, ucc.met))
    ucc = Dict(zip(ucc.element, ucc.bulk))                      # My estimate (bulk Earth)
    tm = Dict(zip(rg.Element, rg.Taylor_and_McLennan_1985))     # Taylor and McLennan
    rg = Dict(zip(rg.Element, rg.This_Study))                   # Rudnick and Gao

    # Convert units to percent for normalization
    for d in (rg, tm)
        for k in keys(rg)
            if units[k] == "percent"
                continue
            elseif units[k] == "ppm"
                d[k] = d[k] / 10_000
            elseif units[k] == "ppb"
                d[k] = d[k] / 10_000_000
            end
        end
    end
    rg = Dict(zip(keys(rg), normalize!(collect(values(rg)))))
    tm = Dict(zip(keys(tm), normalize!(collect(values(tm)))))

    # Normalize to 100% anhydrous
    for d in (ucc, ucc_ign, ucc_sed, ucc_met)
        delete!(d, "Volatiles")
    end
    ucc = Dict(zip(keys(ucc), normalize!(collect(values(ucc)))))
    ucc_ign = Dict(zip(keys(ucc_ign), normalize!(collect(values(ucc_ign)))))
    ucc_sed = Dict(zip(keys(ucc_sed), normalize!(collect(values(ucc_sed)))))
    ucc_met = Dict(zip(keys(ucc), normalize!(collect(values(ucc_met)))))

    # We changed everything to wt.%, but spidergram needs ppm
    REEs = get_REEs()
    for d in (rg, tm, ucc, ucc_ign, ucc_sed, ucc_met)
        for k in REEs
            haskey(d, string(k)) && (d[string(k)] *= 10_000)
        end
    end

    
## --- Spider... gram. Spidergram.
    h = spidergram(tm, label="Taylor and McLennan, 1985 / 1995", 
        markershape=:dtriangle, seriescolor=:olivedrab)
    spidergram!(h, rg, label="Rudnick and Gao, 2014", 
        markershape=:utriangle, seriescolor=:cadetblue)

    spidergram!(h, ucc_ign, label="This Study (Bulk Igneous)",
        markershape=:star5, seriescolor=colors.ign)
    spidergram!(h, ucc_sed, label="This Study (Bulk Sedimentary)",
        markershape=:+, seriescolor=colors.sed)
    spidergram!(h, ucc_met, label="This Study (Bulk Metamorphic)",
        markershape=:x, seriescolor=colors.met)

    spidergram!(h, ucc, label="This Study (Bulk Earth)",
        markershape=:circle, seriescolor=:black)

    display(h)
    savefig("$filepath/spidergram.pdf")


## --- End of file 