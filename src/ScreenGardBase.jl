## --- Load data files
    # Main file
    gard_sample = importdataset("data/gard2019/sample.csv", ',', importas=:Tuple);
    npoints = length(gard_sample.sample_id);

    # Other files
    gard_age = importdataset("data/gard2019/age.csv", ',', importas=:Tuple);
    gard_rgroup = importdataset("data/gard2019/rockgroup.csv", ',', importas=:Tuple);
    gard_major = importdataset("data/gard2019/major.csv", ',', importas=:Tuple);
    gard_trace = importdataset("data/gard2019/trace.csv", ',', importas=:Tuple);
    gard_ref = importdataset("data/gard2019/reference.csv", ',', importas=:Tuple);


## --- Fill in any missing ages using bounds and time periods
    # Note that for lack of a better system, if there is no age and only one bound, the
    # age just becomes that bound.

    # From given ages, set missing ages as the average of the upper and lower bounds
    for i in eachindex(gard_age.age)
        gard_age.age[i] = ifelse(isnan(gard_age.age[i]), 
            nanmean([gard_age.age_max[i], gard_age.age_min[i]]), gard_age.age[i]
        )
    end

    # Get GTS bounds and names
    bounds_GTS = get_GTS_boundaries("data/boundaries_green2022.csv");
    names_GTS = (
        eon = string.(keys(bounds_GTS.eon)),
        era = string.(keys(bounds_GTS.era)),
        period = string.(keys(bounds_GTS.period)),
        epoch = string.(keys(bounds_GTS.epoch)),
        stage = string.(keys(bounds_GTS.stage)),
    )

    # Preallocate
    period_age = Array{Float64}(undef, length(gard_age.age), 1)
    period_age_min = Array{Float64}(undef, length(gard_age.age), 1)
    period_age_max = Array{Float64}(undef, length(gard_age.age), 1)

    # Get lowercased, and figure out which samples we can use
    for k in (:time_period_min, :time_period, :time_period_max)
        gard_age[k] .= lowercase.(gard_age[k])
    end
    
    viable = @. (
        isnan(gard_age.age) && (
            (gard_age.time_period != "") || 
            (gard_age.time_period_min != "") || (gard_age.time_period_max != "")
        )
    );

    # Get ages from names
    for i in eachindex(gard_age.age_id)
        if viable[i]
            # If there is a known min age, use that
            if isnan(gard_age.age_min[i]) && gard_age.time_period_min[i] != ""
                time = assign_GTS_age(gard_age.time_period_min[i], names_GTS, bounds_GTS)
                period_age_min[i] = time[1].val - time[1].err
            else
                period_age_min[i] = gard_age.age_min[i]
            end

            # If there is a known max age, use that
            if isnan(gard_age.age_max[i]) && gard_age.time_period_max[i] != ""
                time = assign_GTS_age(gard_age.time_period_max[i], names_GTS, bounds_GTS)
                period_age_max[i] = time[2].val + time[2].err
            else
                period_age_max[i] = gard_age.age_max[i]
            end

            # If we have an age, use that. Otherwise, get the age from the bounds
            if gard_age.time_period[i] != ""
                time = assign_GTS_age(gard_age.time_period[i], names_GTS, bounds_GTS)
                period_age[i] = nanmean(time).val

                # Replace any missing bounds with bounds from the time period
                # Propagate error in the min and max bounds by drawing randomly from the 
                # age distribution (gaussian with mean age, standard deviation error)
                period_age_min[i] = ifelse(isnan(period_age_min[i]), 
                    (randn()+time[1].val)*time[1].err, period_age_min[i]
                )
                period_age_max[i] = ifelse(isnan(period_age_max[i]), 
                    (randn()+time[2].val)*time[2].err, period_age_max[i]
                )

            else
                period_age[i] = nanmean([period_age_min[i], period_age_max[i]])
            end

        else
            period_age[i] = period_age_min[i] = period_age_max[i] = NaN
        end
    end

    # Put the new ages into the age fields
    gard_age.age_min[viable] .= period_age_min[viable]
    gard_age.age[viable] .= period_age[viable]
    gard_age.age_max[viable] .= period_age_max[viable]

    # Remove ages younger than 0 or older than 4000
    gard_age.age[(0 .> gard_age.age) .| (gard_age.age .> 4000)] .= NaN;


## --- Then, fix all the weird and / or misspelled ages (TO DO)



## --- Map ages to samples using the ageID
    # Preallocate
    age = Array{Float64}(undef, length(gard_sample.age_id), 1)
    age_min = Array{Float64}(undef, length(gard_sample.age_id), 1)
    age_max = Array{Float64}(undef, length(gard_sample.age_id), 1)

    # We can use gard_age.age_id as an index... quite neat!
    @assert gard_age.age_id == collect(1:length(gard_age.age_id)) 
    @turbo for i in eachindex(age)
        age[i] = gard_age.age[Int(gard_sample.age_id[i])]
        age_max[i] = gard_age.age_max[Int(gard_sample.age_id[i])]
        age_min[i] = gard_age.age_min[Int(gard_sample.age_id[i])]
    end


## --- Assign rock classes to each sample
    # Get the names of each rock group
    rockgroup_name = Array{String}(undef, length(gard_rgroup.rgroup_id), 1)
    for i in eachindex(rockgroup_name)
        rockgroup_name[i] = strip(join(
            [gard_rgroup.rock_group[i], gard_rgroup.rock_origin[i], gard_rgroup.rock_facies[i]], 
            " "
        ))
    end
    
    # Remove invalid / malformed UTF characters that are in here for god knows what reason
    # Also, replace() doesn't work so I had to make an extra special function :) <3
    for i in eachindex(gard_sample.rock_name)
        gard_sample.rock_name[i] = replace_malformed_char(gard_sample.rock_name[i])
        gard_sample.sample_description[i] = replace_malformed_char(gard_sample.sample_description[i])
        gard_sample.qap_name[i] = replace_malformed_char(gard_sample.qap_name[i])
    end

    # Match samples and tag metamorphic samples
    gard_cats = match_rocktype(
        gard_sample.rock_name, gard_sample.sample_description, gard_sample.qap_name, 
        Int.(gard_sample.rgroup_id), 
        rockgroup_id = Int.(gard_rgroup.rgroup_id), rockgroup_name = vec(rockgroup_name),
        showprogress=show_progress
    );

    metamorph_cats = find_metamorphics(
        gard_sample.rock_name, gard_sample.sample_description, gard_sample.qap_name, 
        Int.(gard_sample.rgroup_id), 
        rockgroup_id = Int.(gard_rgroup.rgroup_id), rockgroup_name = vec(rockgroup_name),
        showprogress=show_progress
    );

    # Just check what we have 
    # not_matched = RockWeatheringFlux.find_unmatched(gard_cats)
    # println("unmatched: $(count(not_matched))")
    # for k in keys(gard_cats)
    #     println("$(rpad(k, 15)) $(lpad(string(count(gard_cats[k])), 6))")
    # end


## --- Standardize element names to match EarthChem data
    # Parse element IDs
    major_id = Int.(gard_sample.major_id)
    trace_id = Int.(gard_sample.trace_id)

    # Get non-volatile elements
    # We're not going to take data for Ge, Pa, Pm, Rh, Ru and that's ok <3
    majors, minors = get_elements();
    allelements = [majors; minors]
    deleteat!(allelements, findall(x->x==:Volatiles,allelements))    

    # Preallocate
    out = NamedTuple{Tuple(allelements)}(Array{Float64}(undef, npoints) for _ in allelements);

    # For each element, lowercase it, find it, and assign to the new Tuple
    p = Progress(length(allelements), desc="Creating geochemial data array", enabled=show_progress)
    for e in allelements
        if haskey(gard_major, Symbol(lowercase(string(e))))
        # Major elements can be copied as-is
            for i in eachindex(major_id)
                out[e][i] = gard_major[Symbol(lowercase(string(e)))][major_id[i]]
            end

        elseif haskey(gard_trace, Symbol(lowercase(string(e) * "_ppm")))
        # Trace elements should be converted to wt.%
            for i in eachindex(trace_id)
                out[e][i] = gard_trace[Symbol(lowercase(string(e) * "_ppm"))][trace_id[i]] / 10000
            end

        elseif e == :FeOT
        # Manually copy FeOT, convert all species to FeOT
            for i in eachindex(major_id)
                FeO = gard_major.feo[major_id[i]]
                FeO_T = gard_major.feo_tot[major_id[i]]
                Fe2O3 = gard_major.fe2o3[major_id[i]]
                Fe2O3_T = gard_major.fe2o3_tot[major_id[i]]

                out.FeOT[i] = feoconversion(FeO, Fe2O3, FeO_T, Fe2O3_T)
            end
        else
            @warn "No key found for $e"
        end
        next!(p)
    end


## --- Get volatiles and carbonates as temporary non-Tuple arrays
    # Preallocate
    H2O = similar(out.SiO2)
    CO2, SO3, LOI = similar(H2O), similar(H2O), similar(H2O)
    CaCO3, MgCO3 = similar(H2O), similar(H2O)

    for i in eachindex(H2O)
        H2O[i] = gard_major.h2o_tot[major_id[i]]
        CO2[i] = gard_major.co2[major_id[i]]
        SO3[i] = gard_major.so3[major_id[i]]
        LOI[i] = gard_major.loi[major_id[i]]
        CaCO3[i] = gard_major.caco3[major_id[i]]
        MgCO3[i] = gard_major.mgco3[major_id[i]]
    end


## --- Do you want to make 77 histograms? (Visualize outliers)
    # If not, this is not the code cell for you. Better luck next time!
    # using Plots, StatsBase
    
    # # Different rock classes need visualized separately.
    # temp_cats = deepcopy(gard_cats);
    # include_minor!(temp_cats);

    # # Major element oxides
    # for m in (:SiO2, :Al2O3, :FeOT, :TiO2, :MgO, :CaO, :Na2O, :K2O)
    #     t = @. 0 < out[m] <= 100;
    #     p = percentile(out[m][t], 99)
    #     psed = percentile(out[m][t .& temp_cats.sed], 99)
    #     pign = percentile(out[m][t .& temp_cats.ign], 99)
    #     pmet = percentile(out[m][t .& temp_cats.met], 99)

    #     xmax = min(100, round(Int, max(percentile(out[m][t], 99.99),
    #         percentile(out[m][t .& temp_cats.sed], 99.99),
    #         percentile(out[m][t .& temp_cats.ign], 99.99),
    #         percentile(out[m][t .& temp_cats.met], 99.99),))
    #     )

    #     h = Plots.histogram(
    #         title="$m; 99th bulk = $(round(p, digits=2)) wt.%",
    #         framestyle=:box,
    #         xlims=(-0.1, xmax),
    #         yticks=false, # I don't care about the actual counts
    #         legend=:topright
    #     )

    #     # Igneous rocks 
    #     c, n = bincounts(out[m][temp_cats.ign], 0, xmax, xmax*2)
    #     n = float(n) ./ nansum(float(n) .* step(c))
    #     Plots.plot!(c, n, seriestype=:bar, color=colors.ign, alpha=0.5, linecolor=:match, label="")
    #     vline!(h, [pign], label="$(round(pign, digits=2))", color=colors.ign, linewidth=2)

    #     # Metamorphic rocks
    #     c, n = bincounts(out[m][temp_cats.met], 0, xmax, xmax*2)
    #     n = float(n) ./ nansum(float(n) .* step(c))
    #     Plots.plot!(c, n, seriestype=:bar, color=colors.met, alpha=0.5, linecolor=:match, label="")
    #     vline!([pmet], label="$(round(pmet, digits=2))", color=colors.met, linewidth=2)

    #     # Sedimentary rocks
    #     c, n = bincounts(out[m][temp_cats.sed], 0, xmax, xmax*2)
    #     n = float(n) ./ nansum(float(n) .* step(c))
    #     Plots.plot!(c, n, seriestype=:bar, color=colors.sed, alpha=0.5, linecolor=:match, label="")
    #     vline!([psed], label="$(round(psed, digits=2))", color=colors.sed, linewidth=2)

    #     vline!([p], color=:black, label="$(round(p, digits=2))", linewidth=1)
    #     display(h)
    # end

    # # Double check sedimentary iron... it's chert time
    # minorsed = get_rock_class()[2];
    # h = Plots.histogram(
    #     framestyle=:box,
    #     xlims=(-0.1, 100),
    #     yticks=false, # I don't care about the actual counts
    #     legend=:topright
    # )
    # for k in minorsed
    #     c, n = bincounts(out.FeOT[temp_cats[k]], 0, 100, 200)
    #     n = float(n) ./ nansum(float(n) .* step(c))
    #     Plots.plot!(h, c, n, colors=colors[k], seriestype=:bar, alpha=0.5, linecolor=:match, label="$k")
    # end
    # display(h)

    # # Minor elements and element oxides
    # for m in minors
    #     # Create a ppm dataset because functions only like integers
    #     # Convert all the printouts and ticks back to wt.% though
    #     m_ppm = out[m] .* 10_000;

    #     # Now proceed as normal
    #     t = @. 0 < m_ppm <= 100*10_000;
    #     p = percentile(m_ppm[t], 99)
    #     psed = percentile(m_ppm[t .& temp_cats.sed], 99)
    #     pign = percentile(m_ppm[t .& temp_cats.ign], 99)
    #     pmet = percentile(m_ppm[t .& temp_cats.met], 99)

    #     xmax = min(20*10_000, round(Int, max(percentile(m_ppm[t], 99),
    #         percentile(m_ppm[t .& temp_cats.sed], 99),
    #         percentile(m_ppm[t .& temp_cats.ign], 99),
    #         percentile(m_ppm[t .& temp_cats.met], 99),)*1.1)
    #     )

    #     h = Plots.histogram(
    #         title="$m; 99th bulk = $(round(p/10_000, sigdigits=2)) wt.%",
    #         framestyle=:box,
    #         xlims=(-0.0001, xmax),
    #         xticks=(0:xmax÷5:xmax, string.(round.((0:xmax÷5:xmax) ./ 10_000, digits=4))),
    #         yticks=false, # I don't care about the actual counts
    #         legend=:topright,
    #     )

    #     # Igneous rocks 
    #     c, n = bincounts(m_ppm[temp_cats.ign], 0, xmax, 100)
    #     length(c) != 100 && (c = cntr(range(start=0, stop=xmax, length=101)))  # Sometimes roundoff error is bad
    #     n₁ = float(n) ./ nansum(float(n) .* step(c))
 
    #     # Metamorphic rocks
    #     c, n = bincounts(m_ppm[temp_cats.met], 0, xmax, 100)
    #     n₂ = float(n) ./ nansum(float(n) .* step(c))
    #     length(c) != 100 && (c = cntr(range(start=0, stop=xmax, length=101)))
        
    #     # Sedimentary rocks
    #     c, n = bincounts(m_ppm[temp_cats.sed], 0, xmax, 100)
    #     n₃ = float(n) ./ nansum(float(n) .* step(c))
    #     length(c) != 100 && (c = cntr(range(start=0, stop=xmax, length=101)))

    #     # # Kill the first bin if most data is really small
    #     # if maximum(n₁) == n₁[1] || maximum(n₂) == n₂[1] || maximum(n₃) == n₃[1]
    #     #     n₁[1] = n₂[1] = n₃[1] = 0.0
    #     # end

    #     # Plot things
    #     Plots.plot!(c, n₁, seriestype=:bar, color=colors.ign, alpha=0.5, linecolor=:match, label="")
    #     vline!(h, [pign], label="$(round(pign/10_000, sigdigits=2))", color=colors.ign, linewidth=2)

    #     Plots.plot!(c, n₂, seriestype=:bar, color=colors.met, alpha=0.5, linecolor=:match, label="")
    #     vline!([pmet], label="$(round(pmet/10_000, sigdigits=2))", color=colors.met, linewidth=2)

    #     Plots.plot!(c, n₃, seriestype=:bar, color=colors.sed, alpha=0.5, linecolor=:match, label="")
    #     vline!([psed], label="$(round(psed/10_000, sigdigits=2))", color=colors.sed, linewidth=2)

    #     vline!([p], color=:black, label="$(round(p/10_000, sigdigits=2))", linewidth=1)

    #     display(h)
    #     # println("$m=(0,$(round(p, sigdigits=3))),")
    # end

    # # Volatiles and carbonates 
    # temp = (H2O=H2O, CO2=CO2, SO3=SO3, LOI=LOI, CaCO3=CaCO3, MgCO3=MgCO3);
    # for m in keys(temp)
    #     t = @. 0 < temp[m] <= 100;
    #     p = percentile(temp[m][t], 99)
    #     h = Plots.histogram(temp[m][t], 
    #         title="$m; 99th pct = $(round(p, digits=2)) wt.%",
    #         framestyle=:box, label="", linecolor=:match,
    #         xlims=(-1, min(100, maximum(temp[m][t]) + maximum(temp[m][t])*0.01)),
    #         yticks=false, # I don't care about the actual counts
    #     )
    #     vline!([p], label="")
    #     display(h)
    # end


## --- Screen outliers (e.g. remove negative data and physically impossible data)
    # Screen volatile outliers
    LOI[.!(0 .< LOI .< 50)] .= NaN    # Nothing normal is above 50 (looked at histogram)
    H2O[.!(0 .< H2O .< 50)] .= NaN    # 99th percentile is 71 but I don't trust like that
    CO2[.!(0 .< CO2 .< 53)] .= NaN    # Max is pure MgCO₃; 99th percentile is 44 
    SO3[.!(0 .< SO3 .< 60)] .= NaN    # Max is pure CaSO₄; 99th percentile is 63

    # Screen carbonate outliers 
    CaCO3[.!(0 .< CaCO3 .< 100)] .= NaN
    MgCO3[.!(0 .< MgCO3 .< 100)] .= NaN     # 99th percentile is 41%

    # annnnd... everything else
    screen_outliers!(out, gard_cats);


## --- Calculate volatiles
    # Define some conversion factors 
    CaCO3_to_CaO = (40.078+15.999)/(40.078+15.999*3+12.01)
    CaCO3_to_CO2 = (15.999*2+12.01)/(40.078+15.999*3+12.01)
    MgCO3_to_MgO = (24.305+15.999)/(24.305+15.999*3+12.01)
    MgCO3_to_CO2 = (15.999*2+12.01)/(24.305+15.999*3+12.01)

    CaO_to_SO3 = (32.065+3*15.999)/(40.078+15.999)            # Gypsum / Anhydrite CaSO₄
    CaO_to_H2O = 2*(2*1.00784+15.999)/(40.078+15.999)         # Gypsum CaSO₄ ⋅ 2H₂O

    CaO_to_CO2 = (15.999*2+12.01)/(40.078+15.999)             # Assumes CaCO₃ decomposition

    K2O_to_Al2O3 = (26.98*2+15.999*3)/(39.098*2+15.999)       # K-spar
    Na2O_to_Al2O3 = (26.98*2+15.999*3)/(22.990*2+15.999)      # Na-plag
    Al2O3_to_CaO = (40.078+15.999)/(26.98*2+15.999*3)         # Ca-plag

    # Preallocate
    CO2_CaCO3, CaO_CaCO3 = similar(out.SiO2), similar(out.SiO2)
    CO2_MgCO3, MgO_MgCO3 = similar(out.SiO2), similar(out.SiO2)
    volatiles = similar(out.SiO2)
    
    # Compute CaO / MgO and CO2 from CaCO3 and MgCO3
    # Only replace values if there isn't already a value there
    converted = falses(npoints)
    @turbo for i in eachindex(CO2_CaCO3)
        out.CaO[i] = ifelse(isnan(out.CaO[i]), CaCO3[i] * CaCO3_to_CaO, out.CaO[i])
        out.MgO[i] = ifelse(isnan(out.MgO[i]), MgCO3[i] * MgCO3_to_MgO, out.MgO[i])
        CO2[i] = ifelse(isnan(CO2[i]), (CaCO3[i]*CaCO3_to_CO2 + MgCO3[i]*MgCO3_to_CO2), CO2[i])

        converted[i] = (CaCO3[i] * CaCO3_to_CaO) > 0        
    end

    # [CARBONATES] Compute CO2 from CaO, assuming CaO is from CaCO3 
    # Only replace values if there isn't already a value there, and we didn't convert already
    target = @. gard_cats.carb & !converted & isnan(CO2)
    @turbo for i in eachindex(CO2)
        CO2[i] = ifelse(target[i], out.CaO[i]*CaO_to_CO2, CO2[i])
    end

    # [CLASTICS] Compute CO2 from excess CaO
    # Some Ca is present in feldspars; only convert if there's more CaO than one would 
    # expect from Ca plagioclase
    CaO_excess = similar(CO2)
    @turbo for i in eachindex(CaO_excess)
        CaO_plag = ca_plagioclase(out.K2O[i], out.Al2O3[i], out.Na2O[i], 
            K2O_to_Al2O3, Na2O_to_Al2O3, Al2O3_to_CaO
        )
        CaO_excess[i] = out.CaO[i] - CaO_plag
    end
    
    target = @. (gard_cats.siliciclast .| gard_cats.shale) & !converted & isnan(CO2)
    @turbo for i in eachindex(CO2)
        CO2[i] = ifelse(target[i] & (CaO_excess[i] > 0), CaO_excess[i]*CaO_to_CO2, CO2[i])
    end

    # [EVAPORITES] Compute H2O [GYPSUM] and SO3 [GYPSUM, ANHYDRITE] from CaO
    for i in eachindex(SO3)
        if gard_cats.evap[i] && gard_sample.rock_name[i] == "gypsum"
            SO3[i] = ifelse(isnan(SO3[i]), out.CaO[i]*CaO_to_SO3, SO3[i])
            H2O[i] = ifelse(isnan(H2O[i]), out.CaO[i]*CaO_to_H2O, H2O[i])
        elseif gard_cats.evap[i] && gard_sample.rock_name[i] == "anhydrite"
            SO3[i] = ifelse(isnan(SO3[i]), out.CaO[i]*CaO_to_SO3, SO3[i])
        end
    end

    # Calculate reported volatiles as the larger of (CO₂ + H₂O + SO₃) or LOI
    @turbo volatiles .= vec(nansum([CO2 H2O SO3], dims=2))
    @turbo for i in eachindex(volatiles)
        volatiles[i] = ifelse(LOI[i] > volatiles[i], LOI[i], volatiles[i])
    end


## --- Check that element / element oxide pairs record the same data
    # C and CO2: if within 1%, assume it's the same data and delete C value
    C_as_CO2 = out.C .* (molarmass["C"] + molarmass["O"]*2)/molarmass["C"];
    for i in eachindex(CO2)
        # Delete C value if it's within 1% of CO2 value
        if (CO2[i]-(CO2[i]*0.01) < C_as_CO2[i] < CO2[i]+(CO2[i]*0.01))
            out.C[i] = NaN
        end
    end

    # S and SO3, same thing 
    S_as_SO3 = out.S .* (molarmass["S"] + molarmass["O"]*3)/molarmass["S"];
    for i in eachindex(SO3)
        # Delete S value if it's within 1% of SO3 value
        if (SO3[i]-(SO3[i]*0.01) < S_as_SO3[i] < SO3[i]+(SO3[i]*0.01))
            out.S[i] = NaN
        end
    end

    # Convert BaO and SrO [wt.%] to Ba and Sr 
    Ba = gard_major.bao[major_id] .* (molarmass["Ba"]/(molarmass["Ba"] + molarmass["O"]));
    Sr = gard_major.sro[major_id] .* (molarmass["Sr"]/(molarmass["Sr"] + molarmass["O"]));
    for i in eachindex(out.Ba)
        # Add Ba from BaO to trace Ba, as long as values aren't within 1% of each other
        if !(out.Ba[i]-(out.Ba[i]*0.01) < Ba[i] < out.Ba[i]+(out.Ba[i]*0.01))
            out.Ba[i] = nanadd(Ba[i], out.Ba[i])
        end

        # Same deal with Sr and SrO 
        if !(out.Sr[i]-(out.Sr[i]*0.01) < Sr[i] < out.Sr[i]+(out.Sr[i]*0.01))
            out.Sr[i] = nanadd(Sr[i], out.Sr[i])
        end
    end

    # Actually, don't do this, because we lose samples this way.
    # This might be because some of the data wasn't actually ppm? It's not super clear
    # # Convert major elements to major element oxides
    # oxides = (:SiO2, :Al2O3, :TiO2, :MgO, :CaO, :Na2O, :K2O, :Cr2O3, :MnO, :NiO, :P2O5);
    # sourcekey = (:si_ppm, :al_ppm, :ti_ppm, :mg_ppm, :ca_ppm, :na_ppm, :k_ppm, :cr_ppm, 
    #     :mn_ppm, :ni_ppm, :p_ppm);
    # conversionfactor = (
    #     (molarmass["Si"]   + molarmass["O"]*2)/molarmass["Si"],
    #     (molarmass["Al"]*2 + molarmass["O"]*3)/molarmass["Al"],
    #     (molarmass["Ti"]   + molarmass["O"]*2)/molarmass["Ti"],
    #     (molarmass["Mg"]   + molarmass["O"]  )/molarmass["Mg"],
    #     (molarmass["Ca"]   + molarmass["O"]  )/molarmass["Ca"],
    #     (molarmass["Na"]*2 + molarmass["O"]  )/molarmass["Na"],
    #     (molarmass["K"] *2 + molarmass["O"]  )/molarmass["K"],
    #     (molarmass["Cr"]*2 + molarmass["O"]*3)/molarmass["Cr"],
    #     (molarmass["Mn"]   + molarmass["O"]  )/molarmass["Mn"],
    #     (molarmass["Ni"]   + molarmass["O"]  )/molarmass["Ni"],
    #     (molarmass["P"]*2  + molarmass["O"]*5)/molarmass["P"],
    # );
    # for i in eachindex(oxides)
    #     # Add computed oxide to existing oxide, unless they're within 1% of each other
    #     computed = gard_trace[sourcekey[i]][trace_id] ./ 10000 .* conversionfactor[i]
    #     lower = out[oxides[i]] .- (out[oxides[i]] .* 0.01)
    #     upper = out[oxides[i]] .+ (out[oxides[i]] .* 0.01)

    #     overlap = trues(npoints);
    #     for j in eachindex(computed)
    #         if !(lower[j] < computed[j] < upper[j])
    #             out[oxides[i]][j] = nanadd(out[oxides[i]][j], computed[j])
    #             overlap[j] = false
    #         end
    #     end

    #     total = count(.!isnan.(out[oxides[i]]) .& .!isnan.(computed));
    #     pct = round(count(overlap)/total*100, digits=2)
    #     println("$(rpad(oxides[i], 6)) $(rpad(count(overlap), 8)) of $(rpad(total, 6)) ($pct%)")
    # end


## --- [commented out] Look at C / CO2
    # C_as_CO2 = out.C .* (molarmass["C"] + molarmass["O"]*2)/molarmass["C"];
    # overlap = falses(npoints)
    # for i in eachindex(CO2)
    #     overlap[i] = (CO2[i]-(CO2[i]*0.01) < C_as_CO2[i] < CO2[i]+(CO2[i]*0.01))
    # end

    # # Run total analyzed wt.% code first!
    # # If we subtract the overlapping samples, more samples are in the acceptable wt.% range
    # count(overlap)                                          # 66
    # count(overlap .& (bulkweight .> 104))                   # 3
    # count((bulkweight[overlap] .- CO2[overlap]) .> 104)     # 0
    # count((bulkweight[overlap] .- CO2[overlap]) .< 84)      # 1
    # count(bulkweight[overlap] .< 84)                        # 0

    # # Look at the values we're replacing 
    # [C_as_CO2[overlap] CO2[overlap]]


## --- [commented out] Look at oxide / element.
    # # Run oxide / sourcekey / conversion factor code first!
    # for i in eachindex(oxides)
    #     # Add computed oxide to existing oxide, unless they're within 1% of each other
    #     computed = gard_trace[sourcekey[i]][trace_id] ./ 10000 .* conversionfactor[i]
    #     lower = out[oxides[i]] .- (out[oxides[i]] .* 0.01)
    #     upper = out[oxides[i]] .+ (out[oxides[i]] .* 0.01)

    #     overlap = falses(npoints);
    #     for j in eachindex(computed)
    #         if (lower[j] < computed[j] < upper[j])
    #             overlap[j] = true
    #         end
    #     end
        
    #     # Total number of samples with values in the oxide and element
    #     total = count(.!isnan.(out[oxides[i]]) .& .!isnan.(computed));
    #     pct = round(count(overlap)/total*100, digits=2)
    #     println("$(rpad(oxides[i], 6)) $(rpad(count(overlap), 8)) of $(rpad(total, 6)) ($pct%)")
    # end


## --- Compute total wt.% analyzed for all samples
    # Exclude samples OIBs and samples below the continental shelf break (-140 m)
    isOIB = findOIBs(gard_sample.latitude, gard_sample.longitude);

    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, gard_sample.latitude, gard_sample.longitude)
    abovesea = elev .> -140;

    # Compute bulk analyzed weight for rocks above sea level
    p = Progress(npoints ÷ 10, desc="Calculating wt.% ...", enabled=show_progress)
    bulkweight = Array{Float64}(undef, npoints)
    @inbounds for i in eachindex(bulkweight)
        if abovesea[i] & !isOIB[i]
            bulkweight[i] = nansum([out[j][i] for j in allelements]) + volatiles[i]
        else
            bulkweight[i] = NaN
        end
        i%10==0 && next!(p)
    end


## --- Remove some samples for being ontologically bad 
    # Preallocate a BitVector for the screened samples 
    s = trues(npoints)

    # Remove stream sediment 
    s .&= .!(gard_sample.material .== "stream sediment")

    # Get a list of authors so we can remove some of them    
    author = [gard_ref.author[Int(gard_sample.ref_id[i])] for i in eachindex(gard_sample.ref_id)]

    # These are reported incorrectly (Brenhin / Blair checked the pub)
    s .&= .!(author .== "kylander-clark, andrew r c; coleman, drew s; glazner, allen f; bartley, john m, 2005");

    # These are for zirons, not whole-rock samples - from Brenhin / Blair
    s .&= .!(author .== "lanphere, marvin a; baadsgaard, h, 2001");

    # These are all typos (U>100 ppm) - from Brenhin / Blair
    s .&= .!(author .== "Donskaya, T. V.; Gladkochub, D. P. & Mazukabzov, A. M.");


## --- Correct for assumed unmeasured volatile loss
    # Preallocate
    volatiles_assumed = Array{Float64}(undef, npoints)

    # Consider all sedimentary samples, and only those samples
    minorsed = get_rock_class()[2];
    for type in minorsed
        gard_cats.sed .|= gard_cats[type]
    end

    # If total analyzed wt.% below 100%, assume the unreported data are volatiles
    @turbo for i in eachindex(bulkweight)
        volatiles_assumed[i] = ifelse(gard_cats.sed[i] & (bulkweight[i] < 100), (100 - bulkweight[i]), 0.0)
    end

    # Restrict data to between 84 and 104 wt.%
    t = @. 84 <= bulkweight <= 104; tᵢ = count(t)           # Without assumed volatiles
    t = @. 84 <= bulkweight .+ volatiles_assumed <= 104;    # With assumed volatiles
    t .&= s                                                 # Nothing ontologically bad
    
    @turbo volatiles .+= volatiles_assumed
    isgypsum = vec(gard_sample.rock_name .== "gypsum");

    screenvolatiles!(vec(t), volatiles, isevap=gard_cats.evap, isgypsum=isgypsum,
        general=dol, evaporite=bas, gypsum=gyp    
    );

    # Print to terminal
    up = count(t) - tᵢ
    @info """Saving $(count(t)) samples ($(round(count(t)/length(t)*100, digits=2))%)
    Assuming volatiles increased count from $tᵢ to $(count(t))
    Total increase = $up
    """


## --- Restrict to in-bounds only and normalize composition 
    # This is inefficient, but is not sensitive to changes in the order of any elements
    out = merge(out, (
            Volatiles=volatiles,
            Latitude = gard_sample.latitude,
            Longitude = gard_sample.longitude,
            Loc_Prec = gard_sample.loc_prec,
            Age = age,
            Age_Max = age_max,
            Age_Min = age_min,
            Sample_ID = collect(1:npoints)
        )
    )
    out = NamedTuple{keys(out)}([out[i][t] for i in keys(out)])

    # Normalize to 100%
    p = Progress(length(out.SiO2) ÷ 10, desc="Normalizing compositions ...", enabled=show_progress)
    contributing = [allelements; :Volatiles]             # Need to re-include volatiles!
    for i in eachindex(out.SiO2)
        sample = [out[j][i] for j in contributing]       # Get it
        normalize!(sample)                               # Normalize it
        for j in eachindex(contributing)                 # Put it back
            out[contributing[j]][i] = sample[j]
        end
        i%10==0 && next!(p)
    end

    # While we're here, exclude locations that aren't possible 
    out.Latitude[(-90 .> out.Latitude) .| (out.Latitude .> 90)] .= NaN
    out.Longitude[(-180 .> out.Longitude) .| (out.Longitude .> 180)] .= NaN


 ## --- Get or restrict metadata
    # References
    references = Array{String}(undef, length(gard_sample.ref_id))
    for i in eachindex(references)
        j = Int(gard_sample.ref_id[i])
        year = isnan(gard_ref.year[j]) ? "" : string(Int(gard_ref.year[j]))

        references[i] = join([gard_ref.author[j], gard_ref.title[j], gard_ref.journal[j],
            year, gard_ref.doi[j]], " | "
        )
    end

    # Rock class
    gard_kittens = NamedTuple{keys(gard_cats)}(gard_cats[k][t] for k in keys(gard_cats))
    for k in keys(gard_kittens)
        if show_progress && count(gard_kittens[k]) < 10
            @warn "After filtering, type \"$k\" has $(count(gard_kittens[k])) samples."
        end
    end
    metamorph_kittens = NamedTuple{keys(metamorph_cats)}(metamorph_cats[k][t] for k in keys(metamorph_cats))

     
## --- Save all data to file
    fid = h5open(fileout, "w")
    data = create_group(fid, "bulk")
    text = create_group(fid, "bulktext")
    class = create_group(fid, "bulktypes")

    # Data 
    allkeys = collect(keys(out))
    write(data, "header", string.(allkeys))
    writebulk = create_dataset(data, "data", Float64, (count(t), length(allkeys)))
    for i in eachindex(allkeys)
        writebulk[:,i] = out[allkeys[i]]
    end

    # References and rock descriptions
    write(text, "Reference", references[t])
    write(text, "Rock_Name", gard_sample.rock_name[t])
    write(text, "Sample_Description", gard_sample.sample_description[t])
    write(text, "QAP_Name", gard_sample.qap_name[t])

    # Rock class
    a = Array{Int64}(undef, length(gard_kittens[1]), length(gard_kittens))
    for i in eachindex(keys(gard_kittens))
        for j in eachindex(gard_kittens[i])
            a[j,i] = ifelse(gard_kittens[i][j], 1, 0)
        end
    end

    b = similar(a)
    for i in eachindex(keys(metamorph_kittens))
        for j in eachindex(metamorph_kittens[i])
            b[j,i] = ifelse(metamorph_kittens[i][j], 1, 0)
        end
    end

    class["bulk_cats"] = a
    class["metamorphic_cats"] = b
    class["bulk_cats_head"] = string.(collect(keys(gard_kittens)))
    
    close(fid)


## --- End of file