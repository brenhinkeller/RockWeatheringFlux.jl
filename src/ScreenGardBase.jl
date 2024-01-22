## --- Load data files
    # Main file
    gard_sample = importdataset("data/gard2019/sample.csv", ',', importas=:Tuple)
    npoints = length(gard_sample.sample_id)

    # Other files
    gard_age = importdataset("data/gard2019/age.csv", ',', importas=:Tuple)
    gard_rgroup = importdataset("data/gard2019/rockgroup.csv", ',', importas=:Tuple)
    gard_major = importdataset("data/gard2019/major.csv", ',', importas=:Tuple)
    gard_trace = importdataset("data/gard2019/trace.csv", ',', importas=:Tuple)
    gard_ref = importdataset("data/gard2019/reference.csv", ',', importas=:Tuple)


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
    gard_age.age[(0 .> gard_age.age) .| (gard_age.age .> 4000)] .= NaN


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
    # Also, replace doesn't work so I had to make an extra special function :) <3
    for i in eachindex(gard_sample.rock_name)
        gard_sample.rock_name[i] = replace_malformed_char(gard_sample.rock_name[i])
        gard_sample.sample_description[i] = replace_malformed_char(gard_sample.sample_description[i])
        gard_sample.qap_name[i] = replace_malformed_char(gard_sample.qap_name[i])
    end

    # Match samples
    gard_cats = match_rocktype(
        gard_sample.rock_name, gard_sample.sample_description, gard_sample.qap_name, 
        Int.(gard_sample.rgroup_id), 
        rockgroup_id = Int.(gard_rgroup.rgroup_id), rockgroup_name = vec(rockgroup_name)
    )

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
    majors, minors = get_elements();
    allelements = [majors; minors]
    deleteat!(allelements, findall(x->x==:Volatiles,allelements))    

    # Preallocate
    out = NamedTuple{Tuple(allelements)}(Array{Float64}(undef, npoints) for _ in allelements);

    # # I want a list of all the elements I'm not copying over
    # allkeys = vcat(collect(keys(gard_major)), collect(keys(gard_trace)))
    # allkeys = NamedTuple{Tuple(allkeys)}([false] for _ in allkeys)

    # For each element, lowercase it, find it, and assign to the new Tuple
    p = Progress(length(allelements), desc="Creating geochemial data array", enabled=show_progress)
    for e in allelements
        if haskey(gard_major, Symbol(lowercase(string(e))))
        # Major elements can be copied as-is
            for i in eachindex(major_id)
                out[e][i] = gard_major[Symbol(lowercase(string(e)))][major_id[i]]
            end
            # allkeys[Symbol(lowercase(string(e)))][1] = true

        elseif haskey(gard_trace, Symbol(lowercase(string(e) * "_ppm")))
        # Trace elements should be converted to wt.%
            for i in eachindex(trace_id)
                out[e][i] = gard_trace[Symbol(lowercase(string(e) * "_ppm"))][trace_id[i]] / 10000
            end
            # allkeys[Symbol(lowercase(string(e) * "_ppm"))][1] = true

        elseif e == :FeOT
        # Manually copy FeOT
            for i in eachindex(major_id)
                out.FeOT[i] = gard_major.feo_tot[major_id[i]]
            end
            # allkeys.feo_tot[1] = true
        else
            @warn "No key found for $e"
        end
        next!(p)
    end

    # # Print the list of everything that wasn't copied
    # uncopied = ""
    # for k in keys(allkeys)
    #     if allkeys[k][1] == false
    #         uncopied *= " $k,"
    #     end
    # end
    # println(uncopied)
    # major_id, fe2o3, fe2o3_tot, feo, sro, h2o_plus, h2o_minus, h2o_tot, co2, so3, bao, 
    # caco3, mgco3, loi, trace_id, br_ppm, h_ppm, n_ppm, p_ppm, al_ppm, ca_ppm, cr_ppm, 
    # fe_ppm, ge_ppm, k_ppm, mg_ppm, mn_ppm, na_ppm, ni_ppm, pa_ppm, pm_ppm, rh_ppm, ru_ppm, 
    # si_ppm, ti_ppm,
    
    # Of these, the non-volatiles I maybe want to think about are:
    # Ge, Pa, Pm, Rh, Ru


## --- Check that element / element oxide pairs record the same data
    # Preallocate temporary comparison arrays 
    oxide = Array{Float64}(undef, npoints)
    elem = Array{Float64}(undef, npoints)

    # Check that oxides convert to their elements. Conversion is [1] to [2]
    pairs = [(:SiO2, :si_ppm, (28.085 + 15.999*2)/28.085), 
        (:Al2O3, :al_ppm, (26.98*2 + 15.999*3)/26.98), 
        (:TiO2, :ti_ppm, (47.867 + 15.999*2)/47.867), 
        (:MgO, :mg_ppm, (24.305 + 15.999)/24.305), 
        (:CaO, :ca_ppm, (40.078 + 15.999)/40.078), 
        (:Na2O, :na_ppm, (23.000*2 + 15.999)/23.000), 
        (:K2O, :k_ppm, (39.098*2 + 15.999)/39.098), 
        (:Cr2O3, :cr_ppm, (51.996*2 + 15.999*3)/51.996), 
        (:MnO, :mn_ppm, (54.938 + 15.999)/54.938), 
        (:NiO, :ni_ppm, (58.693 + 15.999)/58.693), 
        (:P2O5, :p_ppm, (30.974*2 + 15.999*5)/30.974),
        (:Ba, :bao, 137.327/(137.327 + 15.999)),
        (:Sr, :sro, 87.62/(87.62 + 15.999))
    ];
    
    for p in pairs
        # Get data into the comparison arrays, in the same units
        oxide .= out[p[1]]
        if p[1] == :Ba || p[1] == :Sr
            for i in eachindex(elem)
                elem[i] = gard_major[p[2]][major_id[i]] / 10000 * p[3]
            end
        else
            for i in eachindex(elem)
                elem[i] = gard_trace[p[2]][trace_id[i]] / 10000 * p[3]
            end
        end

        # Check that overlapping element data is not within 1% of the oxide
        overlap = @. !isnan(oxide) && !isnan(elem);
        overlap_value = falses(length(overlap));
        for i in eachindex(oxide)
            if overlap[i] && (oxide[i]-(oxide[i]*0.01) < elem[i] < oxide[i]+(oxide[i]*0.01))
                overlap_value[i] = true
            end
        end

        # Terminal output
        double_vals = count(overlap)
        duplicates = count(overlap_value)
        pct = round(duplicates/double_vals*100, digits=2)
        # println("Element | Duplicate | Total Overlapping (% Duplicate)")
        println("$(rpad(p[1], 6)) $(rpad(duplicates, 8)) of $(rpad(double_vals, 6)) ($pct%)")

        # Add the element-as-element-oxide values to the element oxide, unless the values
        # are negative, which is stupid
    end
    
    # Treat iron species separately. Make sure Fe₂O₃ + FeO = FeO_T = Fe₂O₃_T
    Fe2O3_to_FeO = (55.845 + 15.999)/(55.845*2 + 15.999*3)

    Fe2O3_as_FeO = similar(out.FeOT)
    Fe2O3_T_as_FeO = similar(out.FeOT)
    FeO = similar(out.FeOT)
    expected_FeO_T = similar(out.FeOT)

    for i in eachindex(out.FeOT)
        # Get the data out of the major file
        Fe2O3_as_FeO[i] = gard_major.fe2o3[major_id[i]] * Fe2O3_to_FeO
        Fe2O3_T_as_FeO[i] = gard_major.fe2o3_tot[major_id[i]] * Fe2O3_to_FeO
        FeO[i] = gard_major.feo[major_id[i]]

        # Calculate expected FeO_T from FeO and Fe₂O₃
        expected_FeO_T[i] = nanadd(Fe2O3_as_FeO[i], FeO[i])
    end

    overlap = @. !isnan(out.FeOT) && !isnan(expected_FeO_T);

    # TO DO: I assume I want to do something with the overlap array??
    # uhhhhhh numbers??
    [out.FeOT[overlap] expected_FeO_T[overlap]]

    overlap .&= .!isnan.(Fe2O3_T_as_FeO)
    [out.FeOT[overlap] expected_FeO_T[overlap] Fe2O3_T_as_FeO[overlap] FeO[overlap]]    


## --- Get volatiles and carbonates as temporary non-Tuple arrays
    # Preallocate
    H2O = similar(out.SiO2)
    CO2 = similar(out.SiO2)
    SO3 = similar(out.SiO2)
    LOI = similar(out.SiO2)
    CaCO3 = similar(out.SiO2)
    MgCO3 = similar(out.SiO2)

    for i in eachindex(H2O)
        H2O[i] = gard_major.h2o_tot[major_id[i]]
        CO2[i] = gard_major.co2[major_id[i]]
        SO3[i] = gard_major.so3[major_id[i]]
        LOI[i] = gard_major.loi[major_id[i]]
        CaCO3[i] = gard_major.caco3[major_id[i]]
        MgCO3[i] = gard_major.mgco3[major_id[i]]
    end


## --- Screen outliers (e.g. remove negative data and physically impossible data)
    # TO DO: where do these lower bounds come from?
    # Screen volatile outliers
    LOI[(0.005 .>= LOI) .| (LOI .>= 22)] .= NaN    # TO DO: this feels very low
    H2O[(0.005 .>= H2O) .| (H2O .>= 30)] .= NaN    
    CO2[(0.005 .>= CO2) .| (CO2 .>= 53)] .= NaN    # Max is pure MgCO₃
    SO3[(0.005 .>= SO3) .| (SO3 .>= 60)] .= NaN    # Max is pure CaSO₄ (32.065 + 15.999*3) / (40.078 + 32.065 + 15.999*4)

    # Screen carbonate outliers 
    CaCO3[(0 .>= CaCO3) .| (CaCO3 .>= 100)] .= NaN
    MgCO3[(0 .>= MgCO3) .| (MgCO3 .>= 100)] .= NaN

    # annnnd... everything else
    screen_outliers!(out)


## --- Calculate volatiles
    # Define some conversion factors 
    CaCO3_to_CaO = (40.078+15.999)/(40.078+15.999*3+12.01)
    CaCO3_to_CO2 = (15.999*2+12.01)/(40.078+15.999*3+12.01)
    MgCO3_to_MgO = (24.305+15.999)/(24.305+15.999*3+12.01)
    MgCO3_to_CO2 = (15.999*2+12.01)/(24.305+15.999*3+12.01)

    CaO_to_SO3 = (32.065+3*15.999)/(40.078+15.999)            # Gypsum / Anhydrite CaSO₄
    CaO_to_H2O = 2*(2*1.00784+15.999)/(40.078+15.999)         # Gypsum CaSO₄ ⋅ 2H₂O

    CaO_to_CO2 = (15.999*2+12.01)/(40.078+15.999)             # Assumes CaCO₃ decomposition

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

    # [CLASTICS, CARBONATES] Compute CO2 from CaO, assuming CaO is from CaCO3 
    # Only replace values if there isn't already a value there, and we didn't convert already
    target = gard_cats.siliciclast .| gard_cats.carb .| gard_cats.shale;
    @turbo for i in eachindex(CO2)
        CO2[i] = ifelse(target[i] & !converted[i] & isnan(CO2[i]), out.CaO[i]*CaO_to_CO2, CO2[i])
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


## --- Compute total wt.% analyzed for all samples
    # Exclude samples below the continental shelf break (-140 m)
    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, gard_sample.latitude, gard_sample.longitude)
    abovesea = elev .> -140;

    # Compute bulk analyzed weight for rocks above sea level
    p = Progress(npoints ÷ 10, desc="Calculating wt.% ...", enabled=show_progress)
    bulkweight = Array{Float64}(undef, npoints)
    @inbounds for i in eachindex(bulkweight)
        if abovesea[i]
            bulkweight[i] = nansum([out[j][i] for j in allelements]) + volatiles[i]
        else
            bulkweight[i] = NaN
        end
        i%10==0 && next!(p)
    end


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

    # Restrict data
    t = @. 84 <= bulkweight <= 104; tᵢ = count(t)           # Without assumed volatiles
    t = @. 84 <= bulkweight .+ volatiles_assumed <= 104;    # With assumed volatiles
    
    @turbo volatiles .+= volatiles_assumed
    isgypsum = vec(gard_sample.rock_name .== "gypsum");

    screenvolatiles!(vec(t), volatiles, isevap=gard_cats.evap, isgypsum=isgypsum,
        general=dol, evaporite=bas, gypsum=gyp    
    );

    # Print to terminal
    @info """Saving $(count(t)) samples ($(round(count(t)/length(t)*100, digits=2))%)
    Assuming volatiles increased count from $tᵢ to $(count(t))
    Total increase = $(count(t) - tᵢ)
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
            Age_Min = age_min
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


## --- Restrict metadata 
    # References
    references = Array{String}(undef, length(out.SiO2))
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

    # References 
    write(text, "references", references)

    # Rock class
    a = Array{Int64}(undef, length(gard_kittens[1]), length(gard_kittens))
    for i in eachindex(keys(gard_kittens))
        for j in eachindex(gard_kittens[i])
            a[j,i] = ifelse(gard_kittens[i][j], 1, 0)
        end
    end
    class["bulk_cats"] = a
    class["bulk_cats_head"] = string.(collect(keys(gard_kittens)))
    
    close(fid)


## --- End of file