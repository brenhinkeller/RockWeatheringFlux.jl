### --- Load data files (~10 minutes)
    combined = importdataset("data/combined.tsv", '\t', importas=:Tuple)
    npoints = length(combined.SiO2)

    show_progress && @info """Database loaded $(Dates.format(now(), "HH:MM"))"""


## --- Convert rock names to lithologic classes and add metamorphic tags to each sample
    @time combo_cats, meta_cats = match_rocktype(combined.Rock_Group, combined.Rock_Subgroup, 
        combined.Rock_Composition, combined.Rock_Facies, combined.Rock_Name
    ); 
    include_minor!(combo_cats);


## --- Remove improbable ages 
    # We assume missing ages have been filled in from the min / max bounds, if available, 
    # and that ages older than 3800 and younger than 0 have been excluded.

    # Here we exclude ages which are older than the maximum age of the craton as recorded 
    # by the TC1 tectonic age model of Artemieva 2006, doi: 10.1016/j.tecto.2005.11.022
    tc1_age, tc1_min, tc1_max = find_tc1_age(combined.Latitude, combined.Longitude)

    # Get maximum ages. Add 10% to all TC1 maximum and geochemial minimum ages. (Note 
    # absolute minimum TC1 age of 3500 Ma.) If bounds are not known, default to 20% of 
    # recorded age.
    tc1_max *= 1.1
    tc1_max[isnan.(tc1_max)] .= tc1_age[isnan.(tc1_max)] * 1.2
    combined_min = combined.Age_Min * 0.9
    combined_min[isnan.(combined_min)] .= combined.Age[isnan.(combined_min)] * 0.8

    # NaN out any ages that are older than the cratonization age, allowing 10% error tolerance
    invalidage = @. !isnan(tc1_age) && (combined.Age > tc1_age) && (combined_min > tc1_max);
    # difference = abs.(combined.Age - tc1_age)
    # difference[.!invalidage] .= NaN 
    # println(count(invalidage))
    # h = histogram(difference, yaxis=:log10)
    # display(h)

    combined.Age[invalidage .& combo_cats.sed] .= NaN;
    

## --- Convert all units to wt.% and create geochemical data Tuple 
    # Preallocate 
    majors, minors = get_elements()
    isotopes = get_isotopes()
    allelements = [majors; minors; isotopes]
    deleteat!(allelements, findall(x->x==:Volatiles,allelements))
    
    out = NamedTuple{Tuple(allelements)}(Array{Float64}(undef, npoints) for _ in allelements);

    # Define elements not listed in ppm
    wtpct = (majors..., :Cr2O3,:NiO,:P2O5,:MnO);

    # Assign each element to the new Tuple 
    p = Progress(length(allelements), desc="Creating geochemial data array", enabled=show_progress)
    for k in allelements 
        # Convert to wt.% if necessary
        if !(k in wtpct)
            out[k] .= combined[k] ./ 10_000
        else
            out[k] .= combined[k]
        end
        next!(p)
    end


## --- Stoichiometric volatiles calculation 
    # Preallocate 
    CO2 = zeros(npoints)
    H2O = zeros(npoints)
    SO3 = zeros(npoints)
    CaO_excess = zeros(npoints)
    volatiles = zeros(npoints)

    # Stoichiometric conversion factors 
    CaO_to_SO3 = (32.065+3*15.999)/(40.078+15.999)            # Gypsum / Anhydrite CaSO₄
    CaO_to_H2O = 2*(2*1.00784+15.999)/(40.078+15.999)         # Gypsum CaSO₄ ⋅ 2H₂O
    CaO_to_CO2 = (15.999*2+12.01)/(40.078+15.999)             # Assumes CaCO₃ decomposition

    K2O_to_Al2O3 = (26.98*2+15.999*3)/(39.098*2+15.999)       # K-spar
    Na2O_to_Al2O3 = (26.98*2+15.999*3)/(22.990*2+15.999)      # Na-plag
    Al2O3_to_CaO = (40.078+15.999)/(26.98*2+15.999*3)         # Ca-plag

    # [CARBONATES] Compute CO2 from CaO, assuming CaO is from CaCO3 
    # Replace values if there is no reported CO2 value. If the values are different, pick 
    # the larger of the two
    @turbo for i in eachindex(CO2)
        CO2[i] = ifelse(combo_cats.carb[i], out.CaO[i]*CaO_to_CO2, combined.CO2[i])
    end
    @turbo for i in eachindex(CO2)
        CO2[i] = nanmax(CO2[i], combined.CO2[i])
    end

    # [CLASTICS] Compute CO2 from non-feldspar CaO
    # Calculate the amount of CaO expected from Ca-plagioclase (assume all K, Na, and Al 
    # from alkali feldspar). Use the "excess" CaO to calculate CO2. Replace values as above.
    target = @. (combo_cats.siliciclast .| combo_cats.shale)
    @turbo for i in eachindex(CaO_excess)
        CaO_excess[i] = out.CaO[i] - ca_plagioclase(out.K2O[i], out.Al2O3[i], out.Na2O[i], 
            K2O_to_Al2O3, Na2O_to_Al2O3, Al2O3_to_CaO)
    end
    @turbo for i in eachindex(CO2)
        CO2[i] = ifelse(target[i] & (CaO_excess[i] > 0),        # Is clastic with excess CaO
            nanmax(CaO_excess[i]*CaO_to_CO2, combined.CO2[i]),  # Larger of calculated / reported 
            CO2[i]                                              # Else keep existing values
        )
    end

    # [EVAPORITES] Compute H2O [GYPSUM] and SO3 [GYPSUM, ANHYDRITE] from CaO
    # Take SO3 and H2O as the  larger of the calculated and reported values
    for i in eachindex(SO3)
        if combo_cats.evap[i] && combined.Rock_Name[i] == "gypsum"
            SO3[i] = nanmax(out.CaO[i]*CaO_to_SO3, combined.SO3[i])
            H2O[i] = nanmax(out.CaO[i]*CaO_to_H2O, combined.H2O[i])
        elseif combo_cats.evap[i] && combined.Rock_Name[i] == "anhydrite"
            SO3[i] = nanmax(out.CaO[i]*CaO_to_SO3, combined.SO3[i])
        end
    end

    # Calculate reported volatiles as the larger of (CO₂ + H₂O + SO₃) or LOI
    @turbo volatiles .= vec(nansum([CO2 H2O SO3], dims=2))
    @turbo for i in eachindex(volatiles)
        volatiles[i] = ifelse(combined.LOI[i] > volatiles[i], combined.LOI[i], volatiles[i])
    end


## --- Compute total wt.% analyzed for all samples
    # Exclude OIBs and samples on the oceanic crust 
    isOIB = findOIBs(combined.Latitude, combined.Longitude);

    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, combined.Latitude, combined.Longitude)
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


## --- Stoichiometric correction for unreported volatiles 
    # Preallocate 
    volatiles_assumed = Array{Float64}(undef, npoints)

    # Consider all sedimentary samples, and only those samples
    minorsed = get_rock_class()[2];
    for type in minorsed
        combo_cats.sed .|= combo_cats[type]
    end

    # If total analyzed wt.% below 100%, assume the unreported data are volatiles
    @turbo for i in eachindex(bulkweight)
        volatiles_assumed[i] = ifelse((combo_cats.sed[i] & (bulkweight[i] < 100)), (100 - bulkweight[i]), 0.0)
    end

    # Restrict data to between 84 and 104 wt.%
    t = @. 84 <= bulkweight <= 104; tᵢ = count(t)           # Without assumed volatiles
    t = @. 84 <= bulkweight .+ volatiles_assumed <= 104;    # With assumed volatiles
    
    @turbo volatiles .+= volatiles_assumed
    isgypsum = vec(combined.Rock_Name .== "gypsum");

    screenvolatiles!(vec(t), volatiles, isevap=combo_cats.evap, isgypsum=isgypsum,
        general=dol, evaporite=bas, gypsum=gyp    
    );

    # Print to terminal
    up = count(t) - tᵢ
    @info """Saving $(count(t)) samples ($(round(count(t)/length(t)*100, digits=2))%)
    Assuming volatiles increased count from $tᵢ to $(count(t))
    Total increase = $up ($(round(up/count(t)*100, sigdigits=2)) %)
    """


## --- Filter data and normalize to 100% 
    # Combine data into full output array (inefficient but insensitive to element order)
    out = merge(out, (
            Volatiles=volatiles,
            Latitude = combined.Latitude,
            Longitude = combined.Longitude,
            Loc_Prec = combined.Loc_Prec,
            Age = combined.Age,
            Age_Max = combined.Age_Max,
            Age_Min = combined.Age_Min,
            Sample_ID = collect(1:npoints)
        )
    )
    out = NamedTuple{keys(out)}([out[k][t] for k in keys(out)])

    # Normalize to 100%
    p = Progress(length(out.SiO2) ÷ 10, desc="Normalizing compositions ...", enabled=show_progress)
    contributing = [majors; minors; :Volatiles]          # Include volatiles, exclude isotopes!
    TOC = combined.TOC[t]                                # Include TOC in normalization...
    for i in eachindex(out.SiO2)
        sample = [out[j][i] for j in contributing]       # Get it
        normalize!([sample; TOC[i]])                     # Normalize it
        for j in eachindex(contributing)                 # Put it back
            out[contributing[j]][i] = sample[j]
        end
        i%10==0 && next!(p)
    end

    # Remove impossible locations 
    out.Latitude[(-90 .> out.Latitude) .| (out.Latitude .> 90)] .= NaN;
    out.Longitude[(-180 .> out.Longitude) .| (out.Longitude .> 180)] .= NaN;


## --- Get and restrict metadata 
    # Rock classes 
    combo_kittens = NamedTuple{keys(combo_cats)}(combo_cats[k][t] for k in keys(combo_cats))
    for k in keys(combo_cats)
        if show_progress && count(combo_cats[k]) < 10
            @warn "After filtering, type \"$k\" has $(count(combo_cats[k])) samples."
        end
    end
    meta_kittens = NamedTuple{keys(meta_cats)}(meta_cats[k][t] for k in keys(meta_cats));


## --- Save all data to new file 
    fid = h5open(fileout, "w")
    data = create_group(fid, "bulk")
    bulktext = create_group(fid, "bulktext")
    class = create_group(fid, "bulktypes")

    # Data 
    allkeys = collect(keys(out))
    write(data, "header", string.(allkeys))
    writebulk = create_dataset(data, "data", Float64, (count(t), length(allkeys)))
    for i in eachindex(allkeys)
        writebulk[:,i] = out[allkeys[i]]
    end

    # References and rock descriptions
    write(bulktext, "Reference", combined.Reference[t])
    write(bulktext, "Data_Source", combined.Data_Source[t])
    write(bulktext, "Methods", combined.Methods[t])
    write(bulktext, "Rock_Name", combined.Rock_Name[t])

    # Rock class
    a = Array{Int64}(undef, length(combo_kittens[1]), length(combo_kittens))
    for i in eachindex(keys(combo_kittens))
        for j in eachindex(combo_kittens[i])
            a[j,i] = ifelse(combo_kittens[i][j], 1, 0)
        end
    end

    b = similar(a)
    for i in eachindex(keys(meta_kittens))
        for j in eachindex(meta_kittens[i])
            b[j,i] = ifelse(meta_kittens[i][j], 1, 0)
        end
    end

    class["bulk_cats"] = a
    class["metamorphic_cats"] = b
    class["bulk_cats_head"] = string.(collect(keys(combo_kittens)))
    
    close(fid)


## --- Print relevant methods data to terminal 
    # Total number usable samples (descriptive lithology)
    exclude_minor!(combo_kittens);
    unmatched = RockWeatheringFlux.find_unmatched(combo_kittens)
    nondescriptive = combo_kittens.sed .| combo_kittens.ign .| combo_kittens.met .| combo_kittens.cover
    not_used = count(unmatched) + count(nondescriptive)

    include_minor!(combo_kittens)
    # invalidvolc = round(count(combo_kittens.volc .& invalidage[t]) / count(combo_kittens.volc)*100, sigdigits=3)
    # invalidplut = round(count(combo_kittens.plut .& invalidage[t]) / count(combo_kittens.plut)*100, sigdigits=3)
    invalidseds = round(count(combo_kittens.sed .& invalidage[t]) / count(combo_kittens.sed)*100, sigdigits=3)

    @info """Functional dataset size: $(length(combo_kittens.sed) - not_used)
    Nondescriptive: $not_used
    Unused lithology: $(count(unmatched))
    Total samples: $(length(combo_kittens.sed))
    Invalid age included: $(count(invalidage[t]))
        Sedimentary (set to NaN): $invalidseds%
    """


## --- End of File     