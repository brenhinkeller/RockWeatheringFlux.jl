## --- Setup
    # Calculate the flux of each element at each point, as well as the composition of 
    # total produced sediment (bulk and by lithologic class).
    
    # Also calculate the total surficial abundance of each lithologic class, for comparison
    # between relative contribution to mass flux of eroded material and surficial abundance.

    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, Dates
    using StatsBase


## --- Load data
    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    recorded_type  = vec(fid[:,2])
    t = @. gchem_ind != 0
    
    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)

    # Location
    fid = h5open(macrostrat_io, "r")
    rocklat = read(fid["vars"]["rocklat"])[t]
    rocklon = read(fid["vars"]["rocklon"])[t]
    close(fid)


## --- Calculate erosion rate at each coordinate point of interest	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point, exclude values above 1000 m/km
    rockslope = movingwindow(srtm15_slope, rocklat, rocklon, srtm15_sf, n=5)
    rockslope = Measurements.value.(rockslope)
    rockslope[rockslope .>= 1000] .= NaN;

    # Calculate all erosion rates [mm/kyr] (exclude erosion > 10_000 mm/kyr)
    rock_ersn = emmkyr.(rockslope);
    rock_ersn[rock_ersn .> 10_000] .= NaN


## --- Calculate bulk erosion and erosion by element at each point 
    # Definitions 
    crustal_density = 2750                                  # kg/m³
    npoints = length(mbulk.SiO2)                            # Number of *matched* points
    unit_sample_area = (133_367_840 * 1000000) / npoints    # Area of continents / npoints (m²)
        # Area considered = Total - (Antarctica + Greenland)
        # 149,733,926 - (14,200,000 + 2,166,086) = 133_367_840

    # Get element names
    majors, minors = get_elements()
    allelements = [majors; minors]
    nelements = length(allelements)
    
    # Bulk (total / undifferentiated) erosion at each point
    # Multiply by 1e-6 to convert from kg/Myr to kg/yr
    erosion_bulk = [rock_ersn[i] * unit_sample_area * crustal_density * 1e-6 for i = 1:npoints];

    # Contribution of each element to bulk erosion at each point
    # Multiply by 1e-2 to convert from wt.% to fraction
    erosion_element = (;
        vals = Array{Float64}(undef, npoints, nelements),
        errs = Array{Float64}(undef, npoints, nelements),
    )
    for i in eachindex(allelements)
        elementflux = [erosion_bulk[j] * mbulk[allelements[i]][j] * 1e-2 for j = 1:npoints]
        erosion_element.vals[:,i] .= Measurements.value.(elementflux)
        erosion_element.errs[:,i] .= Measurements.uncertainty.(elementflux)
    end

    # Save to file
    fid = h5open(eroded_out, "w")

    g = create_group(fid, "erosion_bulk")
    write(g, "values", Measurements.value.(erosion_bulk))
    write(g, "errors", Measurements.uncertainty.(erosion_bulk))

    g = create_group(fid, "erosion_element")
    write(g, "header", string.(allelements))
    write(g, "values", erosion_element.vals)
    write(g, "errors", erosion_element.errs)

    close(fid)


## --- Load data from file 
    fid = h5open("$eroded_out", "r")

    # Bulk erosion at each point [kg/yr]
    erosion_bulk = (;
        vals = read(fid["erosion_bulk"]["values"]),
        errs = read(fid["erosion_bulk"]["errors"]),
    )

    # Erosion by element at each point [kg/yr]
    header = Tuple(Symbol.(read(fid["erosion_element"]["header"])))
    vals = read(fid["erosion_element"]["values"])
    errs = read(fid["erosion_element"]["errors"])
    erosion_element = (;
        vals = NamedTuple{header}([vals[:,i] for i in eachindex(header)]),
        errs = NamedTuple{header}([errs[:,i] for i in eachindex(header)]),
    )

    close(fid)


## --- Calculate the total (global) mass of eroded sediment and of each element
    # Convert point data from kg/yr to Gt/yr
    kg_gt = 1e12

    # Total mass of eroded sediment (one number for whole earth)
    global_erosion_bulk = (;
        vals = nansum(erosion_bulk.vals ./ kg_gt),
        errs = sqrt(nansum((erosion_bulk.errs ./ kg_gt).^2)),
    )

    # Total mass of each eroded element (one number for each element)
    elem = keys(erosion_element.vals)
    global_erosion_element = (;
        vals = NamedTuple{elem}([nansum(erosion_element.vals[k]) / kg_gt for k in elem]),
        errs = NamedTuple{elem}([sqrt(nansum((erosion_element.errs[k] ./ kg_gt).^2)) for k in elem]),
    );

    # Print to terminal 
    global_sum = (;
        val = round(global_erosion_bulk.vals, sigdigits=3),
        err = round(2*global_erosion_bulk.errs, sigdigits=3),
    )
    global_element_sum = (;
        val = round(nansum(values(global_erosion_element.vals)), sigdigits=3),
        err = round(sqrt(2*nansum(values(global_erosion_element.errs).^2)), sigdigits=3),
    )

    # Print to terminal, warn if sum of elements does not equal bulk erosion rate
    if isapprox(global_sum.val, global_element_sum.val, atol = max(global_sum.err, global_element_sum.err))
        @info "Mass of global sediment production ± 2σ s.d.: $(global_sum.val) ± $(global_sum.err) Gt/yr"
    else
        diff = abs((global_sum.val ± global_sum.err) - (global_element_sum.val ± global_element_sum.err))
        @warn """
        Mass of global sediment production ± 2σ s.d.: $(global_sum.val) ± $(global_sum.err) Gt/yr

        Mass of eroded sediment and sum of eroded elements are not within 2σ s.d. error.
        \t Eroded sediment:        $(global_sum.val) ± $(global_sum.err) Gt/yr
        \t Sum of eroded elements: $(global_element_sum.val) ± $(global_element_sum.err) Gt/yr
        Difference: $diff Gt/yr
        """
    end


## --- Preallocate and set switches for export
    # [SWITCH] lithologic class filter 
    classfilter = megaclass
    classes = keys(classfilter)

    # Elements of interest 
    majors, minors = get_elements()
    nmajors = length(majors)
    allelements = [majors; minors]

    # N samples account for 90% of matches
    npoints = unique_sample(mbulk.Sample_ID, 90)

    # Preallocate results array (element row, rock class column)
    rows = vcat(string.(collect(keys(erosion_element.vals))), "Total")
    cols = hcat("elem", reshape(string.(collect(keys(classfilter))), 1, :));

    result = (;
        vals = Array{Float64}(undef, length(rows), length(cols) - 1),
        errs = Array{Float64}(undef, length(rows), length(cols) - 1),       # 1σ s.d.
    )


## --- Save to file: Absolute mass of eroded elements by rock class  
    # The absolute mass of eroded material which erodes from a given rock class is the sum
    # of the mass of each point of that class 
    for i in eachindex(classes)
        for j in eachindex(allelements)
            rockclass = classfilter[classes[i]]     # For each rock class...
            element = allelements[j]                # For each element...
            
            # Sum all values which are matched to the rock class of interest
            result.vals[j,i] = nansum(erosion_element.vals[element][rockclass])
            result.errs[j,i] = sqrt(nansum(erosion_element.vals[element][rockclass]).^2)

            # Convert sum from kg/yr to Gt/yr
            result.vals[j,i] /= kg_gt
            result.errs[j,i] /= kg_gt
        end
    end

    # The bulk mass eroded from each rock class is the sum of the mass of all elements 
    result.vals[end,:] .= vec(nansum(result.vals[1:end-1,:], dims=1))
    result.errs[end,:] .= vec(sqrt.(nansum((result.errs[1:end-1,:]).^2, dims=1)))

    # Replace all 0's with NaNs!!
    nanzero!(result.vals)
    nanzero!(result.errs)

    # Save to file
    writedlm(erodedabs_out, vcat(cols, hcat(rows, result.vals)))
    writedlm(erodedabs_out_err, vcat(cols, hcat(rows, result.errs)))

    # Save major element values to print to terminal 
    erodedmass_to_terminal = (;
        vals = round.(result.vals[1:nmajors, end], digits=1),
        errs = round.(result.errs[1:nmajors, end]*2, digits=1),
    )


## --- Save to file: Fractional contribution of each rock class and element to total erosion
    # That is, the fraction of total erosion represented by each element and each class.
    # Rows are summative rather than columns; e.g. a sedimentary SiO2 value of 0.61 means 
    # 61% of SiO2 erodes from sedimentary rocks---NOT that 61% of sedimentary eroded 
    # eroded material is SiO2 (the latter is composition, computed later). 

    # The contribution of major classes will be the sum of the contribution of their 
    # subclasses. Keep in mind that Because we don't consider metamorphic rocks as a 
    # descriptive class (e.g., during matching, all rocks are assigned to be a subclass 
    # of sed / volc / plut rocks), the sum of the fractional contributions of 
    # sed + volc + plut = 1
    
    # Calculate contribution 
    result_contribution = (;
        vals = result.vals ./ result.vals[:,end],
        errs = result.errs ./ result.vals[:,end],
    )

    # Save to file 
    writedlm(erodedrel_out, vcat(cols, hcat(rows, result_contribution.vals)))
    writedlm(erodedrel_out_err, vcat(cols, hcat(rows, result_contribution.errs)))

    # Format total contribution for terminal printout 
    contribution = NamedTuple{classes}(result_contribution.vals[end,:].*100)


## --- Save to file: Composition of eroded material 
    # Calcuating the composition of eroded material from the total eroded mass means 
    # infrequently-measured elements (e.g. REEs) will be skewed downward due to missing 
    # data. 

    # TO DO: propagate uncertainty correctly :(
    @warn "Composition of eroded material does not correctly propagate uncertainties"

    # Compute the average total erosion at each point [kg/yr] for each rock class
    erosion_bulk_average = (;
        vals = NamedTuple{classes}(nanmean(erosion_bulk.vals[classfilter[k]]) for k in classes),
        errs = NamedTuple{classes}(nanmean(erosion_bulk.errs[classfilter[k]]) for k in classes),
    );

    # Compute the average erosion at each point [kg/yr] for each element for each rock class 
    erosion_element_average = (;
        vals = NamedTuple{classes}([NamedTuple{Tuple(allelements)}(
            nanmean(erosion_element.vals[e][classfilter[k]]) for e in allelements) for k in classes]
        ),
        errs = NamedTuple{classes}([NamedTuple{Tuple(allelements)}(
            nanmean(erosion_element.errs[e][classfilter[k]]) for e in allelements) for k in classes]
        ),
    );

    # The average composition of each class is average element mass / average total mass
    result_composition = deepcopy(result)
    for i in eachindex(classes)
        for j in eachindex(allelements)
            element = allelements[j]                # For each element...
            
            # This uncertainty propogates correctly!
            val = erosion_element_average.vals[classes[i]][element] / erosion_bulk_average.vals[classes[i]]
            err = erosion_element_average.errs[classes[i]][element] / erosion_bulk_average.vals[classes[i]]
            
            # Convert to wt.%
            result_composition.vals[j,i] = val * 100
            result_composition.errs[j,i] = err * 100
        end
    end

    # Save to file
    writedlm(erodedcomp_out, vcat(cols, hcat(rows, result_composition.vals)))
    writedlm(erodedcomp_out_err, vcat(cols, hcat(rows, result_composition.errs .* 100)))

    # Format for terminal printout 
    composition = NamedTuple{classes}([(
        vals = round.([result_composition.vals[:,k][i] for i in eachindex(majors)], sigdigits=3),
        errs = round.([result_composition.errs[:,k][i] ./ sqrt(npoints) .*2 for i in eachindex(majors)], sigdigits=1)
    ) for k in eachindex(classes)])


## --- Print to terminal: major element mass and composition vibe check
    @info """ Annual mass flux of denuded material:
    Total mass flux [Gt/yr] ± 2σ s.d.
      $(join(rpad.(majors, 8), " "))
      $(join(rpad.(erodedmass_to_terminal.vals, 8), " "))
    ± $(join(rpad.(erodedmass_to_terminal.errs, 8), " "))
    
    Composition [wt.%] of bulk eroded material ± 2 s.e.:
      $(join(rpad.(majors, 8), " "))
      $(join(rpad.(composition.bulk.vals, 8), " "))
    ± $(join(rpad.(composition.bulk.errs, 8), " "))
    """

## --- Print to terminal: composition, formatted for LaTeX / Excel workflow
    # Preallocate
    target = (:sed, :volc, :plut, :bulk)
    out = fill("", length(majors) + 1)

    # Compute anhydrous composition of eroded material 
    @assert majors[end] == :Volatiles   # Assumes volatiles listed at the end of majors
    anhydrous_sum = 100 - composition.bulk.vals[end]
    composition_anhydrous = (;
        vals = round.(composition.bulk.vals[1:end-1] ./ anhydrous_sum .* 100, sigdigits=3),
        errs = round.(composition.bulk.vals[1:end-1] ./ anhydrous_sum .* 100, sigdigits=1),
    )

    # Format major element composition and sum of major elements 
    for key in target
        for i in eachindex(majors) 
            out[i] *= "\$ $(composition[key].vals[i]) \\pm $(composition[key].errs[i]) \$; "
        end 
        out[end] *= "$(round(sum(composition[key].vals), sigdigits=4)); "
    end
    
    # Format anhydrous composition and sum of major elements 
    out_anhydrous = fill("", length(majors)+1)
    for i in eachindex(majors[1:end-1])
        out_anhydrous[i] *= "\$ $(composition_anhydrous.vals[i]) \\pm $(composition_anhydrous.errs[i]) \$"
    end
    out_anhydrous[end] = string(round(sum(composition_anhydrous.vals), sigdigits=4))

    # Print to terminal 
    @info "Composition of eroded material: $target + anhydrous"
    for i in eachindex(out)
        println("$(out[i] * out_anhydrous[i])")
    end


## --- Save to file / Print to terminal: Surficial abundance and fractional contribution
    # This repeats data found in the fractional contribution to eroded material file, but
    # simplified to only bulk data (e.g. undifferentiated by element)

    # To compare surficial abundance and contribution to total erosion, we must use the 
    # lithologies assigned to each sample during matching. This means the total surficial 
    # abundance of sed + ign = 100%. Metamorphic abundances should be considered as subsets 
    # of "descriptive" lithologies: e.g. 7% of sedimentary rocks are metasedimentary. The 
    # "met" class is all metamorphic rocks (metased + metaign + undifferentiated)

    # Calculate surficial abundance 
    matched_surficial = NamedTuple{classes}(
        count(classfilter[k])/length(classfilter[k])*100 for k in keys(classfilter)
    )

    # Save to file 
    cols = ["lithology" "surficial abundance" "fractional contribution"]
    writedlm(surfacelith_calc_out, vcat(cols, hcat(
        collect(string.(classes)), 
        collect(values(matched_surficial)), 
        collect(values(contribution))
    )))

    # Print data for target classes to terminal 
    @info "Surficial abundance / fractional contribution to erosion: $target"
    println(
        """$(join([round(matched_surficial[k], sigdigits=3) for k in target], "; "))
        $(join([round(contribution[k], sigdigits=3) for k in target], "; "))
        """
    )


## --- Save to file: Surficial abundance as mapped
    # Load Macrostrat lithologic classes (matched samples only)
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    fid = h5open(macrostrat_io, "r")
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    macro_cats = delete_cover(macro_cats)
    close(fid)

    # Definitions 
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:end];
    npoints = count(t)

    # Find the location (continent) of each sample and the relative contribution of each 
    # continent to total samples. Exclude NA
    continent = find_geolcont(rocklat, rocklon)
    continent = NamedTuple{Tuple([Symbol.(continents[1:end-1]); :Global])}([
        [continent .== i for i in eachindex(continents[1:end-1])]; [trues(npoints)]
    ])
    dist_cont = merge(NamedTuple{keys(continent)[1:end-1]}(normalize!(
            [count(continent[k])/npoints for k in keys(continent)[1:end-1]])
        ), (Global=100.0,)
    )

    # Preallocate 
    region = collect(keys(continent))
    ncols = (length(minorsed)+2) + (length(minorplut)+2) + (length(minorvolc)+2) + 3 + 4;
    result = Array{Float64}(undef, ncols, length(region));

    # Calculate surficial abundance by region
    for i in eachindex(region)
        # Filter for the region
        s = continent[region[i]]
        nregion = count(s)

        # Calculate abundances of constituent classes in this region, normalized to 100%
        include_minor!(macro_cats);
        ign_undiff = .!(macro_cats.volc .| macro_cats.plut .| macro_cats.carbonatite);
        dist2 = NamedTuple{(:sed, :volc, :plut, :carbonatite, :ign_undiff, :met_undiff)}(
            normalize!([
                count(macro_cats.sed .& s),                 # Sedimentary
                count(macro_cats.volc .& s),                # Volcanic
                count(macro_cats.plut .& s),                # Plutonic
                count(macro_cats.carbonatite .& s),         # Carbonatite
                count(macro_cats.ign .& ign_undiff .& s),   # Undifferentiated ign
                count(megaclass.met_undiff .& s)            # Undifferentiated met
            ]./nregion.*100
        ))
        dist2_ign = dist2.volc + dist2.plut + dist2.carbonatite + dist2.ign_undiff

        # Count frequency of subclasses within constituent classes. 
        exclude_minor!(macro_cats)
        sed = float.([[count(macro_cats[k] .& s) for k in minorsed]; count(macro_cats.sed .& s)])
        volc = float.([[count(macro_cats[k] .& s) for k in minorvolc]; count(macro_cats.volc .& s)])
        plut = float.([[count(macro_cats[k] .& s) for k in minorplut]; count(macro_cats.plut .& s)])

        # Calculate percentage abundance of each class from counts 
        sed .= sed ./ nregion * 100
        volc .= volc ./ nregion * 100
        plut .= plut ./ nregion * 100
        
        # Normalize minor class percentages to the abundance of the constituent class
        sed .= sed ./ sum(sed) .* dist2.sed
        volc .= volc ./ sum(volc) .* dist2.volc
        plut .= plut ./ sum(plut) .* dist2.plut
        
        # Calculate the percentage of metased and metaign rocks as the fraction of total 
        # sed / ign rocks tagged as metamorphic 
        include_minor!(macro_cats);
        metased = count(megaclass.metased .& s) / count(macro_cats.sed .& s)
        metaign = count(megaclass.metaign .& s) / count(macro_cats.ign .& s)
        metased *= dist2.sed 
        metaign *= dist2_ign

        met_total = metased + metaign + dist2.met_undiff

        # Convert data from 100% of region to % of total surface area 
        area_frac = dist_cont[region[i]] / 100

        ign_out = [dist2_ign, dist2.carbonatite, dist2.ign_undiff] .* area_frac
        sed_out = [sum(sed); sed] .* area_frac
        volc_out = [sum(volc); volc] .* area_frac
        plut_out = [sum(plut); plut] .* area_frac
        met_out = [met_total, metased, metaign, dist2.met_undiff] .* area_frac

        # Save all data to the results array in the order [total; subtypes] 
        # Be SURE to check that this is in the same order as labeled in the output 🤦
        result[:,i] .= [sed_out; volc_out; plut_out; ign_out; met_out]
    end
    
    # Normalize each set of results to the global total 
    for i = 1:size(result)[1]
        normresult = result[i, 1:end-1] ./ sum(result[i, 1:end-1]) * result[i,end]
        result[i, 1:end-1] .= normresult
    end

    # Define labels (down here so it's easier to compare to the arrays above)
    sed_label = ["Total Sedimentary"; string.(collect(minorsed)); "Undifferentiated Sedimentary"]
    volc_label = ["Total Volcanic"; string.(collect(minorvolc)); "Undifferentiated Volcanic"]
    plut_label = ["Total Plutonic"; string.(collect(minorplut)); "Undifferentiated Plutonic"]
    ign_label = ["Total Igneous", "Carbonatite", "Undifferentiated Igneous"]
    met_label = ["Total Metamorphic", "Metasedimentary", "Metaigneous", "Undifferentiated Metamorphic"]
    
    cols = ["Lithology" reshape(collect(string.(keys(continent))), 1, :)]

    # Export file, values in percentages
    labels = [sed_label; volc_label; plut_label; ign_label; met_label];
    writedlm(surfacelith_mapped_out, vcat(cols, hcat(labels, result)))

    # Check to make sure the sums work out 
    sums_match = isapprox.(nansum(result[:, 1:end-1], dims=2), result[:, end])
    if count(sums_match) != size(result)[1]
        @warn "Mapped lithology: discrepancy in values for: $(labels[.!vec(sums_match)])" 
    end

    # Print % undifferentiated lithologies to terminal 
    target = containsi.(labels, "Undifferentiated")
    undiff_label = labels[target]
    undiff_value = round.(result[:,end][target], sigdigits=3)
    @info """ Abundance of:
    $(undiff_label[1]): $(undiff_value[1])%
    $(undiff_label[2]): \t$(undiff_value[2])%
    $(undiff_label[3]): \t$(undiff_value[3])%
    $(undiff_label[4]): \t$(undiff_value[4])%
    $(undiff_label[5]): $(undiff_value[5])%
    """


## --- End of File