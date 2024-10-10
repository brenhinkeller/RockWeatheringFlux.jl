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
    srtm15_slope = h5read(srtm_maxslope, "vars/slope")
    srtm15_sf = h5read(srtm_maxslope, "vars/scalefactor")

    # Get slope at each coordinate point, exclude values above 1000 m/km
    rockslope = movingwindow(srtm15_slope, rocklat, rocklon, srtm15_sf, n=5)
    rockslope = Measurements.value.(rockslope)
    rockslope[rockslope .>= 1000] .= NaN;
    
    # # Optionally save and load data from a tsv if your computer is small and weak 
    # writedlm("output/basins/rockslope_temp.tsv", rockslope)
    # rockslope = readdlm("output/basins/rockslope_temp.tsv")

    # Calculate all erosion rates [mm/kyr] (exclude erosion > 10_000 mm/kyr)
    # Gives us 1-σ standard deviations
    rock_ersn = emmkyr.(rockslope);
    rock_ersn[rock_ersn .> 10_000] .= NaN;


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
    fid = h5open(eroded_out, "r")

    # Bulk erosion at each point [kg/yr]
    erosion_bulk = (;
        vals = read(fid["erosion_bulk"]["values"]),
        errs = read(fid["erosion_bulk"]["errors"]),     # 1σ s.d.
    )

    # Erosion by element at each point [kg/yr]
    header = Tuple(Symbol.(read(fid["erosion_element"]["header"])))
    vals = read(fid["erosion_element"]["values"])
    errs = read(fid["erosion_element"]["errors"])
    erosion_element = (;
        vals = NamedTuple{header}([vals[:,i] for i in eachindex(header)]),
        errs = NamedTuple{header}([errs[:,i] for i in eachindex(header)]),   # 1σ s.d.
    )

    close(fid)


## --- Calculate the total (global) mass of eroded sediment and of each element
    # Propagation/ manipulation of error values is marked with `Propagate error`

    # Convert point data from kg/yr to Gt/yr
    kg_gt = 1e12

    # Total mass of eroded sediment (one number for whole earth)
    # Propagate error: sum, division by exact number
    global_erosion_bulk = (;
        vals = nansum(erosion_bulk.vals ./ kg_gt),
        errs = sqrt(nansum((erosion_bulk.errs ./ kg_gt).^2)),
    )

    # Total mass of each eroded element (one number for each element)
    # Propagate error: sum, division by exact number
    elem = keys(erosion_element.vals)
    global_erosion_element = (;
        vals = NamedTuple{elem}([nansum(erosion_element.vals[k]) / kg_gt for k in elem]),
        errs = NamedTuple{elem}([sqrt(nansum((erosion_element.errs[k] ./ kg_gt).^2)) for k in elem]),
    );

    # Print to terminal (these don't get used anywhere else)
    # Propagate error for global_element_sum: sum
    global_sum = (;
        val = round(global_erosion_bulk.vals, sigdigits=3),
        err = round(2*global_erosion_bulk.errs, sigdigits=3),                               # 2σ s.d.
    )
    global_element_sum = (;
        val = round(nansum(values(global_erosion_element.vals)), sigdigits=3),
        err = round(2*sqrt(nansum(values(global_erosion_element.errs).^2)), sigdigits=3),   # 2σ s.d.
    )

    # Print to terminal, warn if sum of elements does not equal bulk erosion rate
    # TODO: s.d or s.e.m.?
    if isapprox(global_sum.val, global_element_sum.val, atol = max(global_sum.err, global_element_sum.err))
        @info "Mass of global sediment production ± 2σ s.d.: $(global_sum.val) ± $(global_sum.err) Gt/yr"
    else
        diff = abs((global_sum.val ± global_sum.err) - (global_element_sum.val ± global_element_sum.err))
        @warn """
        Mass of global sediment production ± 2σ s.d.: $(global_sum.val) ± $(global_sum.err) Gt/yr

        Mass of eroded sediment and sum of eroded elements are not within 2σ s.d.
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
            # Propagate error: sum
            # Note: result.errs was previously using erosion_element.vals. I assume this was a typo?
            result.vals[j,i] = nansum(erosion_element.vals[element][rockclass])
            result.errs[j,i] = sqrt(nansum(erosion_element.errs[element][rockclass]).^2)

            # Convert sum from kg/yr to Gt/yr
            # Propagate error: division by exact number
            result.vals[j,i] /= kg_gt
            result.errs[j,i] /= kg_gt
        end
    end

    # The bulk mass eroded from each rock class is the sum of the mass of all elements 
    # Propagate error: sum
    result.vals[end,:] .= vec(nansum(result.vals[1:end-1,:], dims=1))
    result.errs[end,:] .= vec(sqrt.(nansum((result.errs[1:end-1,:]).^2, dims=1)))

    # Replace all 0's with NaNs!!
    nanzero!(result.vals)
    nanzero!(result.errs)

    # Convert 1σ s.d. to 2σ s.e.m. 
    # Propagate error: division and multiplication by exact number 
    result.errs ./= sqrt(npoints)   # s.d -> s.e.m
    result.errs .*= 2               # 1σ  -> 2σ

    # Save to file
    writedlm(erodedmass_out, vcat(cols, hcat(rows, result.vals)))
    writedlm(erodedmass_out_err, vcat(cols, hcat(rows, result.errs)))

    # Save to publication-formatted file 
    pubcols = hcat(cols[1], mesh(cols[:,2:end], fill("+/- 2 SEM", size(cols))[:,2:end]))
    pubresult = mesh(result.vals, result.errs)
    writedlm(erodedmass_out_csv, vcat(pubcols, hcat(rows, pubresult)), ',')

    # Save major element values to print to terminal 
    erodedmass_to_terminal = (;
        vals = round.(result.vals[1:nmajors, end], digits=1),
        errs = round.(result.errs[1:nmajors, end], digits=1),
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

    # For context, we also want to print the total abundance of that lithology in our 
    # dataset. To compare surficial abundance and contribution to total erosion, we use
    # the lithology assigned to each sample during matching. This means the total surficial 
    # abundance of sed + ign = 100%. Metamorphic abundances are considered subsets 
    # of "descriptive" lithologies: e.g. 7% of sedimentary rocks are metasedimentary. The 
    # "met" class is all metamorphic rocks (metased + metaign + undifferentiated)

    
    # Calculate contribution 
    # Propagate error: division by exact number
    result_contribution = (;
        vals = result.vals ./ result.vals[:,end],
        errs = result.errs ./ result.vals[:,end],    # Error is already 2 s.e.m.
    )

    # Calculate surficial abundance 
    matched_surficial = NamedTuple{classes}(
        count(classfilter[k])/length(classfilter[k]) for k in keys(classfilter)
    )
    surf = reshape(collect(matched_surficial), 1, :)
    surf_err = fill("", size(surf))                     # Empty errors

    # Check column order will print correctly
    @assert collect(string.(keys(matched_surficial))) == reshape(cols, :, 1)[2:end] ":("

    # Save to file 
    rows2 = [rows[1:end-1]; "Undifferentiated"; "Lithology Abundance"]
    writedlm(frac_contributed, vcat(cols, hcat(rows2, vcat(result_contribution.vals, surf))))
    writedlm(frac_contributed_out_err, vcat(cols, hcat(rows2, vcat(result_contribution.errs, surf))))

    # Save to publication-formatted file 
    pubcols = hcat(cols[1], mesh(cols[:,2:end], fill("+/- 2 SEM", size(cols))[:,2:end]))
    pubresult = mesh(result_contribution.vals, result_contribution.errs)
    pubsurf = mesh(surf, surf_err)
    writedlm(frac_contributed_out_csv, vcat(pubcols, hcat(rows2, vcat(pubresult, pubsurf))), ',')

    # Print to terminal
    target = (:sed, :volc, :plut, :bulk)
    contribution = NamedTuple{classes}(result_contribution.vals[end,:].*100)
    @info "Surficial abundance / fractional contribution to erosion: $target"
    println(
        """$(join([round(matched_surficial[k], sigdigits=3) for k in target], "; "))
        $(join([round(contribution[k], sigdigits=3) for k in target], "; "))
        """
    )


## --- Save to file: Composition of eroded material 
    # Calcuating the composition of eroded material from the total eroded mass means 
    # infrequently-measured elements (e.g. REEs) will be skewed downward due to missing 
    # data. 

    # Compute average total erosion at each point [kg/yr] for each rock class (2 s.e.m errors) 
    # Propagate error: mean (sum, division by an exact number)
    # Propagate error: s.d. -> s.e.m. division and multiplication by exact number 
    erosion_bulk_average = (;
        vals = NamedTuple{classes}(nanmean(erosion_bulk.vals[classfilter[k]]) for k in classes),
        errs = NamedTuple{classes}(
            ((sqrt(nansum(erosion_bulk.errs[classfilter[k]]).^2) / count(classfilter[k]))
                / sqrt(npoints)    # s.d -> s.e.m
                * 2                # 1σ  -> 2σ
            )
        for k in classes),
    );

    # Compute the average erosion at each point [kg/yr] for each element for each rock class 
    # Propagate error: mean (sum, division by an exact number)
    # Propagate error: s.d. -> s.e.m. division and multiplication by exact number 
    erosion_element_average = (;
        vals = NamedTuple{classes}([NamedTuple{Tuple(allelements)}(
            nanmean(erosion_element.vals[e][classfilter[k]]) for e in allelements) for k in classes]
        ),
        errs = NamedTuple{classes}([NamedTuple{Tuple(allelements)}(
            ((sqrt(nansum(erosion_element.errs[e][classfilter[k]]).^2) / count(classfilter[k])) 
                / sqrt(npoints)    # s.d -> s.e.m
                * 2                # 1σ  -> 2σ
            ) for e in allelements) for k in classes]
        ),
    );

    # Average composition of each class (element mass / total mass)
    result_composition = deepcopy(result)           # We won't use values, just the structure
    for i in eachindex(classes)
        for j in eachindex(allelements)
            element = allelements[j]                # For each element...
            
            # Average element mass / average total mass
            # Propagate error: divison by exact number 
            val = erosion_element_average.vals[classes[i]][element] / erosion_bulk_average.vals[classes[i]]
            err = erosion_element_average.errs[classes[i]][element] / erosion_bulk_average.vals[classes[i]]
            
            # Convert to wt.%
            # Propagate error: multiplication by exact number
            result_composition.vals[j,i] = val * 100
            result_composition.errs[j,i] = err * 100
        end
    end

    # Save to file
    writedlm(comp_eroded, vcat(cols, hcat(rows, result_composition.vals)))
    writedlm(comp_eroded_err, vcat(cols, hcat(rows, result_composition.errs)))

    # Save to publication-formatted file 
    pubcols = hcat(cols[1], mesh(cols[:,2:end], fill("+/- 2 SEM", size(cols))[:,2:end]))
    pubresult = mesh(result_contribution.vals, result_contribution.errs)
    writedlm(comp_eroded_csv, vcat(pubcols, hcat(rows, pubresult)), ',')

    # Format for terminal printout 
    composition = NamedTuple{classes}([(
        vals = round.([result_composition.vals[:,k][i] for i in eachindex(majors)], sigdigits=3),
        errs = round.([result_composition.errs[:,k][i] for i in eachindex(majors)], sigdigits=1)
    ) for k in eachindex(classes)]);


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


## --- End of File