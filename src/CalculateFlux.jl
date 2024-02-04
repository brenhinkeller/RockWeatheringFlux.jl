## --- Setup
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, Dates


## --- Load data
    # Indices of matched samples from SampleMatch.jl
    fid = readdlm(matchedbulk_io)
    matches = Int.(vec(fid[:,1]))
    t = @. matches != 0

    # Rock class of each sample
    match_cats = match_rocktype(string.(vec(fid[:,2]))[t]);
    include_minor!(match_cats)

    # Geochemical data for each sample
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(fid)

    # Lithology of each sample
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
    )
    close(fid)


## --- Calculate erosion rate at each coordinate point of interest	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)


## --- Calculate denundation at each point
    # Area of land:
        # Total: 149,733,926
        # Antarctica: 14,200,000 
        # Greenland: 2,166,086
        # Area in this study = Total - (Antarctica + Greenland) = 133_367_840
    # Declare constants
    const npoints = length(macrostrat.rocklat)                  # Number of points
    const crustal_density = 2750                                # kg/m³
    const unit_sample_area = (133367840 * 1000000) / npoints  # Area of continents / npoints (m²)

    # Create file to save data
    fid = h5open("$eroded_out", "w")

    # Denundation at each point (kg/yr)
    sampleflux = [rock_ersn[i] * unit_sample_area * crustal_density * 1e-6 for i = 1:npoints]

    # Save to file
    sampleflux_val, sampleflux_err = unmeasurementify(sampleflux)
    bulk_denundation = create_group(fid, "bulk_denundation")
    write(bulk_denundation, "values", sampleflux_val)
    write(bulk_denundation, "errors", sampleflux_err)


## --- Calculate flux of each element at each point
    # Define elements
    majors, minors = get_elements()
    allelements = [majors; minors]
    nelements = length(allelements)

    # Preallocate file space
    element_flux = create_group(fid, "element_flux")
    elem_vals = create_dataset(element_flux, "values", Float64, (npoints, nelements))
    elem_errs = create_dataset(element_flux, "errors", Float64, (npoints, nelements))
    write(element_flux, "header", string.(allelements))

    for i in eachindex(allelements)
        # Since each Macrostrat point has a corresponding EarthChem sample, use that 
        # sample to calculate flux of each element at each point
        elementflux = [sampleflux[j] * mbulk[allelements[i]][j] * 1e-2 for j = 1:npoints]

        # Save to file
        elementflux_val, elementflux_err = unmeasurementify(elementflux)
        elem_vals[:,i] += elementflux_val
        elem_errs[:,i] += elementflux_err
    end

    close(fid)


## --- Open the file
    # Redundant, but means we can skip straight to this cell and start analyzing things
    fid = h5open("$eroded_out", "r")

    # Total / bulk denundation at each point in kg/yr
    path = fid["bulk_denundation"]
    # bulk_denundation = read(path["values"]) .± read(path["errors"])
    bulk_denundation = read(path["values"])

    # Global flux of each element at each point in kg/yr
    path = fid["element_flux"]
    vals = read(path["values"])
    errs = read(path["errors"])
    header = Tuple(Symbol.(read(path["header"])))
    # elementflux = NamedTuple{header}([vals[:,i] .± errs[:,i] for i in eachindex(header)])
    elementflux = NamedTuple{header}([vals[:,i] for i in eachindex(header)])

    close(fid)

    # If I actually try to propagate uncertainty, it takes over 2 hours (which is not the 
    # run time, but the time where I gave up and killed the program).
    # Uncertainty comes from: standard deviation of slope at each point, erosion uncertainty
    # that's then propagated foward. We assume that wt.% of each element is exact.
    # For now, just don't... propagate. lol.


## --- Calculate absolute and relative flux
    # Conversion for kg to Gt. Data file is in units of kg/yr
    const kg_gt = 1e12

    # Global denundation in Gt/yr.
    global_denun = nansum(bulk_denundation ./ kg_gt)

    # Total global flux for each element in Gt/yr
    totalelementflux = NamedTuple{keys(elementflux)}([nansum(i) / kg_gt for i in elementflux])

    # Sum of all element fluxes in Gt/yr
    global_elem = nansum(totalelementflux)

    # Print to terminal
    @info """ Results:
    Total global denundation: $(round(global_denun, sigdigits=3)) Gt/yr
    Sum of all element fluxes: $(round(global_elem, sigdigits=3)) Gt/yr
    """
    # if global_elem.val > global_denun.val
    #     diff = global_elem.val - global_denun.val
    #     @warn """
    #     Mass of total eroded material from all elements is greater than bulk eroded mass. 
    #     Difference: $diff
    #     """
    # end
    if global_elem > global_denun
        diff = global_elem - global_denun
        @warn """
        Mass of total eroded material from all elements is greater than bulk eroded mass. 
        Difference: $diff
        """
    end


## --- Export crustal composition results
    # Preallocate (element row, rock type column)
    result = Array{Float64}(undef, length(elementflux), length(match_cats) + 1)

    # Calculate the absolute contribution of each rock type to the denudation of each
    # element. This is the sum of the individual contributions of each point mapped to that
    # rock type.

    # We use the rock classes assigned during lithological / geochemical sample matching, 
    # because this represents the geochemical composition of rocks of that class (e.g., if 
    # something is matched to basalt and rhyolite, and was assigned basalt for the purposes
    # of sample matching, we don't want to contaminate our rhyolite data with that point).
    # Our sample density is high enough that by the law of large numbers, we should get
    # back to the same result. 

    for i in eachindex(keys(match_cats))
        for j in eachindex(keys(elementflux))
            # result[j,i] = nansum(unmeasurementify(elementflux[j])[1][match_cats[i]])
            result[j,i] = nansum(elementflux[j][match_cats[i]])
        end
    end

    # Convert units from kg/yr to gt/yr
    result = result ./ kg_gt

    # Add global total denudation by element to the last row (units already converted)
    # global_element = unmeasurementify(totalelementflux)[1]
    global_element = values(totalelementflux)
    result[:,end] .= global_element

    # Save to file
    rows = string.(collect(keys(elementflux)))
    cols = string.(collect(keys(match_cats)))
    cols = hcat("", reshape(cols, 1, length(cols)), "total")

    writedlm(erodedabs_out, vcat(cols, hcat(rows, result)))                    # Absolute
    writedlm(erodedrel_out, vcat(cols, hcat(rows, result./global_element)))    # Relative / fraction

    # Terminal printout
    majors, minors = get_elements()
    nmajors = length(majors)
    majorcomp = round.(result[1:nmajors, end], digits=1)
    # majorcomp_rel = round.(result[1:nmajors, end]./global_denun.val*100, digits=1)
    majorcomp_rel = round.(result[1:nmajors, end]./global_denun*100, digits=1)

    @info """Composition of bulk eroded material:
    Absolute (gt/yr)
    $(join(keys(elementflux)[1:nmajors], " \t "))
    $(join(majorcomp, " \t "))

    Relative (%)
    $(join(keys(elementflux)[1:nmajors], " \t "))
    $(join(majorcomp_rel, " \t "))
    """

    
## --- End of File