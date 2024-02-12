## --- Setup
    # Calculate the flux of each element at each point, as well as the composition of 
    # total produced sediment (bulk and by lithologic class).

    # Optionally run lines 12-54, and skip to line 122 to load already calculated data.

    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, Dates


## --- Load data
    # Indices of matched samples from SampleMatch.jl
    fid = readdlm(matchedbulk_io)
    matches = Int.(vec(fid[:,1]))
    t = @. matches != 0

    # Rock class of each sample
    # We use the rock classes assigned during lithological / geochemical sample matching, 
    # because this represents the geochemical composition of rocks of that class (e.g., if 
    # something is matched to basalt and rhyolite, and was assigned basalt for the purposes
    # of sample matching, we don't want to contaminate our rhyolite data with that point).
    # Our sample density is high enough that by the law of large numbers, we should get
    # back to the same result. 
    match_cats = match_rocktype(string.(vec(fid[:,2]))[t]);
    include_minor!(match_cats)
    match_cats = delete_cover(match_cats)


    # Geochemical data for each sample
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][matches[t]] for i in eachindex(header)])
    close(fid)

    # Location of each sample
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
    )
    header = Tuple(Symbol.(read(fid["type"]["macro_cats_head"])))
    data = read(fid["type"]["metamorphic_cats"])
    data = @. data > 0
    metamorphic_cats = NamedTuple{header}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    # Metamorphic samples
    include_minor!(metamorphic_cats);
    match_cats.met .|= (metamorphic_cats.sed .| metamorphic_cats.ign .| metamorphic_cats.met)


## --- Calculate erosion rate at each coordinate point of interest	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point
    # Function returns the standard deviation of slope in each window, which we don't
    # actually care about propagating
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )
    rockslope = Measurements.value.(rockslope)

    # Calculate all erosion rates (mm/kyr) (propagate uncertainty)
    rock_ersn = emmkyr.(rockslope)


## --- Calculate denundation at each point
    # Area of land:
        # Total: 149,733,926
        # Antarctica: 14,200,000 
        # Greenland: 2,166,086
        # Area in this study = Total - (Antarctica + Greenland) = 133_367_840
    # Declare constants
    const npoints = length(macrostrat.rocklat)                # Number of points
    const crustal_density = 2750                              # kg/m³
    const unit_sample_area = (133367840 * 1000000) / npoints  # Area of continents / npoints (m²)

    # Create file to save data
    fid = h5open("$eroded_out", "w")

    # Denundation at each point (kg/yr)
    sampleflux = [rock_ersn[i] * unit_sample_area * crustal_density * 1e-6 for i = 1:npoints]

    # Save to file
    sampleflux_val, sampleflux_err = unmeasurementify(sampleflux)
    bulk_denudation = create_group(fid, "bulk_denudation")
    write(bulk_denudation, "values", sampleflux_val)
    write(bulk_denudation, "errors", sampleflux_err)


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
    path = fid["bulk_denudation"]
    bulk_denudation_val = read(path["values"])
    bulk_denudation_err = read(path["errors"])

    # Global flux of each element at each point in kg/yr
    path = fid["element_flux"]
    vals = read(path["values"])
    errs = read(path["errors"])
    header = Tuple(Symbol.(read(path["header"])))
    elementflux_val = NamedTuple{header}([vals[:,i] for i in eachindex(header)])
    elementflux_err = NamedTuple{header}([errs[:,i] for i in eachindex(header)])

    close(fid)


## --- Calculate mass of bulk sediment, and mass of sediment by element
    # Conversion for kg to Gt. Data file will be in units of kg/yr
    const kg_gt = 1e12

    # Global denundation in Gt/yr.
    global_denun_val = nansum(bulk_denudation_val ./ kg_gt)
    global_denun_err = sqrt(nansum((bulk_denudation_err ./ kg_gt).^2))

    # Total global flux for each element in Gt/yr
    totalelementflux_val = NamedTuple{keys(elementflux_val)}([nansum(i) / kg_gt for i in elementflux_val])
    totalelementflux_err = NamedTuple{keys(elementflux_err)}([sqrt(nansum((i ./ kg_gt).^2)) 
        for i in elementflux_err]
    )

    # Sum of all element fluxes in Gt/yr
    global_elem_val = nansum(values(totalelementflux_val))
    global_elem_err = sqrt(nansum(values(totalelementflux_err).^2))

    # Print to terminal
    @info """ Results:
    Total global denundation: $(round(global_denun_val, sigdigits=3)) ± $(round(global_denun_err, sigdigits=3)) Gt/yr
    Sum of all element fluxes: $(round(global_elem_val, sigdigits=3)) ± $(round(global_elem_err, sigdigits=3)) Gt/yr
    """

    if global_elem_val > global_denun_val
        diff = global_elem_val - global_denun_val
        @warn """
        Mass of total eroded material from all elements is greater than bulk eroded mass. 
        Difference: $diff
        """
    end


## --- Export absolute and relative contribution to eroded material
    # Get major elements for terminal printouts
    majors = get_elements()[1]
    nmajors = length(majors)

    # We want an option to filter for all samples and undifferentiated metamorphics 
    class = merge(match_cats, (bulk=trues(length(match_cats[1])),))

    # Preallocate (element row, rock type column)
    result = Array{Float64}(undef, length(elementflux_val)+1, length(class))
    result_err = similar(result)
    rows = vcat(string.(collect(keys(elementflux_val))), "Total")
    cols = hcat("", reshape(string.(collect(keys(class))), 1, :))

    # Calculate the absolute contribution of each rock class to the denudation of each
    # element. This is the sum of the individual contributions of each spatial point
    for i in eachindex(keys(class))
        for j in eachindex(keys(elementflux_val))
            result[j,i] = nansum(elementflux_val[j][class[i]])
            result_err[j,i] = sqrt(nansum((elementflux_err[j][class[i]]).^2))
        end
    end
    result[end,:] .= vec(nansum(result[1:end-1,:], dims=1))
    result_err[end,:] .= vec(sqrt.(nansum((result_err[1:end-1,:]).^2, dims=1)))

    # Convert units from kg/yr to gt/yr.
    # We now have absolute amount of material eroded by lithologic class each year. 
    result ./= kg_gt
    result_err ./= kg_gt
    writedlm(erodedabs_out, vcat(cols, hcat(rows, result)))
    writedlm(erodedabs_out_err, vcat(cols, hcat(rows, result_err)))

    # Get absolute values for printout
    majorcomp = round.(result[1:nmajors, end], digits=1)
    majorcomp_err = round.(result_err[1:nmajors, end], digits=1)

    # Also calculate the fraction of total erosion each class contributes. Divide the
    # absolute contribution by the total amount of each element that erodes each year.
    # Note that major classes will be the sum of their constituent minor classes.
    writedlm(erodedrel_out, vcat(cols, hcat(rows, result./result[:,end])))
    writedlm(erodedrel_out_err, vcat(cols, hcat(rows, result_err./result_err[:,end])))


## --- Calculate and export the composition of eroded material
    # Compute the composition of eroded material [wt.%] for each lithologic class as the
    # fraction of that element out of the total erosion for that class.
    for i in eachindex(keys(class))
        result[:,i] ./= result[end,i]
        result_err[:,i] .= sqrt.((result_err[:,i] ./ result[:,i]).^2)
    end
    writedlm(erodedcomp_out, vcat(cols, hcat(rows, result.*100)))
    writedlm(erodedcomp_out_err, vcat(cols, hcat(rows, result_err.*100)))

    # Terminal printout
    majorcomp_rel = round.(result[1:nmajors, end].*100, digits=1)
    majorcomp_rel_err = round.(result_err[1:nmajors, end].*100, digits=1)

    @info """Composition of bulk eroded material:
    Absolute [Gt/yr]
      $(join(keys(elementflux_val)[1:nmajors], " \t "))
      $(join(majorcomp, " \t "))
    ± $(join(majorcomp_err, " \t "))

    Composition [wt.%]
      $(join(keys(elementflux_val)[1:nmajors], " \t "))
      $(join(majorcomp_rel, " \t "))
    ± $(join(majorcomp_rel_err, " \t "))
    """

    
## --- End of File