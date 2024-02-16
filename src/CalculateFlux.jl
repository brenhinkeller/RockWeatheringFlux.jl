## --- Setup
    # Calculate the flux of each element at each point, as well as the composition of 
    # total produced sediment (bulk and by lithologic class).

    # Optionally run lines 12-54, and skip to line 122 to load already calculated data.

    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, Dates
    using StatsBase


## --- Load data
    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
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

    # Get slope at each coordinate point
    # Function returns the standard deviation of slope in each window, which we don't
    # actually care about propagating
    rockslope = movingwindow(srtm15_slope, rocklat, rocklon, srtm15_sf, n=5)
    rockslope = Measurements.value.(rockslope)

    # Calculate all erosion rates (mm/kyr) (propagate uncertainty)
    rock_ersn = emmkyr.(rockslope);


## --- Calculate denundation at each point
    # Area of land:
        # Total: 149,733,926
        # Antarctica: 14,200,000 
        # Greenland: 2,166,086
        # Area in this study = Total - (Antarctica + Greenland) = 133_367_840
    # Declare constants
    const crustal_density = 2750                            # kg/m³
    npoints = length(mbulk.SiO2)                            # Number of *matched* points
    unit_sample_area = (133_367_840 * 1000000) / npoints    # Area of continents / npoints (m²)

    # Create file to save data
    fid = h5open(eroded_out, "w")

    # Denundation at each point (kg/yr), for 
    sampleflux = [rock_ersn[i] * unit_sample_area * crustal_density * 1e-6 for i = 1:npoints];

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
    @info """ Results ± 2σ s.d.:
    Total global denundation: $(round(global_denun_val, sigdigits=3)) ± $(round(global_denun_err, sigdigits=3)*2) Gt/yr
    Sum of all element fluxes: $(round(global_elem_val, sigdigits=3)) ± $(round(global_elem_err, sigdigits=3)*2) Gt/yr
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
    cols = hcat("elem", reshape(string.(collect(keys(class))), 1, :))

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
    writedlm(erodedrel_out_err, vcat(cols, hcat(rows, result_err./result[:,end])))


## --- Calculate and export the composition of eroded material
    # How many samples explain 90% of the matches?
    npoints = unique_sample(mbulk.Sample_ID, 90)

    # Compute the composition of eroded material [wt.%] for each lithologic class as the
    # fraction of that element out of the total erosion for that class.
    for i in eachindex(keys(class))
        normconst = result[end,i]
        result[:,i] ./= normconst
        result_err[:,i] ./= normconst
    end
    writedlm(erodedcomp_out, vcat(cols, hcat(rows, result.*100)))
    writedlm(erodedcomp_out_err, vcat(cols, hcat(rows, result_err.*100)))

    # Terminal printout
    majorcomp_rel = round.(result[1:nmajors, end].*100, sigdigits=3)
    majorcomp_rel_err = round.(result_err[1:nmajors, end].*100 ./ sqrt(npoints) .*2, sigdigits=1)

    @info """Composition of bulk eroded material:
    Absolute [Gt/yr] ± 1σ s.d.
      $(join(keys(elementflux_val)[1:nmajors], " \t "))
      $(join(majorcomp, " \t "))
    ± $(join(majorcomp_err, " \t "))

    Composition [wt.%] ± 2 s.e.
      $(join(keys(elementflux_val)[1:nmajors], " \t "))
      $(join(rpad.(majorcomp_rel, 8), " "))
    ± $(join(rpad.(majorcomp_rel_err, 8), " "))
    """

    
## --- End of File