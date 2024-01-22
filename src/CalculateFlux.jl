## --- Setup
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5


## --- Load EarthChem data
    # Indices of matched EarthChem samples from SampleMatch.jl
    fid = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    # # Matched types, majors inclusive of minors
    # bulktype = string.(vec(fid[:,2]))
    # macro_cats = match_rocktype(bulktype[t])

    # Load bulk, but just the samples matched to the Macrostrat data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(fid)

    
## --- Load pre-generated Macrostrat data, but only if there's an associated EarthChem sample
    @info "Loading Macrostrat data"

    # Load and match
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
    )
    # header = read(fid["type"]["macro_cats_head"])
    # data = read(fid["type"]["macro_cats"])
    # data = @. data > 0
    # macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    
## --- Matched rock types
    # From the matched samples
    fid = readdlm("$matchedbulk_io")
    bulktype = string.(vec(fid[:,2]))
    out_cats = match_rocktype(bulktype[t])

    # From Macrostrat
    fid = h5open("$macrostrat_io", "r")
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    # Major types include minor types
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorsed
        out_cats.sed .|= out_cats[type]
        macro_cats.sed .|= macro_cats[type]
    end
    for type in minorvolc
        out_cats.volc .|= out_cats[type]
        macro_cats.volc .|= macro_cats[type]
    end
    for type in minorplut
        out_cats.plut .|= out_cats[type]
        macro_cats.plut .|= macro_cats[type]
    end
    for type in minorign
        out_cats.ign .|= out_cats[type]
        macro_cats.ign .|= macro_cats[type]
    end

    # # Figure out how many data points weren't matched
    # known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    # total_known = count(known_rocks)

    # matched = known_rocks .| macro_cats.cover
    # not_matched = .!matched
    # multi_matched = ((macro_cats.sed .& macro_cats.ign) .| (macro_cats.sed .& macro_cats.met) 
    #     .| (macro_cats.ign .& macro_cats.met)
    # )

    # # Print to terminal
    # @info """
    # Macrostrat parsing complete!
    # not matched = $(count(not_matched))
    # """

    # Delete cover
    macro_cats = delete_cover(macro_cats)
    out_cats = delete_cover(out_cats)
    

## --- Calculate erosion rate at each coordinate point of interest	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point with a known EarthChem sample
    # Modify this function to return an error as well
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
        elementflux = [sampleflux[j] * bulk[allelements[i]][j] * 1e-2 for j = 1:npoints]

        # Save to file
        elementflux_val, elementflux_err = unmeasurementify(elementflux)
        elem_vals[:,i] += elementflux_val
        elem_errs[:,i] += elementflux_err
    end

    close(fid)


## --- Open the file
    # Redundant, but means we can skip straight to analyzing things if we want
    fid = h5open("$eroded_out", "r")

    # Total / bulk denundation at each point in kg/yr
    path = fid["bulk_denundation"]
    bulk_denundation = read(path["values"]) .± read(path["errors"])

    # Global flux of each element at each point in kg/yr
    path = fid["element_flux"]
    vals = read(path["values"])
    errs = read(path["errors"])
    header = Tuple(Symbol.(read(path["header"])))
    elementflux = NamedTuple{header}([vals[:,i] .± errs[:,i] for i in eachindex(header)])

    close(fid)


## --- Calculate absolute and relative flux
    # Conversion for kg to Gt. Data file is in units of kg/yr
    const kg_gt = 1e12

    # Global denundation in Gt/yr
    global_denun = nansum(bulk_denundation) / kg_gt
    
    # Total global flux for each element in Gt/yr
    totalelementflux = NamedTuple{keys(elementflux)}([nansum(i) / kg_gt for i in elementflux])

    # Sum of all element fluxes in Gt/yr
    global_elem = nansum(totalelementflux)

    # Print to terminal
    @info """
    Total global denundation: $global_denun Gt/yr
    Sum of all element fluxes: $global_elem Gt/yr
    """
    if global_elem.val > global_denun.val
        diff = global_elem.val - global_denun.val
        @warn """
        Mass of total eroded material from all elements is greater than bulk eroded mass. 
        Difference: $diff
        """
    end


## --- Export crustal composition results
    # Preallocate (element row, rock type column)
    result = Array{Float64}(undef, length(elementflux), length(out_cats) + 1)

    # Calculate the absolute contribution of each rock type to the denudation of each
    # element. This is the sum of the individual contributions of each point mapped to that
    # rock type.

    # We use the randomly assigned, single-matched rock types, since the code is easier to
    # write. By the law of large numbers, this should be equivalent to using multi-matched
    # samples.

    for i in eachindex(keys(out_cats))
        for j in eachindex(keys(elementflux))
            result[j,i] = nansum(unmeasurementify(elementflux[j])[1][out_cats[i]])
        end
    end

    # Convert units from kg/yr to gt/yr
    result = result ./ kg_gt

    # Add global total denudation by element to the last row (units already converted)
    global_element = unmeasurementify(totalelementflux)[1]
    result[:,end] = global_element

    # Save to file
    rows = string.(collect(keys(elementflux)))
    cols = string.(collect(keys(out_cats)))
    cols = hcat("", reshape(cols, 1, length(cols)), "total")

    writedlm(erodedabs_out, vcat(cols, hcat(rows, result)))                    # Absolute
    writedlm(erodedrel_out, vcat(cols, hcat(rows, result./global_element)))    # Relative / fraction

    # Terminal printout
    majors, minors = get_elements()
    nmajors = length(majors)
    majorcomp = round.(result[1:nmajors, end], digits=1)
    majorcomp_rel = round.(result[1:nmajors, end]./global_denun.val*100, digits=1)

    @info """Composition of bulk eroded material:
    Absolute (gt/yr)
    $(join(keys(elementflux)[1:nmajors], " \t "))
    $(join(majorcomp, " \t "))

    Relative (%)
    $(join(keys(elementflux)[1:nmajors], " \t "))
    $(join(majorcomp_rel, " \t "))
    """

    
## --- End of File