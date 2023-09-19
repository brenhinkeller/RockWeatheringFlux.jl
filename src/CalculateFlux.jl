## --- Setup
    # Packages
    using StatGeochem
    using DelimitedFiles
    using ProgressMeter
    using HDF5
    using LoopVectorization
    using Static
    using Measurements

    # Local utilities
    include("utilities/Utilities.jl")


## --- Load EarthChem data
    # Indices of matched samples from SampleMatch.jl
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0     # Exclude samples with missing data

    # Load bulk, but just the samples matched to the Macrostrat data
    bulkfid = h5open("output/bulk.h5", "r")
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(bulkfid)

    
## --- Load pre-generated Macrostrat data, but only if there's an associated EarthChem sample
    @info "Loading Macrostrat data"

    # Load and match
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(macrofid["vars"]["rocklat"])[t],
        rocklon = read(macrofid["vars"]["rocklon"])[t],
        age = read(macrofid["vars"]["age"])[t],
    )
    header = read(macrofid["type"]["macro_cats_head"])
    data = read(macrofid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])

    # Major types include minor types
    minorsed, minorign, minormet = get_minor_types()
    for type in minorsed
        macro_cats.sed .|= macro_cats[type]
    end
    for type in minorign
        macro_cats.ign .|= macro_cats[type]
    end
    for type in minormet
        macro_cats.met .|= macro_cats[type]
    end

    # Figure out how many data points weren't matched
    known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    total_known = count(known_rocks)

    matched = known_rocks .| macro_cats.cover
    not_matched = .!matched
    multi_matched = ((macro_cats.sed .& macro_cats.ign) .| (macro_cats.sed .& macro_cats.met) 
        .| (macro_cats.ign .& macro_cats.met)
    )

    # Print to terminal
    @info """
    Macrostrat parsing complete!
    not matched = $(count(not_matched))
    """

    if count(multi_matched) > 0
        @warn """
        $(count(multi_matched)) conflicting matches present
        sed and ign = $(count(macro_cats.sed .& macro_cats.ign))
        sed and met = $(count(macro_cats.sed .& macro_cats.met))
        ign and met = $(count(macro_cats.ign .& macro_cats.met))
        """
    end
    

## --- Definitions
    # Elements of interest
    majors, minors = get_elements()
    allelements = [majors; minors]
    nelements = length(allelements)

    # Define rock sub-types (do not compute cover)
    npoints = count(t)
    subcats = collect(keys(macro_cats))
    deleteat!(subcats, findall(x->x==:cover,subcats))

    # Get relative proportion of all rock types
    nmajors = count(macro_cats.sed) + count(macro_cats.met) + count(macro_cats.ign)
    rel_proportion = NamedTuple{Tuple(subcats)}(
        float.([count(macro_cats[i]) for i in subcats]) / nmajors
    )
    

## --- Calculate erosion rate at each coordinate point of interest	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point with a known EarthChem sample
    # Modify this function to return an error as well
    rockslope = window_avg(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, halfwidth=5
    )

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)


## --- Calculate denundation at each point
    # Declare constants
    const crustal_density = 2750                                # kg/m³
    const unit_sample_area = (148940000 * 1000000) / npoints    # Area of contients / npoints (m²)

    # Create file to save data
    fid = h5open("$eroded_out", "w")

    # Denundation at each point
    sampleflux = [rock_ersn[i] * unit_sample_area * crustal_density * 1e-6 for i = 1:npoints]

    # Save to file
    sampleflux_val, sampleflux_err = unmeasurementify(sampleflux)
    bulk_denundation = create_group(fid, "bulk_denundation")
    write(bulk_denundation, "values", sampleflux_val)
    write(bulk_denundation, "errors", sampleflux_err)


## --- Calculate flux of each element at each point
    # Preallocate file space
    element_flux = create_group(fid, "element_flux")
    elem_vals = create_dataset(element_flux, "values", Float64, (npoints, nelements))
    elem_errs = create_dataset(element_flux, "errors", Float64, (npoints, nelements))
    write(element_flux, "header", string.(allelements))

    for i in eachindex(allelements)
        # Since each Macrostrat point has a corresponding EarthChem sample, use that 
        # sample to calculate flux of eachelement at each point
        elementflux = [sampleflux[j] * bulk[i][j] * 1e-2 for j = 1:npoints]

        # Save to file
        elementflux_val, elementflux_err = unmeasurementify(elementflux)
        elem_vals[:,i] += elementflux_val
        elem_errs[:,i] += elementflux_err
    end

    close(fid)


## --- Open, the file. Stop, having it be closed.
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
    globalflux = nansum(totalelementflux)

    # Print to terminal
    @info """
    Total global denundation: $global_denun Gt/yr
    Sum of all element fluxes: $globalflux Gt/yr
    """
    if globalflux.val > global_denun.val
        @warn """
        Mass of total eroded material from all elements is greater than bulk eroded mass. 
        """
    end


## --- Export crustal composition results
    # Preallocate
    result = zeros(length(header), length(macro_cats))  # Element row, rock type column
    rows = string.(header)
    cols = hcat("", string.(reshape(subcats, 1, length(subcats))), "global")

    # Absolute contribution of each rock type to element flux, calculated as the sum of
    # the individual contributions of each point
    p = Progress(npoints ÷ 100, desc="Calculating absolute element fluxes...")
    @time for i = 1:npoints
        # Get the rock types matched with the point, and classify as major / minor
        type = get_type(macro_cats, i, all_keys=true)
        mtype = Dict(zip((:sed, :ign, :met), zeros(Int, 3)))
        ntypes = 0

        for j in eachindex(type)
            t = class_up(type[j], minorsed, minorign, minormet)
            t===nothing && continue
            mtype[t] += 1
        end
        for j in keys(mtype)
            ntypes += ifelse(mtype[j]>1, mtype[j]-1, mtype[j])
        end
        
        # The contribution of the point should be split between the number of types matched
        # to that point
        for j in eachindex(header)
            contrib = (elementflux[header[j]][i] / ntypes).val
            contrib = ifelse(isnan(contrib), 0, contrib)

            # Get column number for the results table from subcats 
            col_n = zeros(Int, length(type))
            for t in eachindex(type)
                for s in eachindex(subcats)
                    if type[t]==subcats[s]
                        col_n[t]=s
                        break
                    end
                end
            end

            # Assign!
            t = col_n .!= 0
            for c in col_n[t]
                result[j,c] += contrib
            end
        end

        if i % 100 == 0
            next!(p)
        end
    end

    result = result ./ kg_gt

    # Add the global total to the last row
    global_by_element = unmeasurementify(totalelementflux)[1]
    result[:,end] = global_by_element

    # Absolute contribution
    writedlm("$erodedabs_out", vcat(cols, hcat(collect(rows), result)), ',')

    # Relative contribution
    # total_absolute = nansum(result[1:end,1:end-1], dims=2)
    # writedlm("$erodedrel_out", vcat(cols, hcat(collect(rows), result[:,1:end-1] ./ total_absolute, 
    #     result[:,end] ./ TEFval)), ','
    # )


## --- End of File
