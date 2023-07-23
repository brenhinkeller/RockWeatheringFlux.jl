## --- Setup
    # Packages
    using StatGeochem
    using DelimitedFiles
    using ProgressMeter
    using Measurements
    using HDF5
    using LoopVectorization
    using Static

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
        rocklat = read(macrofid["rocklat"]),
        rocklon = read(macrofid["rocklon"]),
        typecategory = read(macrofid["typecategory"]),
    )
    close(macrofid)
    macro_cats = match_rocktype(macrostrat.typecategory[t])

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
    

## --- Calculate erosion rate at each coordinate point of interest	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point with a known EarthChem sample
    # Modify this function to return an error as well
    rockslope = avg_over_area(srtm15_slope, macrostrat.rocklat[t], macrostrat.rocklon[t], 
        srtm15_sf, halfwidth=7
    )

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)


## --- Get erosion and continental area for each rock type
    # Average erosion by rock type (m/Myr)
    erosion = NamedTuple{Tuple(subcats)}(
        [nanmean(rock_ersn[macro_cats[i]]) for i in subcats]
    )

    # Crustal area (m²), assume proportional distribution of rock types under cover
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    crustal_area = NamedTuple{Tuple(subcats)}(
        [count(macro_cats[i]) / total_known * contl_area for i in subcats]
    )


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
        # sample to calculate flux
        elementflux = [sampleflux[j] * bulk[i][j] * 1e-2 for j = 1:npoints]

        # Save to file
        elementflux_val, elementflux_err = unmeasurementify(elementflux)
        elem_vals[:,i] += elementflux_val
        elem_errs[:,i] += elementflux_err
    end

    close(fid)

## --- End of File
