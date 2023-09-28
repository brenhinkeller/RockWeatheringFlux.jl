# Translate the code from Keller et al., 2015 (https://doi.org/10.1038/nature14584) into 
# Julia, and run it for all rock types (not just igneous). 
# 
# This will give me something to compare SampleMatch results to.

# May run into some problems with minor sedimentary types, and it's unclear if sorting 
# metamorphic rocks by metased / metaign is doing anything useful

## --- Set up
    # Packages
    using StatGeochem
    using DelimitedFiles
    using Measurements
    using HDF5
    using MAT
    using ProgressMeter
    using LoopVectorization
    using Static

    # Local utilities
    include("utilities/Utilities.jl")


## --- Load EarthChem data 
    # Filtered and normalized to 100%
    fid = h5open("output/bulk.h5", "r")

    # Bulk
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Rock type matches
    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)


## --- Definitions we'll need for resampling
    # Major types are inclusive of all minor types
    minorsed, minorign, minormet = get_minor_types()
    for type in minorsed
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end
    for type in minormet
        bulk_cats.met .|= bulk_cats[type]
    end

    # Get the data we'll use in resampling. Exclude type, age, lat/lon
    bulkmatrix = float.(unelementify(bulk)[2:end,1:end-4])    # Also exclude header row
    header = collect(keys(bulk))[1:end-4]

    # Uncertainties based on standardized uncertainties used in Keller et al., 2015
    errVP = matread("data/volcanicplutonic/plutonic.mat")["plutonic"]["err"]
    errVP["In"] = 0.0       # Not listed, assume 0?

    uncert = Array{Float64}(undef, size(bulkmatrix))
    for i in eachindex(header)
        uncert[:,i] .= errVP[string(header[i])]
    end


## --- Resample based on spatial weights
    types = keys(bulk_cats)

    # Everyone's favorite file format!
    fid = h5open("output/resampled/resampled.h5", "w")
    g_main = create_group(fid, "vars")
    g = create_group(g_main, "data")

    for t in types
        t==:cover && continue
        println("Starting $t")

        # Calculate spatial weights, keep ∼1/5 of the data in each resampling
        k = invweight_location(bulk.Latitude[bulk_cats[t]], bulk.Longitude[bulk_cats[t]])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

        # Get all the data we'll use in the resampling
        rockdata = bulkmatrix[bulk_cats[t][:],:]
        rockuncert = uncert[bulk_cats[t][:],:]

        # Final dataset size should be proportional to the starting dataset, just because
        # we have 270K igneous samples and 6 evaporite samples. Minimum 500 samples
        nrows = max(count(bulk_cats[t]) * 5, 500)

        sim = bsresample(rockdata, rockuncert, nrows, p)

        # Save the data to the file, or whatever
        g₀ = create_group(g, "$t")
        g₀["data"] = sim
        g₀["k"] = k
    end

    # May as well just resample the everything too, as a treat
    k = invweight_location(bulk.Latitude, bulk.Longitude)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    sim = bsresample(bulkmatrix, uncert, 1_500_000, p)

    # And save
    g₀ = create_group(g, "bulk")
    g₀["data"] = sim
    g₀["k"] = k

    # Save the header too, we'll want that
    g_main["header"] = string.(header)

    close(fid)


## --- End of file