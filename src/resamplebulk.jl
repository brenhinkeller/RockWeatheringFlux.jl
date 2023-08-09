# Fuck it! Resample every rock type! Brute force your way through it!
#
# This is based **only** on spatial distribution---I want current crustal composition, not
# something I can use as a timeseries

## --- Set up
    # Packages
    using StatGeochem
    using DelimitedFiles
    using Measurements
    using HDF5
    using ProgressMeter
    using LoopVectorization
    using Static

    # Local utilities
    include("utilities/Utilities.jl")

    # Get bulk data (normalized to 100%)
    bulkfid = h5open("output/bulk.h5", "r")
        header = read(bulkfid["bulk"]["header"])
        data = read(bulkfid["bulk"]["data"])
        bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(bulkfid)
    bulk_cats = match_earthchem(bulk.Type, major=false)


## --- Do stuff (resample)
    # Types to resample. Resample alluvium too, why not
    sed, ign, met = get_minor_types()
    types = (sed..., :sed, ign..., :ign, met..., :met, :alluvium)

    for t in types
        # Terminal printout
        println("Saving $t")

        if length(bulk.Latitude[bulk_cats[t]])==0
            @warn "No samples for $t type"
            continue
        end

        # Calculate resampling weights based on spatial distribution only
        # I'm a fool! There's a specific function just for location
        k = invweight_location(bulk.Latitude[bulk_cats[t]], bulk.Longitude[bulk_cats[t]])

        # Keep ~ 1/5 of the data in each resampling
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        s = @. !isnan(p)

        # Assume error ± 0 wt.%
        nrows = 1_000_000
        e = zeros(count(s))
        silica = bsresample(bulk.SiO2[bulk_cats[t]][s], zeros(length(bulk.SiO2[bulk_cats[t]][s])),
            nrows, p[s]
        )         

        # Save data to file
        gap = length(silica) - length(k)
        k = vcat(k, fill(NaN, gap))
        p = vcat(p, fill(NaN, gap))
        writedlm("output/resampled/rs_$t.tsv", vcat(["k" "p" "SiO2"], hcat(k, p, silica)))
    end


## --- End of file