function lh_stripped(lat, lon, bulklat, bulklon, sample_idxs, sample_age, bulkage, geochem_data, bulkgeochem)
    matched_sample = Array{Int64}(undef, length(lat), 1)

    # Find most likely EarthChem sample for each Macrostrat sample
    for i in eachindex(lat)
        # Age (σ = 38 Ma)
        age_diff = abs.(bulkage .- sample_age[i])
        lh_age = -(age_diff.^2)./(38^2)

        # Distance (σ = 1.8 arc degrees)
        dist = arcdistance(lat[i], lon[i], bulklat, bulklon)
        lh_dist = -(dist.^2)./(1.8^2)

        # Subtract vectors
        lh_total = nanadd.(lh_dist, -lh_age)

        # Only look at geochemical data if the samples are reasonably close
        reduced = percentile(lh_total, 25)
        reduced_idx = findall(>=(reduced), lh_total)
        lh_total = lh_total[reduced_idx]

        # Find most geochemically similar rock from this reduced data set
        # lh_geochem = zeros(Float64, length(lh_total))
        # for elem in eachindex(geochem_data)
        #     # # Get all EarthChem data for that rock type, change NaNs to 0s
        #     # # NaN -> 0 means missing data is not excluded, but it is penalized
        #     # geochem_elem = zeronan!(bulkgeochem[elem])

        #     # # Geochemical distance (σ = geochem_data[elem].e)
        #     # geochem_diff = abs.(geochem_elem .- geochem_data[elem].m)
        #     # lh_geochem .-= -(geochem_diff.^2)./(geochem_data[elem].e^2)
        # end

        # Find least negative / most likely
        (val, idx) = findmax((lh_total[findall(!isnan, lh_total)]))
        matched_sample[i] = sample_idxs[reduced_idx][findall(!isnan, lh_total)][idx]
    end

    return matched_sample
end

@code_warntype lh_stripped(lat, lon, bulklat, bulklon, sample_idxs, sample_age, bulkage, geochem_data, bulkgeochem)