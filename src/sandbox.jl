function geochem_likelihood(chem_samples, geochem_data)
    lh_geochem = zeros(Float64, count(chem_samples))
        for elem in eachindex(geochem_data)
            # Get all EarthChem samples for that element and this rock type
            # Replace NaNs with 0
            bulk_geochem = zeronan!(bulk[elem][chem_samples])

            # Assume Macrostrat sample is the average geochemistry for that rock type
            # TO DO: More granular estimate than just by rock type
            geochem_diff = abs.(bulk_geochem .- geochem_data[elem].m)
            
            # Calculate likelihood for that element and add to the other likelihoods
            lh_elem = -(geochem_diff.^2)./(geochem_data[elem].e^2)
            lh_geochem .+= lh_elem
        end

    return lh_geochem
end

function likelihood(chem_samples, sample_idxs, geochem_data, lat, lon, sample_age, bulklat, bulklon)
    matched_sample = Array{Int64}(undef, length(lat), 1)

    # For each Macrostrat sample, calculate
    for j in eachindex(lat)
        # Distance (σ = 1.8 arc degrees)
        dist = arcdistance(lat[j], lon[j], bulklat, bulklon)
        lh_dist = -(dist.^2)./(1.8^2)

        # Age (σ = 38 Ma)
        # TO DO: AgeEst vs Age? Maybe recalculate AgeEst?
        age_diff = abs.(bulk.AgeEst[chem_samples] .- sample_age[j])
        lh_age = -(age_diff.^2)./(38^2)

        # Geochemistry
        # TO DO: only look at ones that are close in age and space
        # lh_geochem = zeros(Float64, count(chem_samples))
        # for elem in eachindex(geochem_data)
        #     # Get all EarthChem samples for that element and this rock type
        #     # Replace NaNs with 0
        #     bulk_geochem = zeronan!(bulk[elem][chem_samples])

        #     # Assume Macrostrat sample is the average geochemistry for that rock type
        #     # TO DO: More granular estimate than just by rock type
        #     geochem_diff = abs.(bulk_geochem .- geochem_data[elem].m)
            
        #     # Calculate likelihood for that element and add to the other likelihoods
        #     lh_elem = -(geochem_diff.^2)./(geochem_data[elem].e^2)
        #     lh_geochem .+= lh_elem
        # end

        lh_geochem = geochem_likelihood(chem_samples, geochem_data)

        # Calculate total likelihood for each EarthChem sample
        # This has to be added in steps because nanadd can only do one array at a time
        lh_total = nanadd(lh_dist, lh_age)
        lh_total = nanadd(lh_total, lh_geochem)

        # Get the index of the most likely EarthChem sample
        # TO DO: the most likely sample has the largest likelihood? because if the 
            # difference is larger than the total likelihood value is more negative...
        # TO DO: take the average of some arbitrary percentile
        (val, idx) = findmax((lh_total[findall(!isnan, lh_total)]))
        matched_sample[j] = sample_idxs[findall(!isnan, lh_total)][idx]
    end

    return matched_sample
end



function lh_stripped(lat, lon, bulklat, bulklon)
    matched_sample = Array{Int64}(undef, length(lat), 1)

    for i in eachindex(lat)
        dist = arcdistance(lat[i], lon[i], bulklat, bulklon)
        lh_dist = -(dist.^2)./(1.8^2)

        lh_total = lh_dist
        (val, idx) = findmax((lh_total[findall(!isnan, lh_total)]))
        matched_sample[j] = sample_idxs[findall(!isnan, lh_total)][idx]
    end

    return lh_dist
end