function ll2(lat::Number, lon::Number, bulklat, bulklon, sampleidxs, sampleage::Number, bulkage, geochemdata, bulkgeochem)
    # TO DO: reuse variable names, reduce allocations
    # TO DO: for missing ages / locations, penalize but do not exclude. age - negative sample age?
    
    # Preallocate
    lh_age = Array{Float64}(undef, length(bulklat), 1)
    lh_dist = Array{Float64}(undef, length(bulklat), 1)
    lh_total = Array{Float64}(undef, length(bulklat), 1)

    # Find most likely EarthChem sample
    
    # Age (σ = 38 Ma)
    age_diff = abs.(bulkage .- sampleage)
    lh_age .= -(age_diff.^2)./(38^2)

    # Distance (σ = 1.8 arc degrees)
    dist = arcdistance(lat, lon, bulklat, bulklon)
    lh_dist .= -(dist.^2)./(1.8^2)

    # Find current likelihoods, only look at geochemical data if the samples are close
    lh_total .= vec(nanadd.(lh_dist, lh_age))
    reduced = percentile(vec(lh_total), 25)
    reduced_idx = findall(>=(reduced), lh_total)
    lh_total = lh_total[reduced_idx]

    # Find most geochemically similar rock from this reduced data set
    lh_geochem = zeros(Float64, length(lh_total))
    @inbounds for elem in eachindex(geochemdata)
        # Get reduced EarthChem data for this element and this rock type
        reduced_bulkgeochem = bulkgeochem[elem][reduced_idx]

        # NaN -> 0 means missing data is not excluded, but it is penalized
        geochem_elem = zeronan!(reduced_bulkgeochem)

        # Geochemical distance (σ = geochemdata[elem].e)
        geochem_diff = abs.(geochem_elem .- geochemdata[elem].m)
        lh_geochem .+= -(geochem_diff.^2)./(geochemdata[elem].e^2)
    end

    # Select a sample, weighted based on likelihood
    matched_sample = rand_prop_liklihood(lh_total)

    return matched_sample
end

function llage1(bulkage::Vector, sampleage::Number)
    # Preallocate
    npoints = length(bulkage)
    lh_age = Array{Float64}(undef, npoints, 1)

    # Age (σ = 38 Ma)
    @. age_diff = bulkage - sampleage
    @. lh_age = -(age_diff^2)/(38^2)

    matched_sample = rand_prop_liklihood(lh_age)
    return matched_sample
end

function llage2(bulkage::Vector, sampleage::Number)
    # Preallocate
    npoints = length(bulkage)
    lh_age = Array{Float64}(undef, npoints, 1)
    lh_age = Array{Float64}(undef, npoints, 1)

    chunks = vcat(collect(1:50:npoints), npoints+1)

    for i in 1:length(chunks[1:end-1])
        j = chunks[i]
        k = chunks[i+1]-1
        @views agechunk = bulkage[j:k]
        
        # Age (σ = 38 Ma)
        # age_diff = agechunk .- sampleage
        @. lh_age[j:k] = -((agechunk .- sampleage)^2)/(38^2)
    end

    matched_sample = rand_prop_liklihood(lh_age)
    return matched_sample
end

llage1(bulkage, sampleage)
llage2(bulkage, sampleage)