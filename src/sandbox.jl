# include("intermed.jl")

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

function llage1(bulkage::Vector, sampleage::Number,
        bulklat::Vector, bulklon::Vector, lat::Number, lon::Number
    )
    # Preallocate
    npoints = length(bulkage)
    ll_age = Array{Float64}(undef, npoints, 1)
    ll_dist = Array{Float64}(undef, npoints, 1)
    ll_total = Array{Float64}(undef, npoints, 1)

    # Replace missing values: this will penalize but not exclude missing data
    @inbounds for i in 1:npoints
        if isnan(bulkage[i])
            bulkage[i] = -sampleage
        end

        # Assume if one coordinate is missing, so is the other one
        if isnan(bulklat[i])
            bulklat[i] = -lat
            bulklon[i] = -lon
        end
    end

    @turbo for i in 1:npoints
        # Age (σ = 38 Ma)
        ll_age[i] = -((bulkage[i] - sampleage)^2)/(38^2)

        # Distance (σ = 1.8 arc degrees)
        ll_dist[i] = -((haversine(lat, lon, bulklat[i], bulklon[i]))^2)/(1.8^2)
    end

    @. ll_total = ll_age + ll_dist

    bulktop = percentile(vec(ll_total), 25)
    reduced_idx = findall(>=(bulktop), ll_total)
    ll_total = ll_total[reduced_idx]

    matched_sample = rand_prop_liklihood(ll_total)
    return matched_sample
end

function llage2(bulkage::Vector, sampleage::Number,
        bulklat::Vector, bulklon::Vector, lat::Number, lon::Number
    )
    # Preallocate
    npoints = length(bulkage)
    ll_age = Array{Float64}(undef, npoints, 1)
    ll_dist = Array{Float64}(undef, npoints, 1)
    ll_total = Array{Float64}(undef, npoints, 1)

    # Replace missing values: this will penalize but not exclude missing data
    @inbounds for i in 1:npoints
        if isnan(bulkage[i])
            bulkage[i] = -sampleage
        end

        # Assume if one coordinate is missing, so is the other one
        if isnan(bulklat[i])
            bulklat[i] = -lat
            bulklon[i] = -lon
        end
    end

    @turbo for i in 1:npoints
        # Age (σ = 38 Ma)
        ll_age[i] = -((bulkage[i] - sampleage)^2)/(38^2)

        # Distance (σ = 1.8 arc degrees)
        ll_dist[i] = -((haversine(lat, lon, bulklat[i], bulklon[i]))^2)/(1.8^2)
    end

    @. ll_total = ll_age + ll_dist

    # Reduce likelihoods to top quartile (75th percentile)
    bulktop = percentile(vec(ll_total), 75)
    @inbounds for i in 1:npoints
        if ll_total[i] < bulktop
            ll_total[i] = NaN
        end
    end
    

    matched_sample = rand_prop_liklihood(ll_total)
    return matched_sample
end

llage1(bulkage, sampleage, bulklat, bulklon, lat, lon)
llage2(bulkage, sampleage, bulklat, bulklon, lat, lon)


# @benchmark llage1($bulkage, $sampleage, $bulklat, $bulklon, $lat, $lon)
# @benchmark llage2($bulkage, $sampleage, $bulklat, $bulklon, $lat, $lon)

#=
to try:
    NEW current project don't remove missing data, just penalize
    reduce lh_total size. test (for V2) if I can directly assign things to lh_total without the age and dist intermeds
    two chunk loops instead of 1
    allocating the correct size in the array
=#

# git log -1 --pretty=format:"%s, %H"
#=
Tried and didn't work:
    Preallocating distance array in V1 cfda79de984a7d4a83332608e2b4dfc71a032cfb
    Use smaller chunk size in V2 2b75bdc3b4afb97066f0849c79578ac4673f8d12
=#

#=
Things that DID work
    looping through bulk one at a time instead of chunks (V2) 860f8d8d2fb011c753046ba0a89ef0f5fe23c891
    @inbounds and rm temp variables (I think?) (V2) a140d88d17a5fc9ff24c94c4f9e6545c8a0aeae2
    @turbo (3x speedup) d144f06055c2f5cc45503b4af47187270981fd42
    remove NaNs 3b6270780eb6a162312346c4a92b4cf7e28fbd79
=#
