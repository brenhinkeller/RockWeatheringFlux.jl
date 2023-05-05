# include("intermed.jl")

function llage1(bulkage::Vector, sampleage::Number,
        bulklat::Vector, bulklon::Vector, samplelat::Number, samplelon::Number,
        bulkgeochem::NamedTuple, samplegeochem::NamedTuple
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
            bulklat[i] = -samplelat
            bulklon[i] = -samplelon
        end
    end

    @turbo for i in 1:npoints
        # Age (σ = 38 Ma)
        ll_age[i] = -((bulkage[i] - sampleage)^2)/(38^2)

        # Distance (σ = 1.8 arc degrees)
        ll_dist[i] = -((haversine(samplelat, samplelon, bulklat[i], bulklon[i]))^2)/(1.8^2)
    end

    @. ll_total = ll_age + ll_dist

    # Reduce likelihoods to top quartile (75th percentile)
    bulktop = percentile(vec(ll_total), 75)
    @inbounds for i in 1:npoints
        if ll_total[i] < bulktop
            ll_total[i] = NaN
        end
    end
    
    # Geochemical log-likelihoods
    for elem in eachindex(bulkgeochem)
        @turbo for i in 1:npoints
            ll_total[i] += -((bulkgeochem[elem][i] - geochemdata[elem].m)^2)/(geochemdata[elem].e^2)
        end
    end

    matched_sample = rand_prop_liklihood(ll_total)
    return matched_sample
end

function llage2(bulkage::Vector, sampleage::Number,
        bulklat::Vector, bulklon::Vector, samplelat::Number, samplelon::Number,
        bulkgeochem::NamedTuple, samplegeochem::NamedTuple
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
            bulklat[i] = -samplelat
            bulklon[i] = -samplelon
        end
    end

    @turbo for i in 1:npoints
        # Age (σ = 38 Ma)
        ll_age[i] = -((bulkage[i] - sampleage)^2)/(38^2)

        # Distance (σ = 1.8 arc degrees)
        ll_dist[i] = -((haversine(samplelat, samplelon, bulklat[i], bulklon[i]))^2)/(1.8^2)
    end

    @. ll_total = ll_age + ll_dist
    
    # Geochemical log-likelihoods
    for elem in eachindex(bulkgeochem)
        @turbo for i in 1:npoints
            ll_total[i] += -((bulkgeochem[elem][i] - samplegeochem[elem].m)^2)/(samplegeochem[elem].e^2)
        end
    end

    matched_sample = rand_prop_liklihood(ll_total)
    return matched_sample
end

llage1(bulkage, sampleage, bulklat, bulklon, lat, lon, bulkgeochem, geochemdata)
llage2(bulkage, sampleage, bulklat, bulklon, lat, lon, bulkgeochem, geochemdata)

samplelat = lat
samplelon = lon
samplegeochem = geochemdata

# @benchmark llage1($bulkage, $sampleage, $bulklat, $bulklon, $lat, $lon, $bulkgeochem, $geochemdata)
# @benchmark llage2($bulkage, $sampleage, $bulklat, $bulklon, $lat, $lon, $bulkgeochem, $geochemdata)

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
    replace too low LLs with NaN 89010d08bf630c11a998b5dc2b95d5a588b772b7
=#
