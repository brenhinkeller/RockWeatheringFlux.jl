# TO DO: Define Measurement{Float64} <: Measurement{Number}

# New NaNStatistics methods for Measurements with NaN values. Converts all NaNs to 0

function NaNStatistics.nanadd(a::Measurement, b::Measurement)
    a = a.val*(a.val==a.val) ± a.err*(a.err==a.err)
    b = b.val*(b.val==b.val) ± b.err*(b.err==b.err)

    return a + b
end

function NaNStatistics.nansum(A::AbstractArray{Measurement{Float64}})
    Tₒ = Base.promote_op(+, eltype(A), Int)
    Σ = zero(Tₒ)
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        Aᵢ = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err)
        Σ += Aᵢ
    end
    return Σ
end

function NaNStatistics.nanmean(A::AbstractArray{Measurement{Float64}})
    Tₒ = Base.promote_op(/, eltype(A), Int)
    n = 0
    Σ = zero(Tₒ)
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ.val==Aᵢ.val || Aᵢ.err==Aᵢ.err
        Aᵢ = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err)
        n += notnan
        Σ += Aᵢ
    end
    return Σ / n
end


# NaNStatistics.nanstd(A; dims=:, dim=:, mean=nothing, corrected=true) = NaNStatistics.sqrt!(_nanvar(mean, corrected, A, dims, dim))

# function NaNStatistics.sqrt!(A::AbstractArray{Measurement{Float64}})
#     @inbounds for i ∈ eachindex(A)
#         A[i] = sqrt(A[i])
#     end
#     return A
# end

NaNStatistics.nanvar(A::AbstractArray{Measurement{Float64}}; dims=:, dim=:, mean=nothing, corrected=true) = 
    NaNStatistics._nanvar(mean, corrected, A::AbstractArray{Measurement{Float64}}, dims, dim)

function NaNStatistics._nanvar(::Nothing, corrected::Bool, A::AbstractArray{Measurement{Float64}}, ::Colon, ::Colon)
    Tₒ = Base.promote_op(/, eltype(A), Int)
    n = 0
    Σ = zero(Tₒ)
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ.val==Aᵢ.val || Aᵢ.err==Aᵢ.err
        Aᵢ = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err)
        n += notnan
        Σ += Aᵢ
    end
    μ = Σ / n

    σ² = ∅ = zero(typeof(μ))
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ.val==Aᵢ.val || Aᵢ.err==Aᵢ.err
        Aᵢ = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err)
        δ = Aᵢ - μ
        σ² += ifelse(notnan, δ * δ, ∅)
    end
    return σ² / max(n-corrected,0)
end

function NaNStatistics._nanvar(μ::Number, corrected::Bool, A::AbstractArray{Measurement{Float64}}, ::Colon, ::Colon)
    n = 0
    σ² = ∅ = zero(typeof(μ))
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ.val==Aᵢ.val || Aᵢ.err==Aᵢ.err
        Aᵢ = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err)
        δ = Aᵢ - μ
        n += notnan
        σ² += ifelse(notnan, δ * δ, ∅)
    end
    return σ² / max(n-corrected, 0)
end


