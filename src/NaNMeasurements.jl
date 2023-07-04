# New NaNStatistics methods for Measurements with NaN values

# Ignore any NaN ± NaN values, and respectively convert val ± NaN and NaN ± err to val ± 0
# and 0 ± err.

using Static

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

NaNStatistics.nanstd(A; dims=:, dim=:, mean=nothing, corrected=true) = 
    NaNStatistics.sqrt!(NaNStatistics._nanvar(mean, corrected, A, dims, dim))

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

"""
```julia
    zeronan!(A::AbstractArray{Measurement{Float64}}, allnans=true)
```

Replace all `NaN`s in A with zeros of the same type.

Any element containing a `NaN` in _either_ the value or error will be replaced with 
0.0 ± 0.0. Optionally specify `allnans` as `false` to replace the `NaN` value with a zero
while retaining the non-`NaN` value in the value or error. 

## Examples
```julia
julia> A = [2.0 ± NaN]
1-element Vector{Measurement{Float64}}:
 2.0 ± NaN

julia> zeronan!(A)
1-element Vector{Measurement{Float64}}:
 0.0 ± 0.0

julia> A = [2.0 ± NaN]
1-element Vector{Measurement{Float64}}:
 2.0 ± NaN

julia> zeronan!(A, allnans=false)
1-element Vector{Measurement{Float64}}:
 2.0 ± 0.0
```
"""
NaNStatistics.zeronan!(A, allnans::Bool) = NaNStatistics.zeronan!(A, static(allnans))

NaNStatistics.zeronan!(A, allnans::True) = NaNStatistics.zeronan!(A)

function NaNStatistics.zeronan!(A::AbstractArray{Measurement{Float64}}, allnans::False)
    ∅ = zero(eltype(A))
    @inbounds for i ∈ eachindex(A)
        Aᵢ = A[i]
        notnan = Aᵢ.val==Aᵢ.val && Aᵢ.err==Aᵢ.err
        if !notnan
            A[i] = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err)
        end
    end
    return A
end


## --- End of file