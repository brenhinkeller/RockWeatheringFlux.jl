# New NaNStatistics methods for Measurements with NaN values

# Ignore any NaN ± NaN and NaN ± err values. Convert val ± NaN to val ± 0

using Static

"""
```julia
    unmeasurementify(A::AbstractArray{Measurement{Float64}})
```

Separate an array `A` of `measurements` into an array of values and an array of errors.

## Example
```julia
julia> A = [1 ± 0.1, 2 ± 0.1]
2-element Vector{Measurement{Float64}}:
 1.0 ± 0.1
 2.0 ± 0.1

julia> val, err = unmeasurementify(A)
([1.0, 2.0], [0.1, 0.1])
```
"""
function unmeasurementify(A::AbstractArray{Measurement{Float64}})
    val = fill(NaN, length(A))
    err = fill(NaN, length(A))
    for i in eachindex(A)
        val[i] = A[i].val
        err[i] = A[i].err
    end
    return val, err
end

function NaNStatistics.nanadd(a::Measurement, b::Measurement)
    a = a.val*(a.val==a.val) ± a.err*(a.err==a.err)
    b = b.val*(b.val==b.val) ± b.err*(b.err==b.err)

    return a + b
end

function NaNStatistics.nansum(A::AbstractArray{Measurement{Float64}})
    Aᵥ, Aₑ = unmeasurementify(A)
    Tₒ = Base.promote_op(+, eltype(Aᵥ), Int)
    Σᵥ = Σₑ = ∅ = zero(Tₒ)

    @inbounds for i ∈ eachindex(A)
        Avalᵢ = Aᵥ[i]
        notnanval = Avalᵢ==Avalᵢ
        Σᵥ += ifelse(notnanval, Avalᵢ, ∅)

        Aerrᵢ = Aₑ[i]
        notnanerr = Aerrᵢ==Aerrᵢ
        Σₑ += ifelse(notnanerr, Aerrᵢ * Aerrᵢ, ∅)
    end
    return Σᵥ ± sqrt(sum(Σₑ))
end

function NaNStatistics.nanmean(A::AbstractArray{Measurement{Float64}})
    Aᵥ, Aₑ = unmeasurementify(A)
    Tₒ = Base.promote_op(/, eltype(Aᵥ), Int)
    n = 0
    Σᵥ = Σₑ = ∅ = zero(Tₒ)
    
    @inbounds for i ∈ eachindex(Aᵥ)
        Avalᵢ = Aᵥ[i]
        notnanval = Avalᵢ==Avalᵢ
        Σᵥ += ifelse(notnanval, Avalᵢ, ∅)

        Aerrᵢ = Aₑ[i]
        notnanerr = Aerrᵢ==Aerrᵢ
        Σₑ += ifelse(notnanerr, Aerrᵢ * Aerrᵢ, ∅)

        n += (notnanval || notnanerr)
    end
    return (Σᵥ / n) ± (sqrt(sum(Σₑ)) / n)
end

NaNStatistics.nanstd(A::AbstractArray{Measurement{Float64}}; mean=nothing, corrected=true) = 
    NaNStatistics.sqrt!(NaNStatistics._nanvar(mean, corrected, A::AbstractArray{Measurement{Float64}}))

NaNStatistics.nanvar(A::AbstractArray{Measurement{Float64}}; mean=nothing, corrected=true) = 
    NaNStatistics._nanvar(mean, corrected, A::AbstractArray{Measurement{Float64}})

function NaNStatistics._nanvar(::Nothing, corrected::Bool, A::AbstractArray{T}) where T <: Measurement{Float64}
    @warn "Unpredictable behavior in nanvar when errors are 0"

    Aᵥ, Aₑ = unmeasurementify(A)
    Tₒ = Base.promote_op(/, eltype(Aᵥ), Int)
    n = 0
    Σᵥ = Σₑ = ∅ = zero(Tₒ)
    @inbounds for i ∈ eachindex(A)
        Avalᵢ = Aᵥ[i]
        notnanval = Avalᵢ==Avalᵢ
        Σᵥ += ifelse(notnanval, Avalᵢ, ∅)

        Aerrᵢ = Aₑ[i]
        notnanerr = Aerrᵢ==Aerrᵢ
        Σₑ += ifelse(notnanerr, Aerrᵢ * Aerrᵢ, ∅)

        n += (notnanval || notnanerr)
    end
    μᵥ = Σᵥ / n
    μₑ = sqrt(sum(Σₑ)) / n

    σ²ᵥ = σ²ₑ = ∅ = zero(typeof(μᵥ))
    @inbounds for i ∈ eachindex(A)
        δᵥ = Aᵥ[i] - μᵥ
        notnanval = δᵥ==δᵥ
        σ²ᵥ += ifelse(notnanval, δᵥ * δᵥ, ∅)

        δₑ = (δᵥ * δᵥ) * sqrt(4 * (Aₑ[i] * Aₑ[i] + μₑ * μₑ) / (δᵥ * δᵥ))
        notnanerr = δₑ==δₑ
        σ²ₑ += ifelse(notnanerr, δₑ * δₑ, ∅)
    end
    return (σ²ᵥ / max(n-corrected,0)) ± (sqrt(sum(σ²ₑ)) / max(n-corrected,0))
end

function NaNStatistics._nanvar(μ::Number, corrected::Bool, A::AbstractArray{T}) where T <: Measurement{Float64}
    @warn "Unpredictable behavior in nanvar when errors are 0"

    Aᵥ, Aₑ = unmeasurementify(A)
    μᵥ = μ.val
    μₑ = μ.err
    σ²ᵥ = σ²ₑ = ∅ = zero(typeof(μᵥ))

    @inbounds for i ∈ eachindex(A)
        δᵥ = Aᵥ[i] - μᵥ
        notnanval = δᵥ==δᵥ
        σ²ᵥ += ifelse(notnanval, δᵥ * δᵥ, ∅)

        δₑ = (δᵥ * δᵥ) * sqrt(4 * (Aₑ[i] * Aₑ[i] + μₑ * μₑ) / (δᵥ * δᵥ))
        notnanerr = δₑ==δₑ
        σ²ₑ += ifelse(notnanerr, δₑ * δₑ, ∅)
    end
    return (σ²ᵥ / max(n-corrected,0)) ± (sqrt(sum(σ²ₑ)) / max(n-corrected,0))
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