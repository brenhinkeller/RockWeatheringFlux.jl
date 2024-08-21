# New NaNStatistics methods for Measurements with NaN values
# Ignore any NaN ± NaN and NaN ± err values. Convert val ± NaN to val ± 0

# ## --- Summary Statistics
#     """
#     ```julia
#     nanadd(A ± B, C ± D)
#     ```

#     As `nanadd`, but if `A` or `C` is `NaN`, ignore the value. Convert any `NaN` errors to zeros.

#     # Examples
#     ```julia-repl
#     julia> nanadd(NaN ± NaN, 3 ± 4)
#     3 ± 4

#     julia> nanadd(2 ± 4, NaN ± 1)
#     2 ± 4

#     julia> nanadd(4 ± 1, 2 ± NaN)
#     6 ± 1

#     julia> nanadd(6 ± NaN, 2 ± NaN)
#     8.0 ± 0.0

#     julia> nanadd(NaN±NaN, NaN±NaN)
#     0.0 ± 0.0
#     ```

#     """
#     function NaNStatistics.nanadd(a::Measurement, b::Measurement)
#         a = a.val*(a.val==a.val) ± ifelse(a.val==a.val, a.err*(a.err==a.err), 0)
#         b = b.val*(b.val==b.val) ± ifelse(b.val==b.val, b.err*(b.err==b.err), 0)

#         return a + b
#     end
#     export nanadd


#     """
#     ```julia
#     nansum(A)
#     ```

#     Calculate the sum of a collection of measurements `A`. Ignore all `NaN` values and convert
#     all `NaN` errors to zero.

#     """
#     function NaNStatistics.nansum(A::AbstractArray{Measurement{Float64}})
#         Tₒ = Base.promote_op(+, eltype(A[1].val), Int)
#         Σ = ∅ = zero(Tₒ)
        
#         @inbounds for i ∈ eachindex(A)
#             Avalᵢ = A[i].val
#             Aerrᵢ = A[i].err*(A[i].err==A[i].err)
            
#             Σ += ifelse(Avalᵢ==Avalᵢ, (Avalᵢ ± ifelse(Aerrᵢ==Aerrᵢ, Aerrᵢ, ∅)), ∅)
#         end

#         return Σ
#     end
#     export nansum


#     """
#     ```julia
#     nanmean(A)
#     ```

#     Calculate the mean of a collection of measurements `A`. Ignore all `NaN` values and convert
#     all `NaN` errors to zero.

#     """
#     function NaNStatistics.nanmean(A::AbstractArray{Measurement{Float64}})
#         Tₒ = Base.promote_op(/, eltype(A[1].val), Int)
#         n = 0
#         Σ =  ∅ = zero(Tₒ)
        
#         @inbounds for i ∈ eachindex(A)
#             Avalᵢ = A[i].val
#             Aerrᵢ = A[i].err*(A[i].err==A[i].err)

#             Σ += ifelse(Avalᵢ==Avalᵢ, Avalᵢ ± Aerrᵢ, ∅)
#             n += (Avalᵢ==Avalᵢ)
#         end
#         return (Σ/n)
#     end
#     export nanmean


# ## --- Variance and standard deviation
#     NaNStatistics.nanstd(A::AbstractArray{Measurement{Float64}}; 
#         mean=nothing, corrected::Bool=true) = 
#         NaNStatistics.sqrt!(_nanvar(mean, corrected, A))

#     NaNStatistics.nanvar(A::AbstractArray{Measurement{Float64}}; 
#         mean=nothing, corrected::Bool=true) = _nanvar(mean, corrected, A)

#     _nanvar(::Nothing, corrected, A) = _nanvar(nanmean(A), corrected, A)

#     function _nanvar(μ::Measurement{Float64}, corrected::Bool, A::Vector{Measurement{Float64}})
#         n = 0
#         σ² = ∅ = zero(typeof(μ))
#         @inbounds for i ∈ eachindex(A)
#             δ = nanadd(A[i], - μ)
#             notnan = A[i].val==A[i].val
#             n += notnan
#             σ² += ifelse(notnan, δ * δ, ∅)
#         end
#         return σ² / max(n-corrected, 0)
#     end
#     export nanvar, nanstd


# ## --- NaNs in arrays 
    """
    ```julia
    zeronan!(A, [allnans])
    ```

    Replace all `NaN`s in A with zeros of the same type.
    
    Any element containing a `NaN` in _either_ the value or error will be replaced with 
    `0.0 ± 0.0`. To replace only `NaN`s and keep all non-`NaN` values, set 
    `allnans=false`.

    # Examples
    ```julia-repl
    julia> A = [2.0 ± NaN];

    julia> zeronan!(A)
    1-element Vector{Measurement{Float64}}:
    0.0 ± 0.0

    julia> A = [2.0 ± NaN];

    julia> zeronan!(A, allnans=false)
    1-element Vector{Measurement{Float64}}:
    2.0 ± 0.0
    ```
    """
    NaNStatistics.zeronan!(A::AbstractArray{Measurement{Float64}}; allnans::Bool=true) = 
        _zeronan!(A, static(allnans))

    function _zeronan!(A::AbstractArray{Measurement{Float64}}, allnans::True)
        ∅ = zero(eltype(A))
        @inbounds for i ∈ eachindex(A)
            A[i] = ifelse(A[i] == A[i], A[i], ∅)
        end
        return A
    end

    function _zeronan!(A::AbstractArray{Measurement{Float64}}, allnans::False)
        @inbounds for i ∈ eachindex(A)
            Aᵢ = A[i]
            A[i] = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err) 
        end
        return A
    end
    export zeronan!


    """
    ```julia
    zeronan(A; [allnans])
    ```

    Non-mutating variant of `zeronan!`.

    """
    zeronan(A::AbstractArray{Measurement{Float64}}; allnans::Bool=true) = 
        _zeronan(A, static(allnans))

    function _zeronan(A::AbstractArray{Measurement{Float64}}, allnans::True)
        out = similar(A)
        ∅ = zero(eltype(A))
        @inbounds for i ∈ eachindex(A)
            out[i] = ifelse(A[i] == A[i], A[i], ∅)
        end
        return out
    end

    function _zeronan(A::AbstractArray{Measurement{Float64}}, allnans::False)
        out = similar(A)
        @inbounds for i ∈ eachindex(A)
            Aᵢ = A[i]
            out[i] = Aᵢ.val*(Aᵢ.val==Aᵢ.val) ± Aᵢ.err*(Aᵢ.err==Aᵢ.err) 
        end
        return out
    end
    export zeronan

    function zeronan(A::AbstractArray{T}) where T
        out = similar(A)
        ∅ = zero(T)
        @inbounds for i ∈ eachindex(A)
            Aᵢ = A[i]
            out[i] = ifelse(Aᵢ==Aᵢ, Aᵢ, ∅)
        end
        return out
    end
    export zeronan


    # """
    # ```julia
    # onenan!(A)
    # ```
    # Replace all `NaN`s in A with ones of the same type.
    # """
    # function onenan!(A::StridedArray{T}) where T<:Number
    #     ∅ = one(T)
    #     @turbo for i ∈ eachindex(A)
    #         Aᵢ = A[i]
    #         A[i] = ifelse(Aᵢ==Aᵢ, Aᵢ, ∅)
    #     end
    #     return A
    # end

    # function onenan!(A::AbstractArray{T}) where T
    #     ∅ = one(T)
    #     @inbounds for i ∈ eachindex(A)
    #         Aᵢ = A[i]
    #         if isnan(Aᵢ)
    #             A[i] = ∅
    #         end
    #     end
    #     return A
    # end
    # export onenan!


    """
    ```julia
    nanunzero!(A, val)
    ```

    Replace `NaN`s and zeros in `A` with a value specified by `val`.
    """
    function nanunzero!(A::AbstractArray{T}, val) where T
        @inbounds for i ∈ eachindex(A)
            Aᵢ = A[i]
            A[i] = ifelse(Aᵢ==Aᵢ, ifelse(Aᵢ==0.0, val, Aᵢ), val)
        end
        return A
    end
    export nanunzero!


    """
    ```julia
    nanzero!(A)
    ```

    The reverse of `zeronan!`: replace all zeros with `NaN`s.
    """
    function nanzero!(A::AbstractArray{Float64})
        ∅ = 0.0
        @inbounds for i ∈ eachindex(A)
            Aᵢ = A[i]
            if Aᵢ == ∅
                A[i] = NaN
            end
        end
        return A
    end
    export nanzero!


## --- End of file