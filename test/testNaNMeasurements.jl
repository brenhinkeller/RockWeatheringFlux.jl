## -- NaNMeasurements.jl
    # Packages
    using StatGeochem
    using Measurements
    using Test

    # Create test arrays
    A = rand(10) .± rand(10)    # Some values will be NaN
    B = copy(A)                 # NaN values of A removed or replaced with 0
    S = collect(1:10)           # Possible indices

    i_nanerr = rand(S)
    deleteat!(S, findall(x->x==i_nanerr, S))
    A[i_nanerr] = A[i_nanerr].val ± NaN
    B[i_nanerr] = B[i_nanerr].val ± 0

    i_nanval = rand(S)
    deleteat!(S, findall(x->x==i_nanval, S))
    A[i_nanval] = NaN ± A[i_nanval].err
    B[i_nanval] = 0 ± B[i_nanval].err

    i_nanall = rand(S)
    deleteat!(S, findall(x->x==i_nanall, S))
    A[i_nanall] = NaN ± NaN
    deleteat!(B, i_nanall)


## --- nanadd


## --- nansum


## --- nanmean


## --- nanvar


## --- nanstd


## --- End of file