## -- NaNMeasurements.jl
    using StatGeochem
    using Measurements
    using Statistics
    using Test

## --- Create test arrays
    A = collect(1:10) .± reverse!(collect(1:10))
    B = copy(A)
    B_znan_all = copy(B)
    B_znan_par = copy(B)

    A[1] = A[1].val ± NaN
    B[1] = B[1].val ± 0.0
    B_znan_all[1] = 0.0 ± 0.0
    B_znan_par[1] = B[1].val ± 0.0

    A[3] = NaN ± A[3].err
    B[3] = 0.0 ± B[3].err
    B_znan_all[3] = 0.0 ± 0.0
    B_znan_par[3] = 0.0 ± B[3].err

    A[10] = NaN ± NaN
    deleteat!(B, 10)
    B_znan_all[10] = 0.0 ± 0.0
    B_znan_par[10] = 0.0 ± 0.0


## --- Summary statistics
    @test nanadd(A[1], A[3]) == B[1] + B[3]
    @test nansum(A) == sum(B)
    @test nanmean(A) == mean(B)
    @test nanstd(A) ≈ std(B)
    @test nanvar(A) ≈ var(B)


## --- NaN replacement
    testA = copy(A)
    @test zeronan!(testA) == B_znan_all

    testA = copy(A)
    @test zeronan!(testA, true) == B_znan_all

    testA = copy(A)
    @test zeronan!(testA, false) == B_znan_par

## --- End of file