## -- NaNMeasurements.jl
    using StatGeochem
    using Measurements
    using Statistics
    using Test

    using Static
    using LoopVectorization
    include("../src/utilities/NaNMeasurements.jl")

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


## --- Summary statistics: simple cases
    @test nanadd(A[1], A[3]).val == (B[1] + B[3]).val
    @test nansum(A).val == sum(B).val
    @test nanmean(A).val == mean(B).val
    @test nanstd(A).val ≈ std(B).val
    @test nanvar(A).val ≈ var(B).val

    @test nanadd(A[1], A[3]).err == (B[1] + B[3]).err
    @test nansum(A).err == sum(B).err
    @test nanmean(A).err == mean(B).err
    @test nanstd(A).err ≈ std(B).err
    @test nanvar(A).err ≈ var(B).err


## --- NaN replacement
    testA = copy(A)
    @test zeronan!(testA) == B_znan_all

    testA = copy(A)
    @test zeronan!(testA, true) == B_znan_all

    testA = copy(A)
    @test zeronan!(testA, false) == B_znan_par


## --- Zero replacement
    A = [0, 0, 2, 0.0]
    nanzero!(A)
    @test isnan(A[1])
    @test isnan(A[2])
    @test !isnan(A[3])
    @test isnan(A[4])

## --- End of file