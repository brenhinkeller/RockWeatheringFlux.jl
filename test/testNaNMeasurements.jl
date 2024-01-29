# ## --- Basic statistics
    a = NaN ± NaN
    b = rand() ± NaN
    c = NaN ± rand()
    d = rand() ± rand()
    e = rand() ± rand()

#     # nanadd 
#     @test nanadd(a, d) == d
#     @test nanadd(c, d) == d
#     @test nanadd(b, d) == (b+d).val ± d.err

#     # nansum
#     @test nansum([a, a, d]) == d
#     @test nansum([a, c, d, e]) == nanadd(d, e)
#     @test nansum([b, a, d, c]) == (b+d).val ± d.err

#     # nanmean 
#     @test nanmean([d, d]) == (d+d)/2  broken=true
#     @test nanmean([d, e]) == (d+e)/2
#     @test nanmean([c, d]) == d
#     @test nanmean([b, d]) == ((b.val ± 0) + d)/2

#     # Alternatively, do we want the answer to be: (b.val+d.val)/2 ± d.err
#     # That ignores the NaN entirely rather than assuming it's zero. It's harder to do...
#     # I'd also have to modify the other functions to keep that data present. BUT it's 
#     # probably more accurate... 


# ## --- Variance and standard deviation 
#     @test nanvar([a, d, e]) == nanvar([d, e])
#     @test nanvar([d, d]) == 0 ± 0

#     @test nanstd([a, d, e]) == nanstd([d, e])


## --- zeronan!
    A = [a, b, c, d]
    zeronan!(A)
    @test A[1] == 0 ± 0
    @test A[2] == 0 ± 0
    @test A[3] == 0 ± 0
    @test A[4] == d

    A = [a, b, c, d]
    @test zeronan!(A, allnans=true) == zeronan!(A)

    A = [a, b, c, d]
    zeronan!(A, allnans=false)
    @test A[1] == 0 ± 0
    @test A[2] == b.val ± 0
    @test A[3] == 0 ± c.err
    @test A[4] == d


## --- zeronan
    A = [a, b, c, d]

    A₀ = zeronan(A)
    @test A₀ != A
    @test A₀[1] == 0 ± 0
    @test A₀[2] == 0 ± 0
    @test A₀[3] == 0 ± 0
    @test A₀[4] == d

    @test zeronan(A) == zeronan(A, allnans=true)

    A₀ = zeronan(A, allnans=false)
    @test A₀[1] == 0 ± 0
    @test A₀[2] == b.val ± 0
    @test A₀[3] == 0 ± c.err
    @test A₀[4] == d

    B = rand([rand(), NaN], 4)
    while B == B
        global B = rand([rand(), NaN], 4)
    end
    B₀ = zeronan(B)

    @test B₀ != B
    @test !any(isnan, B₀)


## --- nanunzero!
    B = rand([rand(), 0, NaN], 5)
    while (B == B) | !any(==(0), B) | !any(>(0), B)
        global B = rand([rand(), 0, NaN], 5)
    end

    nanunzero!(B, 1.0)
    @test !any(==(0), B)
    @test !any(isnan, B)
    @test any(==(1), B)


## --- nanzero! 
    B = rand([rand(), 0], 4)
    while !any(==(0), B)
        global B = rand([rand(), 0], 4)
    end

    nanzero!(B)
    @test !any(==(0), B)
    @test any(isnan, B)


## --- End of file