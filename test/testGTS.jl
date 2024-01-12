## --- Main function
    eon, era, period, epoch, stage = get_GTS_boundaries()

    @test isa(eon, NamedTuple)
    @test isa(era, NamedTuple)
    @test isa(period, NamedTuple)
    @test isa(epoch, NamedTuple)
    @test isa(stage, NamedTuple)


## --- Supporting functions
    boundaries = importdataset("data/boundaries_green2022.csv", ',', importas=:Tuple)

    # Test GTS to age Tuple converter
    eon = GTS_to_age(boundaries, lowercase.(boundaries.Eon))
    @test first(eon)[1].val == boundaries.Age_Ma[1]
    @test last(eon)[2].val == boundaries.Age_Ma[end]

    # Test bounds of each name
    bound = get_boundaries(boundaries.Age_Stage_above_boundary[2], boundaries, 
        boundaries.Age_Stage_above_boundary
    )
    @test bound[1] == boundaries.Age_Ma[1] ± boundaries.Age_sigma_Ma[1]
    @test bound[2] == boundaries.Age_Ma[2] ± boundaries.Age_sigma_Ma[2]

    bound = get_boundaries(boundaries.Eon[end], boundaries, boundaries.Eon)
    @test bound[1].val == boundaries.Age_Ma[end-1]
    @test bound[2].val == boundaries.Age_Ma[end]

    
## --- End of file