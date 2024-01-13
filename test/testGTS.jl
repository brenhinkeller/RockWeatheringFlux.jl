## --- Main function
    # Get a list of bounds
    bounds_GTS = get_GTS_boundaries();
    @test isa(bounds_GTS.eon, NamedTuple)
    @test isa(bounds_GTS.era, NamedTuple)
    @test isa(bounds_GTS.period, NamedTuple)
    @test isa(bounds_GTS.epoch, NamedTuple)
    @test isa(bounds_GTS.stage, NamedTuple)

    # Assign bounds to a name
    boundaries = importdataset("data/boundaries_green2022.csv", ',', importas=:Tuple)

    bound = RockWeatheringFlux.get_boundaries("holocene", boundaries, lowercase.(boundaries.Epoch))
    @test assign_GTS_age("holocene", names_GTS, bounds_GTS) == Tuple(bound)
    
    bound = RockWeatheringFlux.get_boundaries("early ordovician", boundaries, lowercase.(boundaries.Epoch))
    @test assign_GTS_age("early ordovician", names_GTS, bounds_GTS) == Tuple(bound)


## --- Supporting functions
    # GTS to age Tuple converter
    eon = RockWeatheringFlux.GTS_to_age(boundaries, lowercase.(boundaries.Eon))
    @test first(eon)[1].val == boundaries.Age_Ma[1]
    @test last(eon)[2].val == boundaries.Age_Ma[end]

    stage = RockWeatheringFlux.GTS_to_age(boundaries, lowercase.(boundaries.Age_Stage_above_boundary))
    @test !("" in keys(stage))

    # Query upper and lower bounds of each name
    bound = RockWeatheringFlux.get_boundaries(boundaries.Age_Stage_above_boundary[2], boundaries, 
        boundaries.Age_Stage_above_boundary
    )
    @test bound[1] == boundaries.Age_Ma[1] ± boundaries.Age_sigma_Ma[1]
    @test bound[2] == boundaries.Age_Ma[2] ± boundaries.Age_sigma_Ma[2]

    bound = RockWeatheringFlux.get_boundaries(boundaries.Eon[end], boundaries, boundaries.Eon)
    @test bound[1].val == boundaries.Age_Ma[end-1]
    @test bound[2].val == boundaries.Age_Ma[end]

    # Get time division (era, eon, etc.) for a name
    bounds_GTS = get_GTS_boundaries();
    names_GTS = (
        eon = string.(keys(bounds_GTS.eon)),
        era = string.(keys(bounds_GTS.era)),
        period = string.(keys(bounds_GTS.period)),
        epoch = string.(keys(bounds_GTS.epoch)),
        stage = string.(keys(bounds_GTS.stage)),
    
    )
    @test :eon == RockWeatheringFlux.get_GTS_division(names_GTS.eon[2], names_GTS)
    @test :era == RockWeatheringFlux.get_GTS_division(names_GTS.era[4], names_GTS)
    @test :period == RockWeatheringFlux.get_GTS_division(names_GTS.period[10], names_GTS)
    @test :epoch == RockWeatheringFlux.get_GTS_division(names_GTS.epoch[4], names_GTS)
    @test :stage == RockWeatheringFlux.get_GTS_division(names_GTS.stage[3], names_GTS)
    

## --- End of file