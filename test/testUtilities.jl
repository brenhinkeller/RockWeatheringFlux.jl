## -- Test functions from Utilities.jl
    using Test


## --- Match rock names / types / descriptions to defined rock classes
    rocktype = ["sedimentary rocks", "sedimentary rocks", "sedimentary rocks",
        "flood basalt(s); mafic volcanic rocks", "sedimentary rocks", 
        "crystalline metamorphic rocks", "", "", "intrusive igneous rocks", 
        "crystalline metamorphic rocks", "na"
    ]
    rockname = ["cenozoic sedimentary rocks", "precambrian sedimentary rocks", 
        "paleozoic sedimentary rocks", "mesozoic volcanic rocks", 
        "neoproterozoic sedimentary rocks", "archean crystalline metamorphic rocks", "", 
        "hodgkinson formation - hornfelsed", "archean intrusive rocks", 
        "paleoproterozoic crystalline metamorphic rocks", "na"
    ]
    rockdescrip = ["", "", "", "", "", "", "", "hornfelsed arenite and mudstone", "", "","na"]

    cats = match_rocktype([rocktype[8]], [rockname[8]], [rockdescrip[8]]; source=:macrostrat)
    @test cats.met == [true]
    @test cats.metased == [true]

## --- Points in polygon
    # Define simple shape
    polyx = [0, 10, 10, 0, 0]
    polyy = [0, 0, 10, 10, 0]
    x = [5, 0, 10, 10, 15, -1, 3]
    y = [5, 0, 5, 10, 10, -1, 20]

    xin, yin, = points_in_shape(polyx, polyy, x, y)
    @test xin == [5, 0, 10, 10]
    @test yin == [5, 0, 5, 10]

    # Test for a coordinate system
    polylats = [  49.6,   51.2,   46.7,   45.6]
    polylons = [-118.5, -104.3, -104.3, -112.8]
    lat = [50.3,     80.0,   49.3, 0]
    lon = [-108.6, -117.2, -118.1, 0]

    lon_in, lat_in, = coords_in_shape(polylons, polylats, lon, lat)
    @test lat_in == [50.3, 49.3]
    @test lon_in == [-108.6, -118.1]

    # Crossing the antimeridian
    polylats = [ 45,   45,  -45, -45]
    polylons = [160, -160, -160, 160]
    lat = [ 20,  85,   15, -38,  -30,    0]
    lon = [170, 180, -180, 180, -170, -140]

    lon_in, lat_in, = coords_in_shape(polylons, polylats, lon, lat)
    @test lat_in == [ 20,   15, -38,  -30]
    @test lon_in == [170, -180, 180, -170]

    # Over a pole
    polylats = [80,  80,  80, 80,80,  80,  80, 80]
    polylons = [45, 90, 135, 180, -135, -90, -45, 0]
    lat = [90, 90,   70, 85, 56, -88]
    lon = [0, -80, -180, 36, 97, 156]

    lon_in, lat_in, = coords_in_shape(polylons, polylats, lon, lat)
    @test lat_in == [90, 90, 85, -88]
    @test lon_in == [0, -80, 36, 156]


## --- End of file