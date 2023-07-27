## -- Test functions from Utilities.jl
    using Test

## --- match_earthchem
    codes = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
        2.0, 2.1, 2.3, 
        3.0, 3.1, 3.2
    ]
    testset = rand(codes, 1000);

    catsmj = match_earthchem(testset, major=true);
    cats = match_earthchem(testset, major=false);

    # Sedimentary counts
    @test count(cats.alluvium) == length(findall(==(1.1), testset))
    @test count(cats.siliciclast) == length(findall(==(1.2), testset))
    @test count(cats.shale) == length(findall(==(1.3), testset))
    @test count(cats.carb) == length(findall(==(1.4), testset))
    @test count(cats.chert) == length(findall(==(1.5), testset))
    @test count(cats.evaporite) == length(findall(==(1.6), testset))
    @test count(cats.phosphorite) == length(findall(==(1.7), testset))
    @test count(cats.coal) == length(findall(==(1.8), testset))
    @test count(cats.volcaniclast) == length(findall(==(1.9), testset))

    # Igneous counts
    @test count(cats.volc) == length(findall(==(3.1), testset))
    @test count(cats.plut) == length(findall(==(3.2), testset))

    # Metamorphic counts
    @test count(cats.metased) == length(findall(==(2.1), testset))
    @test count(cats.metaign) == length(findall(==(2.3), testset))

    # Major classifications
    @test catsmj.sed == cats.sed .| cats.metased
    @test catsmj.ign == cats.ign .| cats.metaign
    @test catsmj.met == cats.met .& .!cats.metased .& .!cats.metaign

    # Sub-type classifications
    @test count(cats.sed) == length(findall(==(1.0), testset)) + (count(cats.siliciclast) +
        count(cats.shale) + count(cats.carb) + count(cats.chert) + count(cats.evaporite) + 
        count(cats.phosphorite) + count(cats.coal) + count(cats.volcaniclast)
    )
    @test count(cats.ign) == length(findall(==(3.0), testset)) + (count(cats.volc) +
        count(cats.plut)
    )
    @test count(cats.met) == length(findall(==(2.0), testset)) + (count(cats.metased) + 
        count(cats.metaign)
    )


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