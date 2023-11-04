## -- Test functions from Utilities.jl
    using Test
    include("../src/utilities/Utilities.jl")

## --- Get umbrella class
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    @test class_up(:ign, minorsed, minorign) == :ign 
    @test class_up(:shale, minorsed, minorign) == :sed
    @test class_up(:volc, minorsed, minorign) == :ign
    @test class_up(:basalt, minorsed, minorign) == :ign
    @test class_up(:carbonatite, minorsed, minorign) == :ign

    
## --- Rock name matching
    rocktype = ["sedimentary and volcanic rocks", "major: {limestone},minor: {slate}", "", "",]
    rockname = ["precambrian-phanerozoic sedimentary and volcanic rocks", "rabbitkettle fm", 
        "hodgkinson formation - hornfelsed", "",]
    rockdescrip = ["", "", "hornfelsed arenite and mudstone", "silt, sand, sandstone",]

    # Match Macrostrat names to classes
    cats = match_rocktype(rocktype, rockname, rockdescrip; source=:macrostrat, unmultimatch=false)

    @test cats.siliciclast == [false, false, false, true]
    @test cats.shale == [false, false, false, false, ]
    @test cats.carb == [false, true, false, false, ]
    @test cats.chert == [false, false, false, false, ]
    @test cats.evaporite == [false, false, false, false, ]
    @test cats.coal == [false, false, false, false, ]
    @test cats.phosphorite == [false, false, false, false, ]
    @test cats.volcaniclast == [false, false, false, false, ]
    @test cats.sed == [true, true, false, true,]
    @test cats.volc == [true, false, false, false, ]
    @test cats.plut == [false, false, false, false, ]
    @test cats.ign == [true, false, false, false, ]
    @test cats.metased == [false, false, true, false, ]
    @test cats.metaign == [false, false, false, false, ]
    @test cats.met == [false, false, true, false, ]
    @test cats.cover == [false, false, false, false, ]

    @test all(get_type(cats, 1, all_keys=true) == (:sed, :volc, :ign))
    @test all(get_type(cats, 2, all_keys=true) == (:carb, :sed))
    @test all(get_type(cats, 3, all_keys=true) == (:metased, :met))
    @test all(get_type(cats, 4, all_keys=true) == (:siliciclast, :sed))

    # Match Macrostrat rock names to names
    cats = match_rockname(rocktype, rockname, rockdescrip)

    @test all(get_type(cats, 1, all_keys=true) == (:sediment, :volcanic, :volcan))
    @test all(get_type(cats, 2, all_keys=true) == (:limestone,))
    @test all(get_type(cats, 3, all_keys=true) == (:hornfels,))
    @test all(get_type(cats, 4, all_keys=true) == (:sand, :silt))

    # Match EarthChem rock names to classes
    Rock_Name = ["basalt", "sandstone", ""]
    Type = ["volcanic", "siliciclastic", ""]
    Material = ["igneous", "sedimentary", "exotic"]

    cats = match_rocktype(Rock_Name, Type, Material; source=:earthchem, unmultimatch=false)

    @test cats.siliciclast == [false, true, false]
    @test cats.shale == [false, false, false]
    @test cats.carb == [false, false, false]
    @test cats.chert == [false, false, false]
    @test cats.evaporite == [false, false, false]
    @test cats.coal == [false, false, false]
    @test cats.phosphorite == [false, false, false]
    @test cats.volcaniclast == [false, false, false]
    @test cats.sed == [false, true, false]
    @test cats.volc == [true, false, false]
    @test cats.plut == [false, false, true]
    @test cats.ign == [true, false, true]
    @test cats.metased == [false, false, false]
    @test cats.metaign == [false, false, false]
    @test cats.met == [false, false, false]
    @test cats.cover == [false, false, false]


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