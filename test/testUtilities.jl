## -- Test functions from Utilities.jl
    using Test
    include("../src/utilities/Utilities.jl")


## --- Get umbrella class
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    @test class_up(:ign, minorsed, minorign) == :ign 
    @test class_up(:shale, minorsed, minorign) == :sed
    @test class_up(:volc, minorsed, minorign) == :ign
    @test class_up(:carbonatite, minorsed, minorign) == :ign

    @test class_up(:basalt, minorsed, minorvolc, minorplut, minorign) == :volc
    @test class_up(:gabbro, minorsed, minorvolc, minorplut, minorign) == :plut
    
    small_type = (
        a = ("basalt", "granite", "gabbro", ),
        b = ("granite", "granodiorite", "rhyolite", ),
        c = ("schist", "gneiss", ),
    )
    @test class_up(small_type, "granite") == :a
    @test class_up(small_type, "granite", all_types=true) == (:a, :b)
    @test class_up(small_type, "gneiss") == :c
    

## --- Unmatched vs. unmetamorphosed unmatched types
    cats = get_cats(true, 4)[2]
    cats.sed .= [false, true, true, false]
    cats.met .= [true, true, false, false]

    @test find_unmatched(cats) == [false, false, false, true]
    @test find_unmetamorphosed_unmatched(cats) == .!cats.sed


## --- Macrostrat rock name matching
    # Define test set
    rocktype = [
        "volcanic rocks",                    # 1. Volc
        "major: {limestone},minor: {slate}", # 2. Carb
        "", 
        "",
    ]
    rockname = [
        "volcanic rocks",                    
        "rabbitkettle fm", 
        "hornfelsed",                        # 3. Met (keep looking!)
        "",
    ]
    rockdescrip = ["", 
        "", 
        "hornfelsed arenite and mudstone",  # 3. siliciclast, shale, met
        "silt, sand, sandstone"             # 4. siliciclast
    ]

    # Match Macrostrat names to classes
    cats = match_rocktype(rocktype, rockname, rockdescrip)

    # Test the ID'ed rock types
    @test cats.volc == [true, false, false, false]
    @test cats.carb == [false, true, false, false]
    @test cats.siliciclast == [false, false, true, true]
    @test cats.shale == [false, false, true, false]
    @test cats.met == [false, false, true, false]

    # Make sure there aren't false positives
    for k in keys(cats)
        k in (:volc, :carb, :siliciclast, :shale, :met) && continue
        @test cats[k] == [false, false, false, false]
    end

    # Get the types matched with each rock
    @test all(get_type(cats, 1, all_keys=true) == (:volc,))
    @test all(get_type(cats, 2, all_keys=true) == (:carb,))
    @test all(get_type(cats, 3, all_keys=true) == (:siliciclast, :shale, :met,))
    @test all(get_type(cats, 4, all_keys=true) == (:siliciclast,))


## --- EarthChem rock name matching
    # Define test set
    Rock_Name = [
        "basalt",       # Basalt
        "sandstone",    # Siliciclast
        ""]
    RType = [
        "volcanic",
        "siliciclastic", 
        ""
    ]
    Material = [
        "igneous", 
        "sedimentary", 
        "exotic"        # Igneous
    ]

    # Match
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    cats = match_rocktype(Rock_Name, RType, Material, (minorsed..., :sed,), 
        (minorvolc..., minorplut..., minorign..., :ign)
    )

    # Test the ID'ed rock types
    @test cats.basalt == [true, false, false]
    @test cats.siliciclast == [false, true, false]
    @test cats.ign == [false, false, true]

    # Make sure there aren't false positives
    for k in keys(cats)
        k in (:basalt, :siliciclast, :ign) && continue
        @test cats[k] == [false, false, false]
    end


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