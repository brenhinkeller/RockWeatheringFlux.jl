using Test
include("Utilities.jl")

## --- test match_earthchem()
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
    @test cats.sed == catsmj.sed
    @test cats.ign == catsmj.ign
    @test cats.met == catsmj.met

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

## --- End of file