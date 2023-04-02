## --- test match_earthchem()
    using Test
    include("Utilities.jl")

    codes = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.3, 3.1, 3.2]
    testset = rand(codes, 100)

    bulksed, bulkign, bulkvolc, bulkplut, bulkmet, bulkmetased, bulkmetaign = match_earthchem(testset, major=false);
    bulksed_major, bulkign_major, bulkmet_major = match_earthchem(testset, major=true);

    # Rocks classified correctly
    @test (bulksed .| bulkmetased) == bulksed_major
    @test (bulkign .| bulkvolc .| bulkplut .| bulkmetaign) == bulkign_major
    @test bulkmet == bulkmet_major

    # Counts done correctly
    @test count(bulkign) == length(findall(==(3.0), testset))
    @test count(bulkvolc) == length(findall(==(3.1), testset))
    @test count(bulkplut) == length(findall(==(3.2), testset))

    @test count(bulkmet) == length(findall(==(2.0), testset))
    @test count(bulkmetased) == length(findall(==(2.1), testset))
    @test count(bulkmetaign) == length(findall(==(2.3), testset))

## --- End of file