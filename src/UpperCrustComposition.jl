## --- Set up
    # Packages
    using StatGeochem
    using DelimitedFiles
    using Measurements
    using HDF5

    # Local utilities
    include("Utilities.jl")
    include("NaNMeasurements.jl")

    # Define igneous rock compositions by silica (from Keller and Schoene, 2012)
    ignsilica = (
        fel = (62, 74),      # Felsic (low exclusive, high inclusive)
        int = (51, 62),      # Intermediate
        maf = (43, 51),      # Mafic
    )

    # Indices of matched EarthChem samples from SampleMatch.jl
    bulkidx = Int.(vec(readdlm("output/bulkidx.tsv")))
    t = @. bulkidx != 0     # Exclude samples with missing data


## --- Load data for the matched EarthChem samples
    # Macrostrat
    macrostrat = importdataset("output/pregenerated_responses.tsv", '\t', importas=:Tuple)
    macro_cats = match_rocktype(macrostrat.rocktype[t], macrostrat.rockname[t], 
        macrostrat.rockdescrip[t]
    )
    known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    total_known = count(known_rocks)

    # Earthchem
    bulkfid = h5open("output/bulk.h5", "r")
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(bulkfid)

    # Elements of interest
    majors, minors = get_elements()
    allelements = [majors; minors]

    # Define BitVectors for igneous rocks by silica content
    ign_cats = (
        fel = (@. macro_cats.ign & (ignsilica.fel[1] < bulk.SiO2 <= ignsilica.fel[2])),
        int = (@. macro_cats.ign & (ignsilica.int[1] < bulk.SiO2 <= ignsilica.int[2])),
        maf = (@. macro_cats.ign & (ignsilica.maf[1] < bulk.SiO2 <= ignsilica.maf[2])),
        ultramaf = (@. macro_cats.ign & (ignsilica.maf[1] > bulk.SiO2)),
        ultrafel = (@. macro_cats.ign & (ignsilica.fel[2] < bulk.SiO2)),
    )


## --- Compute and export composition of exposed crust!
    UCC = (
        bulk = [nanmean(bulk[i]) for i in allelements],
        sed = [nanmean(bulk[i][macro_cats.sed]) for i in allelements],
        met = [nanmean(bulk[i][macro_cats.met]) for i in allelements],
        ign = [nanmean(bulk[i][macro_cats.ign]) for i in allelements],
        fel = [nanmean(bulk[i][ign_cats.fel]) for i in allelements],
        int = [nanmean(bulk[i][ign_cats.int]) for i in allelements],
        maf = [nanmean(bulk[i][ign_cats.maf]) for i in allelements],
        ultramaf = [nanmean(bulk[i][ign_cats.ultramaf]) for i in allelements],
        ultrafel = [nanmean(bulk[i][ign_cats.ultrafel]) for i in allelements],
    )

    results = Array{Float64}(undef, (length(allelements), length(UCC)))
    for i in eachindex(keys(UCC))
        results[:,i] = UCC[i]
    end

    rows = string.(allelements)
    cols = hcat("", string.(reshape(collect(keys(UCC)), (1, length(keys(UCC))))))
    writedlm("results/exposedcrust.tsv", vcat(cols, hcat(rows, results)))

    
## --- End of file