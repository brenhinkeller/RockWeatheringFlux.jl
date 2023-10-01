## --- Set up
    # Packages
    using StatGeochem
    using DelimitedFiles
    using Measurements
    using HDF5
    using ProgressMeter
    using LoopVectorization
    using Static

    # Local utilities
    include("utilities/Utilities.jl")

    # Get igneous rock silica definitions
    ignsilica = get_ignsilica()


## --- Load data for the matched EarthChem samples
    # Indices of matched EarthChem samples from SampleMatch.jl
    data = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(data[:,1]))
    t = @. bulkidx != 0

    # Matched types, majors inclusive of minors
    bulktype = string.(vec(data[:,2]))
    macro_cats = match_rocktype(bulktype[t])

    # Macrostrat
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
    )
    close(fid)

    # Earthchem
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(fid)

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
        volc = [nanmean(bulk[i][macro_cats.volc]) for i in allelements],
        plut = [nanmean(bulk[i][macro_cats.plut]) for i in allelements],
        fel = [nanmean(bulk[i][ign_cats.fel]) for i in allelements],
        int = [nanmean(bulk[i][ign_cats.int]) for i in allelements],
        maf = [nanmean(bulk[i][ign_cats.maf]) for i in allelements],
        ultramaf = [nanmean(bulk[i][ign_cats.ultramaf]) for i in allelements],
        ultrafel = [nanmean(bulk[i][ign_cats.ultrafel]) for i in allelements],
    )

    # Save to file
    rows = string.(allelements)
    cols = hcat("", string.(reshape(collect(keys(UCC)), (1, length(keys(UCC))))))
    results = Array{Float64}(undef, (length(allelements), length(UCC)))
    for i in eachindex(keys(UCC))
        results[:,i] = UCC[i]
    end
    writedlm("$ucc_out", vcat(cols, hcat(rows, results)))

    # Terminal printout
    majorcomp = round.([UCC.bulk[i] for i in eachindex(majors)], sigdigits=2)

    @info """ Bulk crustal composition:
    $(join(majors, " \t "))
    $(join(majorcomp, " \t "))

    Total: $(round(nansum(UCC.bulk[1:length(majors)]), sigdigits=3))%
    Total (less CO₂ and H₂O): $(round(nansum(UCC.bulk[1:length(majors)-2]), sigdigits=3))%
    """


## --- End of file