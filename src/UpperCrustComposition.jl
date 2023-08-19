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

    # Indices of matched EarthChem samples from SampleMatch.jl
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0


## --- Load data for the matched EarthChem samples
    # Macrostrat
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(macrofid["rocklat"])[t],
        rocklon = read(macrofid["rocklon"])[t],
        age = read(macrofid["age"])[t],
        type = read(macrofid["type"])[t]
    )
    macro_cats = match_rocktype(macrostrat.type)
    close(macrofid)

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
        volc = [nanmean(bulk[i][macro_cats.volc]) for i in allelements],
        plut = [nanmean(bulk[i][macro_cats.plut]) for i in allelements],
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
    writedlm("$ucc_out", vcat(cols, hcat(rows, results)))

    
## --- Canadian Shield Experiment
    # Outline of shield, get points inside
    # cns = importdataset("data/shield.csv", ',', importas=:Tuple)
    # lats, lons, inshield = coords_in_shape(cns.lon, cns.lat, macrostrat.rocklon, macrostrat.rocklat)

    
## --- Daly Gap problems
    # # Plutonic
    # count(macro_cats.plut)
    # count(macro_cats.plut .& ign_cats.maf)
    # count(macro_cats.plut .& ign_cats.int)
    # count(macro_cats.plut .& ign_cats.fel)

    # # Volcanic
    # count(macro_cats.volc)
    # count(macro_cats.volc .& ign_cats.maf)
    # count(macro_cats.volc .& ign_cats.int)
    # count(macro_cats.volc .& ign_cats.fel)

    # # Uncategorized
    # count(macro_cats.ign .& .!(macro_cats.plut .| macro_cats.volc))
    # count(macro_cats.ign .& .!(macro_cats.plut .| macro_cats.volc) .& ign_cats.maf)
    # count(macro_cats.ign .& .!(macro_cats.plut .| macro_cats.volc) .& ign_cats.int)
    # count(macro_cats.ign .& .!(macro_cats.plut .| macro_cats.volc) .& ign_cats.fel)

    # # All
    # count(macro_cats.ign)
    # count(ign_cats.maf)
    # count(ign_cats.int)
    # count(ign_cats.fel)


## --- More experimentation to see what's going on with igneous rocks
    # allnames = [get_type(name_cats, i, all_keys=true) for i in eachindex(macro_cats.sed)]
    # ignnames = unique(allnames[macro_cats.ign])

    # class = Array{String}(undef, length(rocknames))
    # silica = Array{Float64}(undef, length(rocknames))
    # silica_err = Array{Float64}(undef, length(rocknames))

    # for i in eachindex(rocknames)
    #     class[i] = string(class_up(typelist, rocknames[i]))
    #     silica[i] = bulk_lookup[Symbol(rocknames[i])].SiO2.m 
    #     silica_err[i] = bulk_lookup[Symbol(rocknames[i])].SiO2.e
    # end

    # writedlm("rockname_silica.csv", vcat(["name" "class" "silica" "1 sigma"], 
    #     hcat(collect(rocknames), class, silica, silica_err)), ','
    # )


## --- End of file