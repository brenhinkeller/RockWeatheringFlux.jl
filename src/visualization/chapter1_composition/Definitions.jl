## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots

    # Definitions
    filepath = "results/figures/chapter1_composition"

    
## --- Matched EarthChem and Macrostrat / Burwell
    # Indices of matched samples
    fid = readdlm(matchedbulk_io)
    matches = Int.(vec(fid[:,1]))
    t = @. matches != 0

    # EarthChem
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][matches[t]] 
        for i in eachindex(header)])
    close(fid)

    # Macrostrat
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    # Major classes include minors, delete cover
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorsed
        macro_cats.sed .|= macro_cats[type]
    end
    for type in minorvolc
        macro_cats.volc .|= macro_cats[type]
    end
    for type in minorplut
        macro_cats.plut .|= macro_cats[type]
    end
    for type in minorign
        macro_cats.ign .|= macro_cats[type]
    end
    macro_cats = delete_cover(macro_cats)


## --- Unmatched Earthchem
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    for type in minorsed
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorvolc
        bulk_cats.volc .|= bulk_cats[type]
    end
    for type in minorplut
        bulk_cats.plut .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end

## --- End of file 