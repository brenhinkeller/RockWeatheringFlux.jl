## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots

    # Definitions
    filepath = "results/figures/chapter1_composition"
    filepath_png = "results/figures/chapter1_composition/pngs"

    @info """
    File IO:
    Macrostrat:  $macrostrat_io
    Geochemical: $geochem_fid
    Saving to:   $filepath/
    """

    # Plot settings 
    titlefontsize = 20
    labelfontsize = 18
    tickfontsize = 16
    
    
## --- Matched Geochemical Data and Macrostrat / Burwell
    # Indices and rock classes of matched samples
    fid = readdlm(matchedbulk_io)
    matches = Int.(vec(fid[:,1]))
    t = @. matches != 0

    match_cats = match_rocktype(string.(vec(fid[:,2]))[t]);
    include_minor!(match_cats)
    match_cats = delete_cover(match_cats)

    # Geochemical Data
    fid = h5open(geochem_fid, "r")
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
    header = Tuple(Symbol.(read(fid["type"]["macro_cats_head"])))
    data = read(fid["type"]["metamorphic_cats"])
    data = @. data > 0
    metamorphic_cats = NamedTuple{header}([data[:,i][t] for i in eachindex(header)])
    # data = read(fid["type"]["macro_cats"])
    # data = @. data > 0
    # macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    # Metamorphic samples
    include_minor!(metamorphic_cats)
    match_cats.met .|= (metamorphic_cats.sed .| metamorphic_cats.ign)

    # Major classes include minors, delete cover
    # include_minor!(macro_cats)
    # macro_cats = delete_cover(macro_cats)


## --- Unmatched Geochemical Data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    include_minor!(bulk_cats);


## --- End of file 