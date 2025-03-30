## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots, Colors

    # Definitions
    filepath = "results/figures/composition"
    filepath_png = "results/figures/composition_pngs"

    @info """
    File IO:
    Macrostrat:  $macrostrat_io
    Geochemical: $geochem_fid
    Saving to:   $filepath/
    """

    # Standardized color schemes 
    colors_covariance = (
        a = parse(Colorant, "#5b80bb"),             # Based on managua blues
        b = parse(Colorant, "#bb5b80"),
    )
    colors_contrast = (
        a = parse(Colorant, "#913b5b"),             # Based on managua oranges
        b = parse(Colorant, "#3b8691"),
    )
    p = Plots.palette(:managua, 6)
    colors_source = (
        rudnick = parse(Colorant, "#f0a400"),       # Rudnick and Gao, 2014
        shaw = p[2],                                # Shaw et al., 1967, 1976
        condie = p[3],                              # Condie, 1993
        gao = p[4],                                 # Gao et al., 1998
        pease = p[5],                               # Pease et al., 2023
        glorise = :cadetblue,                       # Muller et al., 2021
    )
    colors_dark = (                                 # Modified USGS colors
        sed = parse(Colorant, "#206341"),
        shale = parse(Colorant, "#44c283"),
        ign = colors.ign,
        volc = parse(Colorant, "#ff6bbb"),
        basalt = parse(Colorant, "#ca895d"),
        plut = parse(Colorant, "#fb3448"),
        granite = parse(Colorant, "#f36e7a"),
    )


## --- Matched Geochemical Data and Macrostrat / Burwell
    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)

    # Macrostrat
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        scale = read(fid["vars"]["scale"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)


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