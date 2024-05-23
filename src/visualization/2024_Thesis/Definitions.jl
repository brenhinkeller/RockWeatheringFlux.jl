## --- Set up
    # Packages 
    using RockWeatheringFlux 
    using HDF5, DelimitedFiles
    using Plots

    # Filepath
    filepath = "results/figures/Wetterhahn"

    # Color palette
    cscheme = :managua
    pal = Plots.palette(cscheme, 10)


## --- Matched geochemical data 
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
    

## --- Lithologic data 
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    close(fid)
    

## --- End of file 