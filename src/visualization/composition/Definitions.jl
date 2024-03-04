## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots

    # Definitions
    filepath = "results/figures/composition"
    filepath_png = "results/figures/composition_pngs"

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
    legendfontsize = 18
    
    
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
    )
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