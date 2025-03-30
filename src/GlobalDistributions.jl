## --- Set up
    # Distributions of different samples, scales, etc.

    # Unique packages
    
    # Load data and base packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots, Colors

    # include("Definitions.jl");


## --- Macrostrat data 
    fid = h5open(macrostrat_io, "r")
        
    # Data
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"]),
        rocklon = read(fid["vars"]["rocklon"]),
        age = read(fid["vars"]["age"]),
        rocktype = read(fid["vars"]["rocktype"]),
        rockname = read(fid["vars"]["rockname"]),
        rockdescrip = read(fid["vars"]["rockdescrip"]),
        scale = read(fid["vars"]["scale"]),
    )

    # Rock type matches
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    close(fid)


## --- Bulk data 
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
    bulk_cats = delete_cover(bulk_cats)
    # bulk_cats = delete_volcaniclast(bulk_cats)      # FIXME: if we run this before include_minor it gets mad
    

## --- Make a cool map 
    t = Int.(rand(1:length(bulk.Sample_ID), 20_000));
    mapplot(bulk.Longitude[t], bulk.Latitude[t],
        zcolor=bulk.Age[t], color=:hawaii,
        markersize=1, msw=0, label="",
        colorbar_title = "\nAge [Ma.]",
        right_margin = (15,:px), left_margin = (15,:px)
    )


## --- Get the proportion of rocks in the dataset 
    # Straight up COPIED from SurfaceLithology.jl 😎😎😎😎😎😎

    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:end];

    # Undifferentiated igneous 
    ign_undiff = .!(bulk_cats.volc .| bulk_cats.plut .| bulk_cats.carbonatite);

    # Undifferentiated metamorphic can't be matched with anything else
    met_undiff = bulk_cats.met
    for k in keys(bulk_cats)    
        k == :met && continue
        met_undiff .&= .!bulk_cats[k]
    end

    # Sed, volc, plut, carbonatite, undiff ign, undiff met = 100
    n = length(bulk.Sample_ID)
    abundance = NamedTuple{(:sed, :volc, :plut, :carbonatite, :ign_undiff, :met_undiff)}(
        normalize!([
            count(bulk_cats.sed),                 # Sedimentary
            count(bulk_cats.volc),                # Volcanic
            count(bulk_cats.plut),                # Plutonic
            count(bulk_cats.carbonatite),         # Carbonatite
            count(bulk_cats.ign .& ign_undiff),   # Undifferentiated ign
            count(met_undiff)                     # Undifferentiated met
        ]./n.*100
    ))
    ign_total = (abundance.volc +                 # All igneous rocks 
        abundance.plut + abundance.carbonatite + 
        abundance.ign_undiff
    )
    abundance = merge(abundance, (; ign_total=ign_total))
    

    # Count the number of occurances of each minor class
    exclude_minor!(macro_cats);
    sed = float.([[count(macro_cats[k]) for k in minorsed]; count(macro_cats.sed)])
    volc = float.([[count(macro_cats[k]) for k in minorvolc]; count(macro_cats.volc)])
    plut = float.([[count(macro_cats[k]) for k in minorplut]; count(macro_cats.plut)])

    # Convert above counts to percentage of the total number of samples
    sed .= sed ./ n * 100
    volc .= volc ./ n * 100
    plut .= plut ./ n * 100
    
    # Calculate the percentage of metased and metaign rocks as the fraction of total 
    # sed / ign rocks tagged as metamorphic 
    include_minor!(macro_cats);
    metased = count(megaclass.metased) / count(macro_cats.sed)
    metaign = count(megaclass.metaign) / count(macro_cats.ign)
    metased *= abundance.sed 
    metaign *= abundance_ign

    met_total = metased + metaign + abundance.met_undiff

    # Normalize those minor class percentages to the abundance of the constituent class
    # This just means the percent abundances of e.g. sedimentary rocks will add to the 
    # total percent abundance of all sedimentary rocks that we calcualted earlier 
    sed .= sed ./ sum(sed) .* abundance.sed
    volc .= volc ./ sum(volc) .* abundance.volc
    plut .= plut ./ sum(plut) .* abundance.plut

    # Slam that bad boy together
    # We don't need to do any of the area shit because we're not separating by continent here
    ign_out = [abundance_ign, abundance.carbonatite, abundance.ign_undiff]
    sed_out = [sum(sed); sed]
    volc_out = [sum(volc); volc]
    plut_out = [sum(plut); plut]
    met_out = [met_total, metased, metaign, abundance.met_undiff]

    # Save all data to the results array in the order [total; subtypes] 
    # Be SURE to check that this is in the same order as labeled in the output 🤦
    result = [sed_out; volc_out; plut_out; ign_out; met_out]

    # Define some cool labels 
    sed_label = ["Total Sedimentary"; string.(collect(minorsed)); "Undifferentiated Sedimentary"]
    volc_label = ["Total Volcanic"; string.(collect(minorvolc)); "Undifferentiated Volcanic"]
    plut_label = ["Total Plutonic"; string.(collect(minorplut)); "Undifferentiated Plutonic"]
    ign_label = ["Total Igneous", "Carbonatite", "Undifferentiated Igneous"]
    met_label = ["Total Metamorphic", "Metasedimentary", "Metaigneous", "Undifferentiated Metamorphic"]   
    cols = ["Lithology" "Abundance"]

    labels = [sed_label; volc_label; plut_label; ign_label; met_label];

    # Wheeeeeeeeeeeeeeeee! We could export this if we wanted 
    bulkabundance = vcat(cols, hcat(labels, result))

    # OK, load the global file and see how far off we are
    surfabundance = readdlm(mapped_surface_lith, ',')[:,end]
    difference = bulkabundance[2:end,2] .- surfabundance[2:end]

    smashed = (hcat(bulkabundance[2:end,2], surfabundance[2:end], difference))
    smashed = round.(smashed, sigdigits=3)
    display(hcat(bulkabundance[:,1], vcat(["Geochem" "Maps" "Difference"], smashed)))



## --- Make some maps 
    # # Everything we got! (Only one point had no response)
    # h = mapplot(macrostrat.rocklon, macrostrat.rocklat, 
    #     label="",
    #     markersize=1, 
    #     msc=:auto,
    #     color=:yellow,
    #     title="All Data",
    # )
    # display(h)
    # savefig(h, "dev/macrostratscales_alldata.png")

    # # By Scale 
    # scales = ("tiny", "small", "medium", "large")
    # colorscales = (:firebrick, :darkorange, :forestgreen, :royalblue)
    # for i in eachindex(scales)
    #     t = macrostrat.scale .== scales[i]
    #     h = mapplot(macrostrat.rocklon[t], macrostrat.rocklat[t], 
    #         label="",
    #         markersize=1, 
    #         msc=:auto,
    #         color=colorscales[i],
    #         title="$(scales[i])",
    #     )
    #     display(h)
    #     savefig(h, "dev/macrostratscales_$(scales[i]).png")
    # end


## --- Small scale as a function of age 
    


## --- End of File 