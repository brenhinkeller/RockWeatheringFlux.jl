# Visualize distributions of data; histograms.
#
# CTRL+K CTR+[ to collapse all sections

## --- Set up
    # Base computational packages
    using StatGeochem
    using HDF5
    using LoopVectorization
    using Measurements
    using Static
    using DelimitedFiles
    using ProgressMeter

    using Plots

    # Local utilities
    include("utilities/Utilities.jl")

    # More definitions
    bins = (ign = (40, 80, 40), sed = (0, 100, 100), met = (25, 100, 75))
    minorsed, minorign, minormet = get_minor_types()


## --- [DATA] Matched EarthChem and Macrostrat samples
    # Get indicies of matched samples
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Matched Macrostrat samples
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocktype = read(macrofid["rocktype"])[t],
        rockname = read(macrofid["rockname"])[t],
        rockdescrip = read(macrofid["rockdescrip"])[t],
        rocklat = read(macrofid["rocklat"])[t],
        rocklon = read(macrofid["rocklon"])[t],
        age = read(macrofid["age"])[t],
        type = read(macrofid["type"])[t]
    )
    close(macrofid)
    macro_cats = match_rocktype(macrostrat.type)

    # Earthchem data
    bulkfid = h5open("output/bulk.h5", "r")
        # Data
        header = read(bulkfid["bulk"]["header"])
        data = read(bulkfid["bulk"]["data"])
        bulktype = read(bulkfid["bulk"]["type"])

        # Metadata
        path = bulkfid["bulktext"]["sampledata"]
        headertext = read(path["header"])
        index = read(path["index"])
        target = ["Rock_Name", "Type", "Material"]
        targetind = [findall(==(i), headertext)[1] for i in target]
        bulktext = NamedTuple{Tuple(Symbol.(target))}(
            [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] 
                for i in eachindex(target)]
        )
    close(bulkfid)

    # Matched samples and all data
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    bulk_cats = match_rocktype(bulktype)
    mbulk_cats = match_rocktype(bulktype[bulkidx[t]])


## --- SiO₂ distribution by rock type
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(macro_cats) - 1)

    # Create all plots
    rocks = collect(keys(macro_cats))
    for i in eachindex(rocks)
        r = rocks[i]
        if r == :cover
            continue
        elseif r == :sed || r == :ign || r == :met
            type = r
        else
            r in minorsed && (type = :sed)
            r in minorign && (type = :ign)
            r in minormet && (type = :met)
        end

        c, n = bincounts(mbulk.SiO2[macro_cats[r]], bins[type]...)
        n = float(n) ./ nansum(float(n) .* step(c))
        h = plot(c, n, seriestype=:bar, framestyle=:box, color=colors[r], linecolor=colors[r],
            label="$(string(r)); n = $(count(macro_cats[r]))",      
            # ylabel="Weight", xlabel="SiO2 [wt.%]",
            ylims=(0, round(maximum(n), digits=2)+0.01) 
        )

        fig[i] = h
    end

    # Put into a layout
    h = plot(fig..., layout=(5,3), size=(2000, 2000))
    display(h)
    savefig(h, "distributions.png")


## --- [FN CALL] modal index and Macrostrat data
    # Igneous
    ignbin = (40, 80, 40)
    iᵢ = sameindex(:ign, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)
    # get_matched_samples(iᵢ, bulkidx[t], macrostrat, filter=macro_cats.ign, desc="Igneous")

    iᵥ = sameindex(:volc, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)
    # iᵢ != iᵥ && get_matched_samples(iᵥ, bulkidx[t], macrostrat, filter=macro_cats.volc, desc="Volcanic")

    iₚ = sameindex(:plut, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)
    # iᵢ != iₚ && get_matched_samples(iₚ, bulkidx[t], macrostrat, filter=macro_cats.plut, desc="Plutonic")

    # Sedimentary
    sedbin = (0, 100, 100)
    iₛ = sameindex(:sed, macro_cats, bulk, bulktext, bulkidx[t], bins=sedbin, hist=:off)     # All sed
    # get_matched_samples(iₛ, bulkidx[t], macrostrat, filter=macro_cats.sed, desc="Sedimentary")

    metbin = (0, 100, 100)
    iₘ = sameindex(:met, macro_cats, bulk, bulktext, bulkidx[t], bins=metbin, hist=:off)     # All met

    # Query Macrostrat at that location:
    # https://macrostrat.org/api/mobile/map_query?lat=LAT&lng=LON&z=11

    # exit()


## --- Igneous SiO₂ by rock age
    t = @. (mbulk.Age > 2500) & macro_cats.ign
    c, n = bincounts(mbulk.SiO2[t], bins.ign...)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = plot(c, n, seriestype=:bar, framestyle=:box, color=colors.ign, linecolor=colors.ign,
        label="Archean; n = $(count(t))",      
        # ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01) 
    )

    t = @. (2500 > mbulk.Age > 541) & macro_cats.ign
    c, n = bincounts(mbulk.SiO2[t], bins.ign...)
    n = float(n) ./ nansum(float(n) .* step(c))
    h2 = plot(c, n, seriestype=:bar, framestyle=:box, color=colors.ign, linecolor=colors.ign,
        label="Proterozoic; n = $(count(t))",      
        # ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01) 
    )

    t = @. (541 > mbulk.Age) & macro_cats.ign
    c, n = bincounts(mbulk.SiO2[t], bins.ign...)
    n = float(n) ./ nansum(float(n) .* step(c))
    h3 = plot(c, n, seriestype=:bar, framestyle=:box, color=colors.ign, linecolor=colors.ign,
        label="Phanerozoic; n = $(count(t))",      
        # ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01) 
    )

    # Compile
    h = plot(h1, h2, h3, layout=(1,3), size=(2000, 450), bottom_margin=(50, :px), 
        xlabel="SiO₂ [wt.%]"
    )
    display(h)
    savefig(h, "ign_ages.png")


## --- [DATA] Resampled EarthChem SiO₂ (spatial density only)
    # Igneous
    rs_ign = importdataset("output/resampled/rs_ign.tsv", '\t', importas=:Tuple)
    rs_volc = importdataset("output/resampled/rs_volc.tsv", '\t', importas=:Tuple)
    rs_plut = importdataset("output/resampled/rs_plut.tsv", '\t', importas=:Tuple)

    # Sedimentary
    rs_sed = importdataset("output/resampled/rs_sed.tsv", '\t', importas=:Tuple)


## --- Resampled SiO₂ distribution by rock type
    # All igneous
    c, n, = bincounts(rs_ign.SiO2, 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=c_rs, linecolor=c_rs,
        label="All Igneous (resample); n = $(length(rs_ign.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", 
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_ign.png")

    # Volcanic
    c, n, = bincounts(rs_volc.SiO2, 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=c_rs, linecolor=c_rs,
        label="Volcanic (resample); n = $(length(rs_volc.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_volc.png")

    # Plutonic
    c, n, = bincounts(rs_plut.SiO2, 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=c_rs, linecolor=c_rs,
        label="Plutonic (resample); n = $(length(rs_plut.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_plut.png")

    # All sedimentary
    c, n, = bincounts(rs_sed.SiO2, 0, 100, 200)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=c_rs, linecolor=c_rs,
        label="Sedimentary (resample); n = $(length(rs_sed.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_sed.png")


## --- End of file