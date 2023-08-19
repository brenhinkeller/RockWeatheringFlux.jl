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
    # All igneous
    c, n = bincounts(mbulk.SiO2[macro_cats.ign], 40, 80, 40)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_ign, linecolor=clr_ign,
        label="All Igneous; n = $(count(macro_cats.ign))",      
        ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01) 
    )
    display(h)
    savefig("c_ign.png")

    # Volcanic
    c, n = bincounts(mbulk.SiO2[macro_cats.volc], 40, 80, 40)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_volc, linecolor=clr_volc,
        label="Volcanic; n = $(count(macro_cats.volc))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", 
        ylims=(0, round(maximum(n), digits=2)+0.01)
    )
    display(h)
    savefig("c_volc.png")

    # Plutonic
    c, n = bincounts(mbulk.SiO2[macro_cats.plut], 40, 80, 40)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_plut, linecolor=clr_plut,
        label="Plutonic; n = $(count(macro_cats.plut))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", 
        ylims=(0, round(maximum(n), digits=2)+0.01),
    )
    display(h)
    savefig("c_plut.png")

    # All sedimentary
    c, n = bincounts(mbulk.SiO2[macro_cats.sed], 0, 100, 100)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_sed, linecolor=clr_sed, 
        label="All Sedimentary; n = $(count(macro_cats.sed))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", legend=:topleft,
        ylims=(0, round(maximum(n), digits=2)+0.01),
    )
    display(h)
    savefig("c_sed.png")


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
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_rs, linecolor=clr_rs,
        label="All Igneous (resample); n = $(length(rs_ign.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", 
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_ign.png")

    # Volcanic
    c, n, = bincounts(rs_volc.SiO2, 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_rs, linecolor=clr_rs,
        label="Volcanic (resample); n = $(length(rs_volc.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_volc.png")

    # Plutonic
    c, n, = bincounts(rs_plut.SiO2, 40, 80, 80)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_rs, linecolor=clr_rs,
        label="Plutonic (resample); n = $(length(rs_plut.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_plut.png")

    # All sedimentary
    c, n, = bincounts(rs_sed.SiO2, 0, 100, 200)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, framestyle=:box, color=clr_rs, linecolor=clr_rs,
        label="Sedimentary (resample); n = $(length(rs_sed.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_sed.png")


## --- [FN CALL] modal index and Macrostrat data
    # Igneous
    ignbin = (40, 80, 40)
    iᵢ = sameindex(:ign, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)
    get_matched_samples(iᵢ, bulkidx[t], macrostrat, filter=macro_cats.ign, desc="Igneous")

    iᵥ = sameindex(:volc, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)
    iᵢ != iᵥ && get_matched_samples(iᵥ, bulkidx[t], macrostrat, filter=macro_cats.volc, desc="Volcanic")

    iₚ = sameindex(:plut, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)
    iᵢ != iₚ && get_matched_samples(iₚ, bulkidx[t], macrostrat, filter=macro_cats.plut, desc="Plutonic")

    # Sedimentary
    sedbin = (0, 100, 100)
    iₛ = sameindex(:sed, macro_cats, bulk, bulktext, bulkidx[t], bins=sedbin, hist=:off)     # All sed
    get_matched_samples(iₛ, bulkidx[t], macrostrat, filter=macro_cats.sed, desc="Sedimentary")

    # Query Macrostrat at that location:
    # https://macrostrat.org/api/mobile/map_query?lat=LAT&lng=LON&z=11

## --- End of file