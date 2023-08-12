## --- Set up
    # CTRL+K CTR+[ to collapse all sections
    # Note that loading Makie packages and Plots packages may cause the code to run
    # incorrectly

    # Base computational packages
    using StatGeochem
    using HDF5
    using LoopVectorization
    using Measurements
    using Static
    using DelimitedFiles
    using ProgressMeter

    # Local utilities
    include("utilities/Utilities.jl")

    colorsch = :jet1

## --- Location and age of EarthChem (bulk) points
    # Packages
    using MAT
    using CairoMakie
    using GeoMakie
    using ImageMagick

    # From wt.% restricted to 84-104 total wt.% 
    bulkfid = h5open("output/bulk.h5", "r")
        header = read(bulkfid["bulk"]["header"])
        data = read(bulkfid["bulk"]["data"])
    close(bulkfid)
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, bulk.Longitude, bulk.Latitude, color=bulk.Age,
        colormap=colorsch, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("bulk_restricted.png", f)

    # All data
    bulk = matread("data/bulk.mat")["bulk"]
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))

    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, vec(bulk.Longitude), vec(bulk.Latitude), color=vec(bulk.Age),
        colormap=colorsch, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("bulk_all.png", f)

## --- Location and age of bulk points matched to Macrostrat samples
    # Packages
    using CairoMakie
    using GeoMakie
    using ImageMagick

    # Get indices of matched samples
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Get matched Earthchem (bulk) samples
    bulkfid = h5open("output/bulk.h5", "r")
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(bulkfid)

    # Plot
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, vec(bulk.Longitude), vec(bulk.Latitude), color=vec(bulk.Age),
        colormap=colorsch, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("bulk_matched.png", f)


## --- Location and age of Macrostrat samples
    # Packages
    using CairoMakie
    using GeoMakie
    using ImageMagick

    # Get indices of matched samples
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Get Macrostrat samples with a known match
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(macrofid["rocklat"])[t],
        rocklon = read(macrofid["rocklon"])[t],
        age = read(macrofid["age"])[t]
    )
    close(macrofid)

    # Plot
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, 
        color=macrostrat.age, colormap=colorsch, markersize = 3
    )
    Colorbar(f[1,2], h, label = "Age [Ma]", height = Relative(0.9))
    display(f)
    save("macrostrat.png", f)


## --- [DATA MISSING] Difference in SiO₂ and location between bulk and Macrostrat
    # # Location data
    # data = importdataset("dist.csv", ',', importas=:Tuple)
    # macrostrat = (rocklat = data.lat_m, rocklon = data.lon_m)
    # dist = data.dist

    # # silica data
    # silica = importdataset("silica.csv", ',', importas=:Tuple)

    # # Distance between sample and matched point
    # f = Figure(resolution = (1200, 600))
    # ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    # h = CairoMakie.scatter!(ax, macrostrat.rocklon, macrostrat.rocklat, color=dist,
    #     markersize = 3, alpha=0.5
    # )
    # Colorbar(f[1,2], h, label = "Arc Degrees from Matched Point", height = Relative(0.9))
    # display(f)

    # # Silica content of matched samples
    # t = silica.SiO2 .> 50

    # f = Figure(resolution = (1200, 600))
    # ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    # h = CairoMakie.scatter!(ax, silica.lon[t], silica.lat[t], color=silica.SiO2[t],
    #     markersize = 3, alpha=0.5
    # )
    # Colorbar(f[1,2], h, label = "SiO₂ [wt/%]", height = Relative(0.9))
    # display(f)


## --- [DATA] Matched EarthChem and Macrostrat samples
    # Get indicies of matched samples
    bulkidx = Int.(vec(readdlm("$matchedbulk_io")))
    t = @. bulkidx != 0

    # Matched Macrostrat samples
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        # rocktype = read(macrofid["rocktype"])[t],
        # rockname = read(macrofid["rockname"])[t],
        # rockdescrip = read(macrofid["rockdescrip"])[t]
        age = read(macrofid["age"])[t],
        type = read(macrofid["typecategory"])[t]
    )
    close(macrofid)
    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
    #     macrostrat.rockdescrip
    # )
    macro_cats = match_rocktype(macrostrat.type)

    # Earthchem for each matched sample
    bulkfid = h5open("output/bulk.h5", "r")
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    close(bulkfid)

    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])


## --- SiO₂ distribution by igneous rock type
    using Plots

    # All igneous
    c, n = bincounts(bulk.SiO2[macro_cats.ign], 40, 80, 40)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="All Igneous; n = $(count(macro_cats.ign))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, 
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_ign.png")

    # Volcanic
    c, n = bincounts(bulk.SiO2[macro_cats.volc], 40, 80, 40)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Volcanic; n = $(count(macro_cats.volc))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box,
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_volc.png")

    # Plutonic
    c, n = bincounts(bulk.SiO2[macro_cats.plut], 40, 80, 40)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Plutonic; n = $(count(macro_cats.plut))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box,
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_plut.png")


## --- SiO₂ distribution by sedimentary rock type 
    # All sedimentary
    c, n = bincounts(bulk.SiO2[macro_cats.sed], 0, 100, 100)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="All Sedimentary; n = $(count(macro_cats.sed))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, 
        ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_sed.png")


## --- [DATA] Resampled EarthChem SiO₂ (spatial density only)
    # Igneous
    rs_ign = importdataset("output/resampled/rs_ign.tsv", '\t', importas=:Tuple)
    rs_volc = importdataset("output/resampled/rs_volc.tsv", '\t', importas=:Tuple)
    rs_plut = importdataset("output/resampled/rs_plut.tsv", '\t', importas=:Tuple)

    # Sedimentary
    rs_sed = importdataset("output/resampled/rs_sed.tsv", '\t', importas=:Tuple)


## --- Resampled SiO₂ distribution by igneous rock type
    # All igneous
    c, n, = bincounts(rs_ign.SiO2, 40, 80, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="All Igneous (resample); n = $(length(rs_ign.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, color=:purple, 
        linecolor=:purple, ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_ign.png")

    # Volcanic
    c, n, = bincounts(rs_volc.SiO2, 40, 80, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Volcanic (resample); n = $(length(rs_volc.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, color=:purple, 
        linecolor=:purple, ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_volc.png")

    # Plutonic
    c, n, = bincounts(rs_plut.SiO2, 40, 80, 160)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Plutonic (resample); n = $(length(rs_plut.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, color=:purple, 
        linecolor=:purple, ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_plut.png")


## --- Resampled SiO₂ distribution by sedimentary rock type
    # All sedimentary
    c, n, = bincounts(rs_sed.SiO2, 0, 100, 100)
    n = float(n) ./ nansum(float(n) .* step(c))
    h = plot(c, n, seriestype=:bar, label="Sedimentary (resample); n = $(length(rs_sed.SiO2))", 
        ylabel="Weight", xlabel="SiO2 [wt.%]", framestyle=:box, color=:purple, 
        linecolor=:purple, ylims=(0, round(maximum(n), digits=2)+0.01))
    display(h)
    savefig("c_rs_sed.png")


## --- [FUNCTION] Find the most common EarthChem sample in the modal bin
    using Plots

    """
    ```julia
    sameindex(type::Symbol, macro_cats, bulk, bulkidx; hist=:on)
    ```

    Count how many samples in the modal bin are from the same EarthChem sample. Optionally
    plot a histogram of the frequency of all selected indices.

    Data in `bulk` __must__ be the unmatched samples, and `bulkidx` must be filtered so
    there are no indices of 0.
    """
    function sameindex(type::Symbol, macro_cats, bulk, bulktext, bulkidx; 
        bins, hist=:on)

        filter = macro_cats[type]

        c, n, = bincounts(bulk.SiO2[bulkidx][filter], bins...,)

        # Find modal bin and all the points in it
        i = findmax(n)[2]
        s = step(c)/2
        tᵢ = @. c[i]-s <= bulk.SiO2[bulkidx][filter] <= c[i]+s

        # Find indices of those points, and count the frequency
        ind = bulkidx[filter][tᵢ]
        unind = unique(ind)
        counts = [count(==(i), ind) for i in unind]

        # What percent of the indices in this bin are the mode?
        f = findmax(counts)[1] 
        j = unind[findmax(counts)[2]]


        # What percent of total indices are this index?
        totalcount = count(==(j), bulkidx[filter])
        totalindex = length(bulkidx[filter])

        # Terminal printout
        @info """
        Type: $type
        Modal bin: $i ($(c[i]-s)-$(c[i]+s) wt.% SiO₂)
        Modal index count: $f of $(length(ind)) ($(round(f/length(ind)*100, digits=2))%)
        Index: $(j)

        This sample is $(round(totalcount/totalindex*100, digits=2))% of all $type matches.
        ---
        Sample information:

        Age: $(bulk.Age[j])
        Lat, Lon: $(bulk.Latitude[j]), $(bulk.Longitude[j])
        SiO₂: $(round(bulk.SiO2[j], digits=2))%

        Name: $(bulktext.Rock_Name[j])
        Type: $(bulktext.Type[j])
        Material: $(bulktext.Material[j])
        """

        # Histogram
        if hist==:on
            c = sort!(unique(bulkidx[filter]))
            n = [count(==(i), bulkidx[filter]) for i in c]
            h = plot(c, n, seriestype=:bar, label="$type", xlabel="Index [all bins]", 
                ylabel="Frequency", framestyle=:box, ylims=(0, maximum(n)+100), color=:black,
            )
            display(h)
        end
    end


## --- [FN CALL, DATA] sameindex() modal bin index counter
    # Load EarthChem data
    bulkfid = h5open("output/bulk.h5", "r")
    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Load EarthChem metadata
    path = bulkfid["bulktext"]["sampledata"]
    header = read(path["header"])
    index = read(path["index"])
    target = ["Rock_Name", "Type", "Material"]
    targetind = [findall(==(i), header)[1] for i in target]
    bulktext = NamedTuple{Tuple(Symbol.(target))}(
        [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] for i in eachindex(target)]
    )
    close(bulkfid)

    # Igneous
    ignbin = (40, 80, 40)
    sameindex(:ign, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)     # All ign
    sameindex(:volc, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)    # Plutonic 
    sameindex(:plut, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off)    # Volcanic

    # Sedimentary
    sedbin = (0, 100, 100)
    sameindex(:sed, macro_cats, bulk, bulktext, bulkidx[t], bins=sedbin, hist=:off)     # All sed

    # Query Macrostrat at that location:
    # https://macrostrat.org/api/mobile/map_query?lat=LAT&lng=LON&z=11

## --- End of file