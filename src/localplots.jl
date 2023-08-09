## --- Set up
    # CTRL+K CTR+[ to collapse all sections

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


## --- [DATA] SiO₂ content by igneous rock type
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


## --- [DATA] Resampled EarthChem SiO₂ (spatial density only)
    # Igneous
    rs_ign = importdataset("output/resampled/ign.tsv", '\t', importas=:Tuple)
    rs_volc = importdataset("output/resampled/volc.tsv", '\t', importas=:Tuple)
    rs_plut = importdataset("output/resampled/plut.tsv", '\t', importas=:Tuple)

    # Sedimentary
    rs_sed = importdataset("output/resampled/sed.tsv", '\t', importas=:Tuple)


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


## --- [FUNCTION] Find the most common EarthChem sample in the modal bin
    using Plots

    """
    ```julia
    sameindex(type::Symbol, macro_cats, bulk, bulkidx; hist=:on)
    ```

    Count how many samples in the modal bin are from the same EarthChem sample. Optionally
    plot a histogram of the frequency of all selected indices.
    """
    function sameindex(type::Symbol, macro_cats, bulk, bulkidx; hist=:on)
        filter = macro_cats[type]

        c, n, = bincounts(bulk.SiO2[filter], 40, 80, 160)

        # Find modal bin and all the points in it
        i = findmax(n)[2]
        s = step(c)/2
        s = @. c[i]-s <= bulk.SiO2[filter] <= c[i]+s

        # Find indices of those points, and count the frequency
        ind = bulkidx[filter][s]
        unind = unique(ind)
        counts = [count(==(i), ind) for i in unind]

        # What percent of the indices are the mode?
        f = findmax(counts)[1] / length(ind)
        i = findmax(counts)[2]

        # Terminal printout
        @info "$(round(f*100, digits=2))% of $type indices are from EarthChem sample i = $(unind[i])"

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


## --- [FN CALL] sameindex() modal bin index counter
    # Igneous
    sameindex(:ign, macro_cats, bulk, bulkidx[t], hist=:off)        # All igneous
    sameindex(:volc, macro_cats, bulk, bulkidx[t], hist=:off)       # Plutonic 
    sameindex(:plut, macro_cats, bulk, bulkidx[t], hist=:off)       # Volcanic


## --- Get EarthChem (meta)data for a given sample index
    # Set sample value
    i = 413791

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

    # Terminal printout
    @info """
    Age: $(bulk.Age[i])
    Lat, Lon: ($(bulk.Latitude[i]), $(bulk.Longitude[i]))
    SiO₂: $(round(bulk.SiO2[i], digits=2))%

    Name: $(bulktext.Rock_Name[i])
    Type: $(bulktext.Type[i])
    Material: $(bulktext.Material[i])
    """

    # Query Macrostrat at that location:
    # https://macrostrat.org/api/mobile/map_query?lat=LAT&lng=LON&z=11


## --- End of File