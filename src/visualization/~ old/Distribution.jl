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
    using StatsPlots

    # Local utilities
    include("../utilities/Utilities.jl")

    # More definitions
    bins = (ign = (40,80,40), sed = (0,100,100), met = (25,100,75))   # xmin, xmax, nbins
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    geochemkeys, = get_elements()


## --- [DATA] Matched EarthChem and Macrostrat samples
    data = readdlm("$matchedbulk_io")

    # Get indicies of matched samples
    bulkidx = Int.(vec(data[:,1]))
    t = @. bulkidx != 0

    wtype = string.(vec(data[:,2]))
    macro_cats = match_rocktype(wtype[t])

    # Make matches nice and inclusive
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

    # Matched Macrostrat samples
    macrofid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocktype = read(macrofid["vars"]["rocktype"])[t],
        rockname = read(macrofid["vars"]["rockname"])[t],
        rockdescrip = read(macrofid["vars"]["rockdescrip"])[t],
        rocklat = read(macrofid["vars"]["rocklat"])[t],
        rocklon = read(macrofid["vars"]["rocklon"])[t],
        age = read(macrofid["vars"]["age"])[t],
    )

    # Type matches are from the list of types in the index file now
    # header = read(macrofid["type"]["macro_cats_head"])
    # data = read(macrofid["type"]["macro_cats"])
    # data = @. data > 0
    # macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])

    # for type in minorsed
    #     macro_cats.sed .|= macro_cats[type]
    # end
    # for type in minorign
    #     macro_cats.ign .|= macro_cats[type]
    # end
    # for type in minormet
    #     macro_cats.met .|= macro_cats[type]
    # end

    close(macrofid)

    # Earthchem data
    bulkfid = h5open("output/bulk.h5", "r")

    header = read(bulkfid["bulk"]["header"])
    data = read(bulkfid["bulk"]["data"])

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

    # Data
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Matches
    header = read(bulkfid["bulktypes"]["bulk_cats_head"])
    data = read(bulkfid["bulktypes"]["bulk_cats"])
    data = @. data > 0

    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    mbulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])

    close(bulkfid)


## --- SiO₂ distribution by rock type
    # Except... I don't want 30 plots. I want to plot the sedimentary rocks where silica
    # distributions are useful, and volcanic / plutonic / igneous rocks. Metamorphic rocks
    # have been reassigned so no more data there.
    get_visualized = (:siliciclast, :shale, :carb, :chert, :sed, :volc, :plut, :ign)
    
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(get_visualized)) 

    for i in eachindex(get_visualized)
        # Figure out what kind of bins we want
        r = get_visualized[i]
        if r==:sed || r in minorsed
            nb = bins.sed
        elseif r==:ign || r in minorign
            nb = bins.ign
        end

        # Lol get visualized 
        c, n = bincounts(mbulk.SiO2[macro_cats[r]], nb...)
        n = float(n) ./ nansum(float(n) .* step(c))
        h = plot(c, n, seriestype=:bar, framestyle=:box, color=colors[r], linecolor=colors[r],
            title="$(string(r)); n = $(count(macro_cats[r]))", label="",     
            ylims=(0, round(maximum(n), digits=2)+0.01) 
        )
        fig[i] = h
    end

    h = plot(fig..., layout=(3,3), size=(1800, 1200))
    display(h)
    # savefig(h, "results/figures/distributions.png")


## --- [FN CALL] Igneous modal index and Macrostrat data
    # Igneous
    iᵢ = sameindex(bulkidx[t][macro_cats.ign], geochemkeys, bins.ign, bulk, bulktext)
    iᵥ = sameindex(bulkidx[t][macro_cats.volc], geochemkeys, bins.ign, bulk, bulktext)
    iₚ = sameindex(bulkidx[t][macro_cats.plut], geochemkeys, bins.ign, bulk, bulktext)

    # get_matched_samples(iᵢ, bulkidx[t], macrostrat, filter=macro_cats.ign, desc="Igneous")
    # get_matched_samples(iᵥ, bulkidx[t], macrostrat, filter=macro_cats.volc, desc="Volcanic")
    # get_matched_samples(iₚ, bulkidx[t], macrostrat, filter=macro_cats.plut, desc="Plutonic")

    # Query Macrostrat at that location:
    # https://macrostrat.org/api/mobile/map_query?lat=LAT&lng=LON&z=11


## --- [FN CALL] Sedimentary + metamorphic modal index
    # Sedimentary
    iₛ = sameindex(bulkidx[t][macro_cats.sed], geochemkeys, bins.sed, bulk, bulktext)
    # get_matched_samples(iₛ, bulkidx[t], macrostrat, filter=macro_cats.sed, desc="Sedimentary")

    sameindex(bulkidx[t][macro_cats.shale], geochemkeys, bins.sed, bulk, bulktext)
    sameindex(bulkidx[t][macro_cats.carb], geochemkeys, bins.sed, bulk, bulktext)
    sameindex(bulkidx[t][macro_cats.chert], geochemkeys, bins.sed, bulk, bulktext)
    sameindex(bulkidx[t][macro_cats.evaporite], geochemkeys, bins.sed, bulk, bulktext)
    sameindex(bulkidx[t][macro_cats.coal], geochemkeys, bins.sed, bulk, bulktext)

    # Metamorphic
    sameindex(bulkidx[t][macro_cats.met], geochemkeys, bins.met, bulk, bulktext)
    

## --- Igneous SiO₂ by rock age
    t = @. (mbulk.Age >= 2500) & macro_cats.ign
    c, n = bincounts(mbulk.SiO2[t], bins.ign...)
    n = float(n) ./ nansum(float(n) .* step(c))
    h1 = plot(c, n, seriestype=:bar, framestyle=:box, color=colors.ign, linecolor=colors.ign,
        label="Archean; n = $(count(t))",      
        # ylabel="Weight", xlabel="SiO2 [wt.%]",
        ylims=(0, round(maximum(n), digits=2)+0.01) 
    )

    t = @. (2500 > mbulk.Age >= 541) & macro_cats.ign
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
    savefig(h, "results/figures/ign_ages.png")


## --- Alternate (Stacked) representation of igneous SiO₂ distribution
    # Important!! The distributions are normalized **before** getting sent to the stacked 
    # bar chart, meaning that this isn't the actual distribution of igneous silica: this
    # shows the relative distribution of silica **for each age distribution**. This is
    # necessary because there are so few Archean and Proterozoic rocks compared to 
    # Phanerozoic rocks. Great Unconformity strikes again.

    t = @. (mbulk.Age >= 2500) & macro_cats.ign
    c, nₐ = bincounts(mbulk.SiO2[t], bins.ign...)
    nₐ = float(nₐ) ./ nansum(float(nₐ) .* step(c))

    t = @. (2500 > mbulk.Age >= 541) & macro_cats.ign
    c, nᵣ = bincounts(mbulk.SiO2[t], bins.ign...)
    nᵣ = float(nᵣ) ./ nansum(float(nᵣ) .* step(c))

    t = @. (541 > mbulk.Age) & macro_cats.ign
    c, nₚ = bincounts(mbulk.SiO2[t], bins.ign...)
    nₚ = float(nₚ) ./ nansum(float(nₚ) .* step(c))

    ticklabel = bins.ign[1]:5:bins.ign[2]
    ticklabel = string.(collect(ticklabel))
    h = groupedbar([nₐ nᵣ nₚ], bar_position=:stack, bar_width=1.0,
        xticks=(1:5:length(c)+1, ticklabel), framestyle=:box, legend=:topright,
        xlabel="SiO₂ [wt.%]", ylabel="Weight", label=["Archean" "Proterozoic" "Phanerozoic"],
        ylims=(0, round(maximum(nₐ.+nᵣ.+nₚ), digits=2)+0.01),
        color=["firebrick" "seagreen" "cyan"],
    )


## --- Spatially resampled EarthChem SiO₂
    # Load file and figure out where the silica is 
    fid = h5open("output/resampled/resampled.h5", "r")
    header = read(fid["vars"]["header"])
    SiO₂ᵢ = findfirst(x -> x=="SiO2", header)

    #  Preallocate
    plts = Array{Plots.Plot{Plots.GRBackend}}(undef, length(macro_cats) - 1)
    bsrbins = (sed=(0,100,200), ign=(40,80,80), met=(25,100,150))   # min, max, nbins

    # Create all plots
    rocks = collect(keys(macro_cats))
    for i in eachindex(rocks)
        # Figure out what type of rock we're looking at
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

        # Get silica data out of the resampled dataset
        silica = read(fid["vars"]["data"]["$r"]["data"])[:,SiO₂ᵢ]

        c, n = bincounts(silica, bsrbins[type]...)
        n = float(n) ./ nansum(float(n) .* step(c))
        h = plot(c, n, seriestype=:bar, framestyle=:box, color=c_rs, linecolor=c_rs,
            label="$(string(r)); n = $(length(silica))",      
            # ylabel="Weight", xlabel="SiO2 [wt.%]",
            ylims=(0, round(maximum(n), digits=2)+0.01) 
        )

        plts[i] = h
    end

    # Put into a layout
    h = plot(plts..., layout=(5,3), size=(2000, 2000))
    display(h)
    savefig(h, "results/figures/resampled_distributions.png")


## --- Major element geochemistry of a given type
    # TEMPORARY: assumes stuff has already been loaded in SampleMatch.jl
    # Also working on the looping stuff, so variable names will need to be changed if this 
    # sticks around

    # https://stackoverflow.com/questions/49168872/increase-space-between-subplots-in-plots-jl

    # type = :shale
    # nplots = length(geochemcopy)
    # plts = Array{Plots.Plot{Plots.GRBackend}}(undef, nplots)

    # for i in eachindex(geochemcopy)
    #     c, n = bincounts(bulk[geochemcopy[i]][bulk_cats[type]], 0, 100, 100)
    #     n = float(n) ./ nansum(float(n) .* step(c))
    #     h = plot(c, n, seriestype=:bar, framestyle=:box, color=:black, linecolor=:black,
    #         label="", ylabel="Weight", xlabel="$(geochemcopy[i]) [wt.%]",
    #         ylims=(0, round(maximum(n), digits=2)+0.01) 
    #     )

    #     plts[i] = h
    # end

    # h = plot(plts..., layout=(nplots, 1), size=(700, nplots*450), 
    #     bottom_margin=(25, :px), left_margin=(150, :px)
    # )
    # display(h)


## --- Distance from matched sample (histogram)
    # TEMPORARY: assumes stuff has already been loaded in SampleMatch.jl
    # dist = haversine.(macrostrat.rocklat[filter], macrostrat.rocklon[filter],
    #     bulk.Latitude[matches], bulk.Longitude[matches]
    # )
    # c, n = bincounts(dist, 0, 180, 90)
    # h = plot(c, n, seriestype=:bar, framestyle=:box,
    #     label="", ylabel="Count", xlabel="Distance [arc degrees]",
    #     ylims=(0, maximum(n)+50) 
    # )


## --- [DATA] Composition of eroded material
    abs = importdataset("$erodedabs_out", ',', importas=:Tuple)
    rel = importdataset("$erodedrel_out", ',', importas=:Tuple)

    # Grab only the data for the extended major elements (majors + P₂O₅)
    majors, minors = get_elements()
    modmajors = [majors; :P2O5]
    row = Array{Int64}(undef, length(modmajors))
    rowname = Symbol.(abs[1])
    for i in eachindex(modmajors)
        for j in eachindex(rowname)
            if modmajors[i]==rowname[j]
                row[i] = j
                break
            end
        end
    end


## --- Composition of eroded material
    ticklabel = string.(modmajors)

    # Absolute
    h = groupedbar([abs.sed[row] abs.ign[row] abs.met[row]], bar_position=:stack, 
        bar_width=0.85, framestyle=:box, legend=:topright, 
        xticks=(1:length(ticklabel), ticklabel), 
        xlabel="Major oxide", ylabel="Absolute Denudation [Gt/yr]", 
        label = ["Sedimentary" "Igneous" "Metamorphic"],
        ylims=(0, round(maximum(abs.sed[row] .+ abs.ign[row] .+ abs.met[row]), digits=2)+1),
        color = [colors.sed colors.ign colors.met]
    )
    savefig(h, "results/figures/absdenudation.png")

    # Relative
    h = groupedbar([rel.sed[row] rel.ign[row] rel.met[row]], bar_position=:stack, 
        bar_width=0.85, framestyle=:box, legend=:topright, 
        xticks=(1:length(ticklabel), ticklabel), 
        xlabel="Major oxide", ylabel="Relative Denudation [Gt/yr]", 
        label = ["Sedimentary" "Igneous" "Metamorphic"],
        ylims=(0, 1.3),
        color = [colors.sed colors.ign colors.met]
    )
    savefig(h, "results/figures/reldenudation.png")


## --- End of file