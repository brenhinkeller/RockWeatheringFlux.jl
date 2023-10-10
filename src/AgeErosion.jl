## --- Model erosion (slope) as a function of rock age
    # TO DO: correlation plot
    # TO DO: PCA?

## -- Set up
    # Packages
    using StatGeochem
    using HDF5
    using DelimitedFiles
    using StatsBase
    using Plots
    using StatsPlots

    using LoopVectorization
    using Static
    using Measurements

    # Local utilities
    include("utilities/Utilities.jl")


## --- Load Macrostrat data
    # Indices of matched EarthChem samples from SampleMatch.jl
    fid = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    # Matched types, majors inclusive of minors
    bulktype = string.(vec(fid[:,2]))
    macro_cats = match_rocktype(bulktype[t])

    # Macrostrat data
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
    )
    close(fid)


## --- Calculate erosion rate at each point of interest
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point with a known EarthChem sample
    # Modify this function to return an error as well
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, n=5
    )

    rockslope, rockslope_uncert = unmeasurementify(rockslope)

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)
    rock_ersn, rock_ersn_uncert = unmeasurementify(rock_ersn)


## --- Resample matched data by major rock type
    # # Set up
    # nsims = Int(1e6)
    # simitemsout = [macrostrat.rocklat, macrostrat.rocklon, macrostrat.age, rockslope]
    # simitemsuncert = (
    #     zeros(length(macrostrat.rocklat)),
    #     zeros(length(macrostrat.rocklon)),
    #     ((macrostrat.agemax .- macrostrat.agemin) ./ 2),
    #     rockslope_uncert,
    # )

    # # Samples without min / max bounds are assigned an uncertainty of 5% of the rock age
    # for i in eachindex(simitemsuncert[3])
    #     simitemsuncert[3][i] = ifelse(isnan(simitemsuncert[3][i]), 0.05 * macrostrat.age[i], 
    #         simitemsuncert[3][i]
    #     )
    # end

    # # Preallocate
    # simout = (
    #     sed = Array{Float64}(undef, nsims, length(simitemsout)),
    #     ign = Array{Float64}(undef, nsims, length(simitemsout)),
    #     met = Array{Float64}(undef, nsims, length(simitemsout)),
    # )

    # # Run simulation for each rock type
    # simtypes = collect(keys(simout))
    # for i in eachindex(simtypes)
    #     # Get data and uncertainty
    #     datafilter = macro_cats[simtypes[i]]
    #     ndata = count(datafilter)

    #     data = Array{Float64}(undef, ndata, length(simitemsout))
    #     uncert = Array{Float64}(undef, ndata, length(simitemsout))
    #     for j in eachindex(simitemsout)
    #         data[:,j] .= simitemsout[j][datafilter]
    #         uncert[:,j] .= simitemsuncert[j][datafilter]
    #     end

    #     test = @. (!isnan(macrostrat.rocklon[datafilter]) & 
    #         !isnan(macrostrat.rocklat[datafilter]) & !isnan(macrostrat.age[datafilter]))
    #     data = data[test[:],:]
    #     uncert = uncert[test[:],:]

    #     # Get resampling weights (spatiotemporal)
    #     k = invweight(macrostrat.rocklat, macrostrat.rocklon, macrostrat.age)
    #     p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    #     # Run simulation and save results
    #     simout[simtypes[i]] .= bsresample(data, uncert, nsims, p)
    # end

    # # Save data
    # fid = h5open("output/resampled/age_slope.h5", "w")
    #     g = create_group(fid, "vars")
    #     g["sed"] = simout.sed
    #     g["ign"] = simout.ign
    #     g["met"] = simout.met
    # close(fid)


## --- If you already have data, load from file
    fid = h5open("output/resampled/age_slope.h5", "r")
        simout = (
            sed = read(fid["vars"]["sed"]),
            ign = read(fid["vars"]["ign"]),
            met = read(fid["vars"]["met"]),
        )
    close(fid)
    

## --- Plot erosion rate as a function of slope, but by rock type
    # Column 3: age
    # Column 4: slope

    c,m,e = binmeans(simout.sed[:,3], simout.sed[:,4], 0, 3800, 38)
    h1 = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue, framestyle=:box,
        label="Sed",
        markershape=:circle, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    c,m,e = binmeans(simout.ign[:,3], simout.ign[:,4], 0, 3800, 38)
    h2 = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red, framestyle=:box,
        label="Ign", ylabel="Hillslope [m/km]", 
        markershape=:circle, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    c,m,e = binmeans(simout.met[:,3], simout.met[:,4], 0, 3800, 38)
    h3 = Plots.plot(c, m, yerror=e, color=:purple, lcolor=:purple, msc=:purple, framestyle=:box,
        label="Met", xlabel="Bedrock Age [Ma]",
        markershape=:circle, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    h = Plots.plot(h1, h2, h3, layout=(3,1), size=(600, 1200), left_margin=(30, :px))
    display(h)
    savefig(h, "results/figures/ageslope.png")


## --- Everything everywhere over 3800 million years?
    # (Plot everything on the same axis)

    h = Plots.plot(xlabel="Bedrock Age [Ma]", ylabel="Hillslope [m/km]", framestyle=:box,
        legend=:topright, yaxis=:log10
    )

    c,m,e = binmeans(simout.sed[:,3], simout.sed[:,4], 0, 3800, 38)
    Plots.plot!(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
        label="Sed", markershape=:diamond,
    )

    c,m,e = binmeans(simout.ign[:,3], simout.ign[:,4], 0, 3800, 38)
    Plots.plot!(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
        label="Ign", markershape=:circle,
    )

    c,m,e = binmeans(simout.met[:,3], simout.met[:,4], 0, 3800, 38)
    Plots.plot!(c, m, yerror=e, color=:purple, lcolor=:purple, msc=:purple,
        label="Met", markershape=:star5,
    )

    display(h)
    savefig(h, "results/figures/ageslope_stack.png")


## --- Try to figure out why Archean samples are eroding so fast
## --- Show that this exists in the real data and is not an artifact of resampling
    c,m,e = binmeans(macrostrat.age[macro_cats.sed], rockslope[macro_cats.sed], 0, 3800, 38)
    h1 = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue, framestyle=:box,
        label="Sed",
        markershape=:circle, yaxis=:log10, legend=:topright,
    )

    c,m,e = binmeans(macrostrat.age[macro_cats.ign], rockslope[macro_cats.ign], 0, 3800, 38)
    h2 = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red, framestyle=:box,
        label="Ign", ylabel="Hillslope [m/km]", 
        markershape=:circle, yaxis=:log10, legend=:topright,
    )

    c,m,e = binmeans(macrostrat.age[macro_cats.met], rockslope[macro_cats.met], 0, 3800, 38)
    h3 = Plots.plot(c, m, yerror=e, color=:purple, lcolor=:purple, msc=:purple, framestyle=:box,
        label="Met", xlabel="Bedrock Age [Ma]",
        markershape=:circle, yaxis=:log10, legend=:topright,
    )

    h = Plots.plot(h1, h2, h3, layout=(3,1), size=(600, 1200), left_margin=(30, :px))
    display(h)
    

## --- Compare Eo - Paleoarchean (higher slope) and Meso - Neoarchean (lower slope)
    # Geologic provinces (Most Archean rocks are shields)
    ep_archean = @. macrostrat.age > 3200;
    mn_archean = @. 3200 > macrostrat.age >= 2500;

    ep_prov = decode_find_geolprov(find_geolprov(macrostrat.rocklat[ep_archean], 
        macrostrat.rocklon[ep_archean])
    )
    mn_prov = decode_find_geolprov(find_geolprov(macrostrat.rocklat[mn_archean], 
        macrostrat.rocklon[mn_archean])
    )

    # Normalized distribution of provinces for each age category
    uniqueprovs = unique(archeanprov)
    ep_prov_count = [count(x -> x==name, ep_prov) for name in uniqueprovs] ./ length(ep_prov)
    mn_prov_count = [count(x -> x==name, mn_prov) for name in uniqueprovs] ./ length(mn_prov)

    x = 1:length(uniqueprovs)
    StatsPlots.groupedbar([ep_prov_count mn_prov_count], bar_position=:dodge,
        framestyle=:box, label=["Eo-Paleoarchean" "Meso-Neoarchean"], ylabel="Proportion",
        xlabel="Geologic Province", xticks=(x, uniqueprovs), legend=:topright, 
        bottom_margin=(30, :px), xrotation = 45, ylims = (0, 1),
    )

    # Average slope of Archean rocks?
    archeanslope = rockslope[oldarchean]
    c, n = bincounts(archeanslope, 0, maximum(archeanslope), 15)
    h = Plots.plot(c, n, seriestype=:bar, framestyle=:box,
        label="", ylabel="Abundance", xlabel="Hillslope [m/km]",
        ylims = (0, maximum(n) + 0.1*maximum(n))
    )

## --- Look at just igneous rocks
    c,m,e = binmeans(macrostrat.age[macro_cats.ign], rockslope[macro_cats.ign], 2500, 3800, 13)
    h = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red, framestyle=:box,
        label="Ign", ylabel="Hillslope [m/km]", xlabel="Age [Ma]", seriestype=:scatter,
        markershape=:circle, yaxis=:log10, legend=:topright,
    )


## --- Average slope of each province
    allprov = decode_find_geolprov(find_geolprov(macrostrat.rocklat, macrostrat.rocklon))
    uniqueprovs = unique(allprov)
    inprov = NamedTuple{Tuple(Symbol.(uniqueprovs))}([allprov .== name for name in uniqueprovs])
    
    avg_slope = [nanmean(rockslope[t]) for t in inprov]
    # TO DO: errors aren't behaving as expected--e.g. maximum error isn't plotted at 317 even though
    # that's the maximum upper error :(
    # lower = avg_slope .- [percentile(rockslope[t], 5) for t in inprov]
    # upper = [percentile(rockslope[t], 95) for t in inprov] .- avg_slope

    # Average slope by geologic province
    x = 1:length(avg_slope)
    h = Plots.plot(x, avg_slope, seriestype=:bar, framestyle=:box, label="", 
        ylabel="Average Slope [m/km]", xlabel="Geologic Province", xticks=(x, uniqueprovs), 
        xrotation = 45, ylims = (0, maximum(upper) + 0.1*maximum(upper)),
        bottom_margin=(30, :px)
    )
    display(h)


## --- Abundance of each geologic eon in each province
    # TO DO: normalize to relative abundance relative to abundance of all rock ages
    # Or something to show how Archean rocks tend to be in shields...
    archean = @. macrostrat.age >= 2500;
    proterozoic = @. 541 <= macrostrat.age < 2500;
    phanerozoic = @. macrostrat.age < 541;

    # Total number of rocks by eon
    # count_archean = 
    # count_proterozoic = 
    # count_phanerozoic = 

    # Abundance by province
    arc_provs = [count(x -> x==name, allprov[archean]) for name in uniqueprovs]
    pro_provs = [count(x -> x==name, allprov[proterozoic]) for name in uniqueprovs]
    pha_provs = [count(x -> x==name, allprov[phanerozoic]) for name in uniqueprovs]
    total_provs = @. arc_provs + pro_provs + pha_provs

    x = 1:length(arc_avg_slope)
    h = StatsPlots.groupedbar([arc_provs pro_provs pha_provs], bar_position=:stack,
        framestyle=:box, label=["Archean" "Proterozoic" "Phanerozoic"],
        ylabel="Abundance", xlabel="Geologic Province", xticks=(x, uniqueprovs), 
        xrotation = 45, ylims = (0, maximum(total_provs) + 0.1*maximum(total_provs)),
        legend=:topright, bottom_margin=(30, :px)
    )

    
## --- End of file 