## --- Model erosion (slope) as a function of rock age
    # TO DO: correlation plot
    # TO DO: PCA?

## -- Set up
    # Packages
    using StatGeochem
    using HDF5
    using DelimitedFiles

    using LoopVectorization
    using Static
    using Plots
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
    # rock_ersn = emmkyr.(rockslope)
    # ersn, = unmeasurementify(rock_ersn)


## --- Resample matched data by major rock type
    # Set up
    nsims = Int(1e6)
    simitemsout = [macrostrat.rocklat, macrostrat.rocklon, macrostrat.age, rockslope]
    simitemsuncert = (
        zeros(length(macrostrat.rocklat)),
        zeros(length(macrostrat.rocklon)),
        ((macrostrat.agemax .- macrostrat.agemin) ./ 2),
        rockslope_uncert,
    )

    # Samples without min / max bounds are assigned an uncertainty of 5% of the rock age
    for i in eachindex(simitemsuncert[3])
        simitemsuncert[3][i] = ifelse(isnan(simitemsuncert[3][i]), 0.05 * macrostrat.age[i], 
            simitemsuncert[3][i]
        )
    end

    # Preallocate
    simout = (
        sed = Array{Float64}(undef, nsims, length(simitemsout)),
        ign = Array{Float64}(undef, nsims, length(simitemsout)),
        met = Array{Float64}(undef, nsims, length(simitemsout)),
    )

    # Run simulation for each rock type
    simtypes = collect(keys(simout))
    for i in eachindex(simtypes)
        # Get data and uncertainty
        ndata = count(macro_cats[simtypes[i]])

        data = Array{Float64}(undef, ndata, length(simitemsout))
        uncert = Array{Float64}(undef, ndata, length(simitemsout))
        for j in eachindex(simitemsout)
            data[:,j] .= simitemsout[j][macro_cats[simtypes[i]]]
            uncert[:,j] .= simitemsuncert[j][macro_cats[simtypes[i]]]
        end

        # Get resampling weights (spatiotemporal)
        k = invweight(macrostrat.rocklat, macrostrat.rocklon, macrostrat.age)
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

        # Run simulation and save results
        simout[simtypes[i]] .= bsresample(data, uncert, nsims, p)
    end


## --- Plot erosion rate vs. rock age
    # t = @. ersn < 600
    t = trues(length(ersn))

    c,m,e = binmeans(macrostrat.age[t], ersn[t], 0, 3800, 38)
    Plots.plot(c, m, yerror=e, seriestype=:scatter,framestyle=:box, 
        color=:red, lcolor=:red, msc=:red,
        label="", ylabel="Erosion rate [m/Myr]", xlabel="Bedrock Age [Ma]",
        yaxis=log10
    )


## --- Plot erosion rate as a function of slope, but by rock type
    c,m,e = binmeans(simout.sed[:,3], simout.sed[:,4], 0, 3800, 38)
    h1 = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue, framestyle=:box,
        label="Sed", xlabel="Bedrock Age [Ma]", ylabel="Hillslope [m/km]", 
        markershape=:circle, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    c,m,e = binmeans(simout.ign[:,3], simout.ign[:,4], 0, 3800, 38)
    h2 = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red, framestyle=:box,
        label="Ign", xlabel="Bedrock Age [Ma]", ylabel="Hillslope [m/km]", 
        markershape=:circle, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    c,m,e = binmeans(simout.met[:,3], simout.met[:,4], 0, 3800, 38)
    h3 = Plots.plot(c, m, yerror=e, color=:purple, lcolor=:purple, msc=:purple, framestyle=:box,
        label="Met", xlabel="Bedrock Age [Ma]", ylabel="Hillslope [m/km]", 
        markershape=:circle, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    Plots.plot(h1, h2, h3, layout=(3,1), size=(600, 1200))


## --- Everything everywhere over 3800 million years?
    # (Plot everything on the same axis)

    h = Plots.plot(xlabel="Bedrock Age [Ma]", ylabel="Hillslope [m/km]", framestyle=:box,
        legend=:topright, yaxis=log10
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

## --- End of file 