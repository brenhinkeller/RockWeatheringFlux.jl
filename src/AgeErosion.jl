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

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)
    ersn, = unmeasurementify(rock_ersn)


## --- Plot erosion rate vs. slope [????]
    # t = @. ersn < 600
    t = trues(length(ersn))

    c,m,e = binmeans(macrostrat.age[t], ersn[t], 0, 3800, 38)
    plot(c, m, yerror=e, seriestype=:scatter,framestyle=:box, 
        color=:red, lcolor=:red, msc=:red,
        label="", ylabel="Erosion rate [m/Myr]", xlabel="Bedrock Age [Ma]",
        yaxis=log10
    )


## --- Plot erosion rate as a function of slope, but by rock type
    # h = plot(xlabel="Bedrock Age [Ma]", ylabel="Erosion rate [m/Myr]", framestyle=:box,
    #     legend=:topright, yaxis=:log10
    # )

    c,m,e = binmeans(macrostrat.age[macro_cats.sed], ersn[macro_cats.sed], 0, 3800, 38)
    h1 = plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue, framestyle=:box,
        label="Sed", xlabel="Bedrock Age [Ma]", ylabel="Erosion rate [m/Myr]", 
        seriestype=:scatter, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    c,m,e = binmeans(macrostrat.age[macro_cats.ign], ersn[macro_cats.ign], 0, 3800, 38)
    h2 = plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red, framestyle=:box,
        label="Ign", xlabel="Bedrock Age [Ma]", ylabel="Erosion rate [m/Myr]", 
        seriestype=:scatter, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    c,m,e = binmeans(macrostrat.age[macro_cats.met], ersn[macro_cats.met], 0, 3800, 38)
    h3 = plot(c, m, yerror=e, color=:purple, lcolor=:purple, msc=:purple, framestyle=:box,
        label="Met", xlabel="Bedrock Age [Ma]", ylabel="Erosion rate [m/Myr]", 
        seriestype=:scatter, yaxis=:log10, legend=:topright, # ylims=(10,500)
    )

    plot(h1, h2, h3, layout=(3,1), size=(600, 1200))

## --- End of file 