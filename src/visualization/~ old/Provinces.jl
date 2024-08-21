#

## --- Set up
    # Packages
    using HDF5
    using ProgressMeter

    using CairoMakie
    using GeoMakie
    using ImageMagick

    # Local utilities
    using Static
    using LoopVectorization
    using Measurements

    include("../utilities/Utilities.jl")


## --- Load geologic provinces
    geolprov = h5read("data/geolprov.h5", "geolprov")

    # Get the lat, lon coordinates for each index
    nrows, ncols = size(geolprov)
    sf_lat = (nrows - 1) / 180
    sf_lon = (ncols - 1) / 360
    lats = [90 - (1/sf_lat*(i-1)) for i = 1:nrows]
    lons = [180 - (1/sf_lon*(i-1)) for i = 1:ncols]

## --- ... Plot geologic provinces



## ---
    # Define colors
    provcolors = (
        10 => "Accreted_Arc",
        11 => "Island_Arc",
        12 => "Continental_Arc",
        13 => "Collisional_Orogen",
        20 => "Extensional",
        21 => "Rift",
        22 => "Plume",
        31 => "Shield",
        32 => "Platform",
        33 => "Basin",
        00 => "No data"
    )

    # Make plot
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")

    @showprogress for i in eachindex(lats)
        for j in eachindex(lons)
            # Don't plot anything if there's no data
            if geolprov[i,j] == 0
                continue
            end

            # Plot a point if there is data
            CairoMakie.scatter!(ax, lats[i], lons[j], color=geolprov[i,j], markersize=3,)

        end
    end
    display(f)


## --- Slope as a function of geologic province
    prov = find_geolprov(macrostrat.rocklat, macrostrat.rocklon)

    # This is a little cursed but it's easy to debug
    provkey = (
        10 => "Accreted_Arc",
        11 => "Island_Arc",
        12 => "Continental_Arc",
        13 => "Collisional_Orogen",
        20 => "Extensional",
        21 => "Rift",
        22 => "Plume",
        31 => "Shield",
        32 => "Platform",
        33 => "Basin",
        # 00 => "No_Data",
    )
    provcodes = [values(provkey)[i].first for i in eachindex(provkey)]
    provnames = [values(provkey)[i].second for i in eachindex(provkey)]
    sample_in_prov = NamedTuple{Tuple(Symbol.(provnames))}([prov .== c for c in provcodes])


## --- Average slope by province
    avg_slope = [rockslope[t] for t in sample_in_prov]


## --- Contribution to top quartile of slope
    rockslope75 = percentile(log10.(rockslope), 75)
    t = @. log10.(rockslope) > rockslope75

    presence = [count(t .& s) / count(t) for s in sample_in_prov]
    x = 1:length(presence)
    h = plot(x, presence, seriestype=:bar, framestyle=:box, label="", 
        ylabel="Relative Abundance in Top Quartile", xlabel="Geologic Province", 
        xticks=(x, provnames), xrotation = 45, 
        ylims = (0, maximum(presence) + 0.1*maximum(presence))
    )
    display(h)


## --- Various experiments
    # Get erosion by province, plot
    ersn_by_prov = [nansum(rock_ersn[t]) for t in sample_in_prov]

    x = 1:length(ersn_by_prov)
    h = plot(x, ersn_by_prov, seriestype=:bar, framestyle=:box, label="", 
        ylabel="Total Erosion [m/Myr]", xlabel="Geologic Province", xticks=(x, provnames), 
        xrotation = 45, ylims = (0, maximum(ersn_by_prov) + 0.1*maximum(ersn_by_prov))
    )
    display(h)

    # Normalize erosive contribution by relative abundance
    abundance_by_prov = [count(sample_in_prov[i]) for i in keys(sample_in_prov)]
    x = 1:length(ersn_by_prov)
    h = plot(x, ersn_by_prov./abundance_by_prov, seriestype=:bar, framestyle=:box, label="", 
        ylabel="Normalized Total Erosion", xlabel="Geologic Province", xticks=(x, provnames), 
        xrotation = 45, ylims = (0, 350)
    )
    display(h)

    # The same plot, but consider means ± standard deviations instead of sums 
    # The error bars aren't working
    mean_slp_prov = [nanmean(rockslope[t]) for t in sample_in_prov]
    upper_slp = [percentile(rockslope[t], 95) for t in sample_in_prov] .- mean_slp_prov
    lower_slp = mean_slp_prov .- [percentile(rockslope[t], 5) for t in sample_in_prov]

    x = 1:length(ersn_by_prov)
    h = plot(x, mean_slp_prov, seriestype=:bar, framestyle=:box, label="", 
        ylabel="Average Slope [m/km]", xlabel="Geologic Province",
        xticks=(x, provnames), xrotation = 45,
        yerror=(lower_slp, upper_slp), ylims = (0, 500)
    )
## --- End of file 