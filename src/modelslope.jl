## --- Set up
    # Packages
    using StatGeochem
    using ProgressMeter
    using Statistics
    using LsqFit: curve_fit
    using HDF5
    using Measurements
    using DelimitedFiles

    # Local utilities
    include("utilities_slope.jl")


## --- Get OCTOPUS basin polygon outlines
    @info "Loading OCTOPUS data"

    # Decompress, load, and recompress 
    run(`gunzip data/octopus/crn_basins_global.kml.gz`);
    (str, isfirstcoord, nbasins, subbasins) = load_octopus_kml("data/octopus/crn_basins_global.kml");
    run(`gzip data/octopus/crn_basins_global.kml`);                                           

    # Parse and export data
    octopusdata = parse_octopus_kml_variables(str)
    exportdataset(octopusdata,"output/octopusdata.tsv",'\t')

    for i in keys(octopusdata)
        if typeof(octopusdata[i]) == Vector{SubString{String}}
            octopusdata[i] = string.(octopusdata[i])
        end
    end
    octopusdata = NamedTuple{Tuple(Symbol.(keys(octopusdata)))}(values(octopusdata))

    # Get basin polygons
    (basin_polygon_n, basin_polygon_lat, basin_polygon_lon) = parse_octopus_polygon_outlines(str,isfirstcoord)


## --- Calculate slope for each basin
    @info "Calculating slope for each basin"

    # Load and parse data 
    # srtm = h5open("data/srtm15plus_aveslope.h5", "r")
    srtm = h5open("data/srtm15plus_maxslope.h5", "r")
    srtm = read(srtm["vars"])
    srtm = NamedTuple{Tuple(Symbol.(keys(srtm)))}(values(srtm))

    # Expected runtime 1.5 hrs
    @timev avgslope, stdslope = get_basin_srtm15plus_aveslope(srtm, nbasins, subbasins, 
        basin_polygon_lat, basin_polygon_lon
    )

    # Save file
    header = ["avg_slope" "err"]
    writedlm("data/basin_srtm15plus_aveslope.tsv", vcat(header, hcat(avgslope, stdslope)))


## --- Alternatively, load pregenerated slope data for each basin
    @info "Loading basin slope data"
    basin_srtm = importdataset("data/basin_srtm15plus_aveslope.tsv", '\t', importas=:Tuple)

    
## --- Fit raw erosion rate as a function of slope (m/km)
    @info "Fitting erosion / slope curve"
        
    # Be data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(octopusdata.ebe_mmkyr) .& (basin_srtm.avg_slope .< 300)
    x = basin_srtm.avg_slope[t]
    y = log10.(octopusdata.ebe_mmkyr[t])

    # Al data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(octopusdata.eal_mmkyr) .& (basin_srtm.avg_slope .< 300)
    x = append!(x, basin_srtm.avg_slope[t])
    y = append!(y, log10.(octopusdata.eal_mmkyr[t]))

    # Fit curve
    p = [0.5, 1/100]
    fobj = curve_fit(linear, x, y, p)

    function emmkyr(slp)
        return 10^(slp * (fobj.param[2]) + (fobj.param[1]))
    end


## --- Plot results
    # # Ideally instead of a line the model would be a ribbon to show error values
    # h = scatter(basin_srtm15plus_aveslope,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", msc=:auto, color=:blue)
    # scatter!(h, basin_srtm15plus_aveslope,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", msc=:auto, color=:orange)
    # plot!(h, 1:500, emmkyr.(1:500), label = "E = 10^(0.056) + 1.00", width=3, color=:black, fg_color_legend=:white)
    # plot!(xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, legend=:topleft, framestyle=:box)
    # savefig(h, "slopeerosion.pdf")


## --- Calculate precipitation for each basin
    # @info "Calculating basin precipitation"

    # prate = ncread("data/prate.sfc.mon.ltm.nc","prate")
    # lat = ncread("data/prate.sfc.mon.ltm.nc","lat")
    # lon = ncread("data/prate.sfc.mon.ltm.nc","lon")

    # precip = mean(prate, dims=3)[:,:,1]
    # imsc(collect(precip'), viridis)

    # data["precip"] = find_precip(data["y_wgs84"],data["x_wgs84"])

    # # Test find_precip to make sure it's not flipped
    # lat = repeat(-90:90,1,361)
    # lon = repeat((-180:180)',181,1)
    # imsc(find_precip(lat,lon),viridis)

## --- End of file