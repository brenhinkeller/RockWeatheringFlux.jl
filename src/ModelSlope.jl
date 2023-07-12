## --- Set up
    # Packages
    using StatGeochem
    using Statistics
    using HDF5
    using Measurements
    using DelimitedFiles
    using Plots
    using NetCDF
    using Isoplot: yorkfit

    # Local utilities
    include("utilities/Utilities.jl")


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


## --- Calculate slope for each basin
    # @info "Calculating slope for each basin"

    # # Get basin polygons
    # (basin_polygon_n, basin_polygon_lat, basin_polygon_lon) = parse_octopus_polygon_outlines(str,isfirstcoord)

    # # Load and parse SRTM15+ data 
    # # srtm = h5open("output/srtm15plus_aveslope.h5", "r")
    # srtm = h5open("output/srtm15plus_maxslope.h5", "r")
    # srtm = read(srtm["vars"])
    # srtm = NamedTuple{Tuple(Symbol.(keys(srtm)))}(values(srtm))

    # # Expected runtime 1.5 hrs
    # @timev avgslope, stdslope = get_basin_srtm15plus_aveslope(srtm, nbasins, subbasins, 
    #     basin_polygon_lat, basin_polygon_lon
    # )

    # # Save file
    # header = ["avg_slope" "err"]
    # writedlm("output/basin_srtm15plus_avg_maxslope.tsv", vcat(header, hcat(avgslope, stdslope)))

    # # File names:
    # #     basin_srtm15plus_avg_maxslope -> average of maximum slopes of each point in basin
    # #     basin_srtm15plus_aveslope     -> average of average slopes of each point in the basin


## --- Alternatively, load pregenerated OCTOPUS and slope data for each basin
    @info "Loading pre-parsed OCTOPUS and basin slope data"
    octopusdata = importdataset("output/octopusdata.tsv",'\t', importas=:Tuple)
    basin_srtm = importdataset("output/basin_srtm15plus_avg_maxslope.tsv", '\t', importas=:Tuple)


## --- Fit erosion rate / slope curve
    # Collect erosion rate and slope data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(basin_srtm.err)
    t_be = t .& .!isnan.(octopusdata.ebe_mmkyr) .& .!isnan.(octopusdata.ebe_err)
    t_al = t .& .!isnan.(octopusdata.eal_mmkyr) .& .!isnan.(octopusdata.eal_err)

    x = (
        v = [basin_srtm.avg_slope[t_be]; basin_srtm.avg_slope[t_al]],
        e = [basin_srtm.err[t_be]; basin_srtm.err[t_al]]
    )
    y = (
        v = [octopusdata.ebe_mmkyr[t_be]; octopusdata.eal_mmkyr[t_al]],
        e = [octopusdata.ebe_err[t_be]; octopusdata.eal_err[t_al]]
    )

    # Get mean in regular space for bins with equal numbers of points
    c, m, ex, ey = binmeans_percentile(x.v, y.v, step=5)

    # Fit slope to means
    fobj = yorkfit(c, ex, log.(m), log.(ey))
    emmkyr(slp) = exp(slp * (fobj.slope) + (fobj.intercept))

    # Plot results
    h = scatter(basin_srtm.avg_slope,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", 
        msc=:auto, color=:blue, alpha=0.5
    )
    scatter!(h, basin_srtm.avg_slope,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", 
        msc=:auto, color=:orange, alpha=0.5
    )
    scatter!(h, c, m, label="Binned Means", msc=:auto, color=:black)

    # Model
    modelval, modelerr = unmeasurementify(emmkyr.(1:600))
    plot!(h, 1:length(modelval), modelval, label="Model", color=:black, width=3)
    plot!(xlabel="SRTM15+ Slope (m/km)", ylabel="Erosion rate (mm/kyr)",
        yscale=:log10, framestyle=:box, legend=:topleft, fg_color_legend=:white,
        
    )
    display(h)


## --- Model erosion as the product of slope and precipitation
    # Get basin precipitation
    precip = ncread("data/prate.sfc.mon.ltm.nc","prate")
    precip = mean(precip, dims=3)[:,:,1]
    basin_precip = find_precip(octopusdata.y_wgs84, octopusdata.x_wgs84, precip)

    # Get slope ⋅ precipitation
    slopeprecip = (basin_srtm.avg_slope .± basin_srtm.err) .* basin_precip
    v, e = unmeasurementify(slopeprecip)
    slopeprecip = (v=v, e=e)

    # Get erosion
    t = .!isnan.(slopeprecip.v) .& .!isnan.(slopeprecip.e)
    t_be = t .& .!isnan.(octopusdata.ebe_mmkyr) .& .!isnan.(octopusdata.ebe_err)
    t_al = t .& .!isnan.(octopusdata.eal_mmkyr) .& .!isnan.(octopusdata.eal_err)

    x = (
        v = [slopeprecip.v[t_be]; slopeprecip.v[t_al]],
        e = [slopeprecip.e[t_be]; slopeprecip.e[t_al]]
    )
    y = (
        v = [octopusdata.ebe_mmkyr[t_be]; octopusdata.eal_mmkyr[t_al]],
        e = [octopusdata.ebe_err[t_be]; octopusdata.eal_err[t_al]]
    )

    # Get bin averages in regular space
    c, m, ex, ey = binmeans_percentile(x.v, y.v, step=5)

    # Log transform **both** x and y and fit model
    fobj = yorkfit(log.(c), log.(ex), log.(m), log.(ey))
    # emmkyr_precip(slp) = exp(slp * (fobj.slope) + (fobj.intercept))

    # Plot results
    h = scatter(slopeprecip.v,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", 
        msc=:auto, color=:blue, alpha=0.5
    )
    scatter!(h, slopeprecip.v,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", 
        msc=:auto, color=:orange, alpha=0.5
    )
    scatter!(h, c, m, label="Binned Means", msc=:auto, color=:black)

    # Model
    modelin = range(start=0, stop=0.4, length=100)
    modelval, modelerr = unmeasurementify(emmkyr_precip.(modelin))
    # plot!(h, 1:length(modelval), modelval, label="Model", color=:black, width=3)
    plot!(xlabel="SRTM15+ Slope (m/km)", ylabel="Erosion rate (mm/kyr)",
        xscale=:log10, yscale=:log10, framestyle=:box, legend=:topleft, 
        fg_color_legend=:white,
    )
    display(h)
    

## --- End of file