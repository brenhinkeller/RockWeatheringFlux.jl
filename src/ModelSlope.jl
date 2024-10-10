## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, Plots
    using Isoplot: yorkfit


## --- Get OCTOPUS basin polygon outlines [Linux / Mac OS only]
    # @info "Loading OCTOPUS data"

    # # Decompress, load, and recompress 
    # run(`gunzip data/octopus/crn_basins_global.kml.gz`);
    # (str, isfirstcoord, nbasins, subbasins) = load_octopus_kml("data/octopus/crn_basins_global.kml");
    # run(`gzip data/octopus/crn_basins_global.kml`);                                           

    # # Parse and export data
    # octopusdata = parse_octopus_kml_variables(str)
    # exportdataset(octopusdata, octopusdata_basins, '\t')

    # for i in keys(octopusdata)
    #     if typeof(octopusdata[i]) == Vector{SubString{String}}
    #         octopusdata[i] = string.(octopusdata[i])
    #     end
    # end
    # octopusdata = NamedTuple{Tuple(Symbol.(keys(octopusdata)))}(values(octopusdata))


## --- Calculate slope for each basin [Linux / Mac OS only]
    # @info "Calculating slope for each basin"

    # # Get basin polygons
    # (basin_polygon_n, basin_polygon_lat, basin_polygon_lon) = parse_octopus_polygon_outlines(str,isfirstcoord)

    # # Load and parse SRTM15+ maximum slope data 
    # srtm = h5open(srtm_maxslope, "r")
    # srtm = read(srtm["vars"])
    # srtm = NamedTuple{Tuple(Symbol.(keys(srtm)))}(values(srtm))

    # # Grab all maximum slopes in each basin. Average them
    # # We're going to use this to calibrate our erosion / slope relationship
    # # Expected runtime 1.5 hrs
    # @timev avgslope, stdslope = get_basin_srtm15plus_aveslope(srtm, nbasins, subbasins, 
    #     basin_polygon_lat, basin_polygon_lon
    # )

    # # Save file
    # header = ["avg_slope" "err"]
    # writedlm(srtm_maxslope_basin_avg, vcat(header, hcat(avgslope, stdslope)))


## --- Alternatively, load pregenerated OCTOPUS and slope data for each basin
    # Windows machines must start at this step
    # Errors are 1-σ standard deviations
    @info "Loading pre-parsed OCTOPUS and basin slope data"
    octopusdata = importdataset(octopusdata_basins,'\t', importas=:Tuple)
    basin_srtm = importdataset(srtm_maxslope_basin_avg, '\t', importas=:Tuple)


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
    # Uses 1-σ standard deviations
    fobj = yorkfit(c, ex, log.(m), log.(ey))
    emmkyr(slp) = exp(slp * (fobj.slope) + (fobj.intercept))
    export emmkyr
    
    
## --- Plot results
    # h = Plots.plot(xlabel="SRTM15+ Slope (m/km)", ylabel="Erosion rate (mm/kyr)",
    #     framestyle=:box, legend=:topleft, fg_color_legend=:white,
    # )

    # # Data
    # Plots.scatter!(basin_srtm.avg_slope,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", 
    #     msc=:auto, color=:blue, yscale=:log10,
    # )
    # Plots.scatter!(h, basin_srtm.avg_slope,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", 
    #     msc=:auto, color=:orange, markershape=:square
    # )
    # Plots.scatter!(h, c, m, label="Binned Means", yerror=ex, color=:black, lcolor=:black, 
    #     msc=:black,
    # )

    # # Model
    # modelval, = unmeasurementify(emmkyr.(1:600))
    # Plots.plot!(h, 1:length(modelval), modelval, label="Model", color=:black, width=3)
    
    # display(h)
    # savefig("results/figures/erosionslope.pdf")


## --- Model erosion as the product of slope and precipitation
    # using NetCDF

    # # Get basin precipitation
    # precip = ncread("data/prate.sfc.mon.ltm.nc","prate")
    # precip = mean(precip, dims=3)[:,:,1]
    # basin_precip = find_precip(octopusdata.y_wgs84, octopusdata.x_wgs84, precip)

    # # Get slope ⋅ precipitation
    # slopeprecip = (basin_srtm.avg_slope .± basin_srtm.err) .* basin_precip
    # v, e = unmeasurementify(slopeprecip)
    # slopeprecip = (v=v, e=e)

    # # Get erosion
    # t = .!isnan.(slopeprecip.v) .& .!isnan.(slopeprecip.e)
    # t_be = t .& .!isnan.(octopusdata.ebe_mmkyr) .& .!isnan.(octopusdata.ebe_err)
    # t_al = t .& .!isnan.(octopusdata.eal_mmkyr) .& .!isnan.(octopusdata.eal_err)

    # x = (
    #     v = [slopeprecip.v[t_be]; slopeprecip.v[t_al]],
    #     e = [slopeprecip.e[t_be]; slopeprecip.e[t_al]]
    # )
    # y = (
    #     v = [octopusdata.ebe_mmkyr[t_be]; octopusdata.eal_mmkyr[t_al]],
    #     e = [octopusdata.ebe_err[t_be]; octopusdata.eal_err[t_al]]
    # )

    # # Get bin averages in regular space
    # c, m, ex, ey = binmeans_percentile(x.v, y.v, step=5)

    # # Log transform **both** x and y and fit model
    # fobj = yorkfit(log.(c), log.(ex), log.(m), log.(ey))
    # # emmkyr_precip(slp) = exp(slp * (fobj.slope) + (fobj.intercept))

    # # Plot results
    # h = scatter(slopeprecip.v,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", 
    #     msc=:auto, color=:blue, alpha=0.5
    # )
    # scatter!(h, slopeprecip.v,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", 
    #     msc=:auto, color=:orange, alpha=0.5
    # )
    # scatter!(h, c, m, label="Binned Means", msc=:auto, color=:black)

    # # Model
    # modelin = range(start=0, stop=0.4, length=100)
    # modelval, modelerr = unmeasurementify(emmkyr_precip.(modelin))
    # # plot!(h, 1:length(modelval), modelval, label="Model", color=:black, width=3)
    # plot!(xlabel="SRTM15+ Slope (m/km)", ylabel="Erosion rate (mm/kyr)",
    #     xscale=:log10, yscale=:log10, framestyle=:box, legend=:topleft, 
    #     fg_color_legend=:white,
    # )
    # display(h)

## --- End of file