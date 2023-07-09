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
    include("UtilitiesSlope.jl")
    include("Utilities.jl")


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
    # # srtm = h5open("data/srtm15plus_aveslope.h5", "r")
    # srtm = h5open("data/srtm15plus_maxslope.h5", "r")
    # srtm = read(srtm["vars"])
    # srtm = NamedTuple{Tuple(Symbol.(keys(srtm)))}(values(srtm))

    # # Expected runtime 1.5 hrs
    # @timev avgslope, stdslope = get_basin_srtm15plus_aveslope(srtm, nbasins, subbasins, 
    #     basin_polygon_lat, basin_polygon_lon
    # )

    # # Save file
    # header = ["avg_slope" "err"]
    # writedlm("data/basin_srtm15plus_avg_maxslope.tsv", vcat(header, hcat(avgslope, stdslope)))

    # File names:
        # basin_srtm15plus_avg_maxslope -> average of maximum slopes of each point in basin
        # basin_srtm15plus_aveslope     -> average of average slopes of each point in the basin


## --- Alternatively, load pregenerated OCTOPUS and slope data for each basin
    @info "Loading pre-parsed OCTOPUS and basin slope data"
    octopusdata = importdataset("output/octopusdata.tsv",'\t', importas=:Tuple)
    basin_srtm = importdataset("data/basin_srtm15plus_avg_maxslope.tsv", '\t', importas=:Tuple)


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


## --- Fit raw erosion rate as a function of slope (m/km)
    @info "Fitting erosion / slope curve"
        
    # Be data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(basin_srtm.err)              # Basin
    t .&= .!isnan.(octopusdata.ebe_mmkyr) .& .!isnan.(octopusdata.ebe_err)      # Erosion

    xval = basin_srtm.avg_slope[t]
    xerr = basin_srtm.err[t]
    yval = log10.(octopusdata.ebe_mmkyr[t])
    yerr = log10.(octopusdata.ebe_err[t])

    # Al data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(basin_srtm.err)              # Basin
    t .&= .!isnan.(octopusdata.eal_mmkyr) .& .!isnan.(octopusdata.eal_err)      # Erosion

    xval = append!(xval, basin_srtm.avg_slope[t])
    xerr = append!(xerr, basin_srtm.err[t])
    yval = append!(yval, log10.(octopusdata.eal_mmkyr[t]))
    yerr = append!(yerr, log10.(octopusdata.eal_err[t]))

    # Old method with no uncertainty
    # a, b = linreg(xval, yval)
    # emmkyr_old(slp) = 10^(slp * b + a)

    # New method with 1σ uncertainty built into the function
    fobj = yorkfit(xval, xerr, yval, yerr)
    emmkyr(slp) = 10^(slp * (fobj.slope) + (fobj.intercept))


# --- Plot results
    # De-measurement model data
    len = 650
    val, err = unmeasurementify(emmkyr.(1:650))
    model = (val, err)

    # Data
    h = scatter(basin_srtm.avg_slope,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", 
        msc=:auto, color=:blue, alpha=0.5
    )
    scatter!(h, basin_srtm.avg_slope,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", 
        msc=:auto, color=:orange, alpha=0.5
    )

    # Model
    plot!(h, 1:len, emmkyr.(1:len), label="Old model", width=3, color=:red)
    plot!(h, 1:len, model.val, ribbon = model.err, label="New model", width=3, color=:cyan)
    
    plot!(h, xlabel="SRTM15+ Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, 
        fg_color_legend=:white, legend=:topleft, framestyle=:box
    )

    display(h)
    # savefig(h, "slopeerosion_avg.pdf")


## --- Calculate precipitation and runoff for each basin
    @info "Calculating basin precipitation"

    precip = ncread("data/prate.sfc.mon.ltm.nc","prate")
    precip = mean(precip, dims=3)[:,:,1]
    basin_precip = find_precip(octopusdata.y_wgs84, octopusdata.x_wgs84, precip)

    runoff = ncread("data/runof.sfc.mon.ltm.nc", "runof")
    runoff = mean(runoff, dims=3)[:,:,1]
    basin_runoff = find_precip(octopusdata.y_wgs84, octopusdata.x_wgs84, runoff)

    # Test find_precip to make sure it's not flipped
    # This could probably be done better with GeoMakie
    imsc(collect(precip'), viridis)
    lat = repeat(-90:90,1,361)
    lon = repeat((-180:180)',181,1)
    imsc(find_precip(lat,lon, precip),viridis)


## --- Model erosion as a function of slope and precipitation
    # Slope-precipitation (kg⋅s/km⋅m ? = 1000kg⋅s/km ??)
    slopeprecip = (basin_srtm.avg_slope .± basin_srtm.err) .* basin_precip
    slopeprecipval = Array{Float64}(undef, length(slopeprecip), 1)
    slopepreciperr = Array{Float64}(undef, length(slopeprecip), 1)
    for i in eachindex(slopeprecip)
        slopeprecipval[i] = slopeprecip[i].val
        slopepreciperr[i] = slopeprecip[i].err
    end
    
    # Be data
    t = vec(.!isnan.(slopeprecipval) .& .!isnan.(slopepreciperr))               # Basin
    t .&= .!isnan.(octopusdata.ebe_mmkyr) .& .!isnan.(octopusdata.ebe_err)      # Erosion

    xval = slopeprecipval[t]
    xerr = slopepreciperr[t]
    yval = log10.(octopusdata.ebe_mmkyr[t])
    yerr = log10.(octopusdata.ebe_err[t])

    # Al data
    t = vec(.!isnan.(slopeprecipval) .& .!isnan.(slopepreciperr))               # Basin
    t .&= .!isnan.(octopusdata.eal_mmkyr) .& .!isnan.(octopusdata.eal_err)      # Erosion

    xval = append!(xval, slopeprecipval[t])
    xerr = append!(xerr, slopepreciperr[t])
    yval = append!(yval, log10.(octopusdata.eal_mmkyr[t]))
    yerr = append!(yerr, log10.(octopusdata.eal_err[t]))

    fobj = yorkfit(xval, xerr, yval, yerr)
    emmkyr_precip(slp) = 10^(slp * (fobj.slope) + (fobj.intercept))


## --- Plot results
    # De-measurement model data
    rng = 0:0.005:0.04
    model = (val = zeros(length(rng)), err = zeros(length(rng)))
    for i = 1:length(rng)
        e = emmkyr_precip(i/300)
        model.val[i] = e.val
        model.err[i] = e.err
    end

    # Data
    h = scatter(slopeprecipval,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", 
        msc=:auto, color=:blue, alpha=0.5
    )
    scatter!(h, slopeprecipval,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", 
        msc=:auto, color=:orange, alpha=0.5
    )

    # Model
    plot!(h, rng, model.val, ribbon = model.err, label="Model", width=3, color=:cyan)
    
    plot!(h, xlabel="SRTM15+ Slope (m/km) ⋅ Precipitation rate (kg/m²/s)", 
        ylabel="Erosion rate (mm/kyr)", yscale=:log10, fg_color_legend=:white, 
        legend=:topleft, framestyle=:box
    )

    display(h)


## --- Experiment with changepoint
    # Sort data by slope
    perm = sortperm(xval)
    ordered_yval = yval[perm]
    ordered_yerr = yerr[perm]

    # Get index of potential change point
    dist = changepoint(ordered_yval, ordered_yerr, 10000, np=1)
    while true
        dist = dist[9000:end]
        sum(dist) > 0 && return dist

        dist = changepoint(ordered_yval, ordered_yerr, 10000)
    end
    point = unique(dist)

    # Restrict model input to data before the change
    t = (xval .< xval[point])

    # Recalculate model
    xval_change = xval[t]
    xerr_change = xerr[t]
    yval_change = yval[t]
    yerr_change = yerr[t]

    fobj = yorkfit(xval_change, xerr_change, yval_change, yerr_change)
    emmkyr_cp(slp) = 10^(slp * (fobj.slope) + (fobj.intercept))


## --- Plot results
    # De-measurement model data
    rng = 0:0.005:0.04
    model = (val = zeros(length(rng)), err = zeros(length(rng)))
    for i = 1:length(rng)
        e = emmkyr_cp(i/300)
        model.val[i] = e.val
        model.err[i] = e.err
    end

    # Data
    h = scatter(slopeprecipval,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", 
        msc=:auto, color=:blue, alpha=0.5
    )
    scatter!(h, slopeprecipval,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", 
        msc=:auto, color=:orange, alpha=0.5
    )

    # Model
    plot!(h, rng, model.val, ribbon = model.err, label="Changepoint Model", width=3, 
        color=:cyan
    )
    
    plot!(h, xlabel="SRTM15+ Slope (m/km) ⋅ Precipitation rate (kg/m²/s)", 
        ylabel="Erosion rate (mm/kyr)", yscale=:log10, fg_color_legend=:white, 
        legend=:topleft, framestyle=:box
    )

    display(h)


## --- End of file