## --- Set up
    # Packages
    using StatGeochem
    using ProgressMeter
    using Statistics
    # using LsqFit: curve_fit
    using HDF5
    using Measurements
    using DelimitedFiles
    using Plots
    using NetCDF
    using Isoplot: yorkfit

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


## --- Alternatively, load pregenerated slope data for each basin
    @info "Loading basin slope data"
    basin_srtm = importdataset("data/basin_srtm15plus_avg_maxslope.tsv", '\t', importas=:Tuple)

    
## --- Fit raw erosion rate as a function of slope (m/km)
    @info "Fitting erosion / slope curve"
        
    # Be data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(basin_srtm.err)              # Basin
    t .&= .!isnan.(octopusdata.ebe_mmkyr) .& .!isnan.(octopusdata.ebe_err)      # Erosion
    t .&= (basin_srtm.avg_slope .< 300)                                         # Slope [TK: why?]

    xval = basin_srtm.avg_slope[t]
    xerr = zeronan!(basin_srtm.err[t])
    yval = log10.(octopusdata.ebe_mmkyr[t])
    yerr = zeronan!(log10.(octopusdata.ebe_err[t]))

    # Al data
    t = .!isnan.(basin_srtm.avg_slope) .& .!isnan.(basin_srtm.err)              # Basin
    t .&= .!isnan.(octopusdata.eal_mmkyr) .& .!isnan.(octopusdata.eal_err)      # Erosion
    t .&= (basin_srtm.avg_slope .< 300)                                         # Slope [TK: why?]

    xval = append!(xval, basin_srtm.avg_slope[t])
    xerr = append!(xerr, zeronan!(basin_srtm.err[t]))
    yval = append!(yval, log10.(octopusdata.eal_mmkyr[t]))
    yerr = append!(yerr, zeronan!(log10.(octopusdata.eal_err[t])))

    # Fit curve
    # In theory, couldn't we use margin_error and confidence_interval?
    # p = [0.5, 1/100]
    # fobj = curve_fit(linear, xval, yval, p)
    
    # function emmkyr(slp)
    #     return 10^(slp * (fobj.param[2]) + (fobj.param[1]))
    # end

    fobj = yorkfit(xval, xerr, yval, yerr)

    # Uncertainties are built into this function
    function emmkyr(slp)
        return 10^(slp * (fobj.slope) + (fobj.intercept))
    end

    # Parameters using maximum slope:
        # 1.0631462800868299
        # 0.003687823054086418

#         plot(data.Age_Bin_Myr, (data.Lower .+ data.Upper) ./ 2, ribbon = (data.Lower .- data.Upper) ./ 2, 
#     fillalpha = 1, color=:skyblue, label = ""
# )
    
## --- Plot results
    # De-measurement model data
    len = 500
    model = (val = zeros(len), err = zeros(len))
    for i in 1:len
        erosion = emmkyr(i)
        model.val[i] = erosion.val
        model.err[i] = erosion.err
    end

    # Data
    h = scatter(basin_srtm.avg_slope,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", msc=:auto, color=:blue, alpha=0.5)
    scatter!(h, basin_srtm.avg_slope,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", msc=:auto, color=:orange, alpha=0.5)

    # Model
    # plot!(h, 1:500, ribbon = model.err)
    plot!(h, 1:500, model.val, label = "E = ", width=3, color=:black, fg_color_legend=:white)
    plot!(xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, legend=:topleft, framestyle=:box)

    savefig(h, "slopeerosion_avg.pdf")


## --- Calculate precipitation and runoff for each basin
    @info "Calculating basin precipitation"

    precip = ncread("data/prate.sfc.mon.ltm.nc","prate")
    precip = mean(precip, dims=3)[:,:,1]
    basin_precip = find_precip(octopusdata.y_wgs84, octopusdata.x_wgs84, precip)

    runoff = ncread("data/runof.sfc.mon.ltm.nc", "runof")
    runoff = mean(runoff, dims=3)[:,:,1]
    basin_runoff = find_precip(octopusdata.y_wgs84, octopusdata.x_wgs84, runoff)

    # # Test find_precip to make sure it's not flipped
    # imsc(collect(precip'), viridis)
    # lat = repeat(-90:90,1,361)
    # lon = repeat((-180:180)',181,1)
    # imsc(find_precip(lat,lon),viridis)


## --- See if adding in a precipitation term makes this work any better
    # New x term
    newx = @. basin_srtm.avg_slope      # No precipitation term because not really working

    # Be data
    t = .!isnan.(newx) .& .!isnan.(octopusdata.ebe_mmkyr) .& (basin_srtm.avg_slope .< 300)
    x = newx[t]
    y = log10.(octopusdata.ebe_mmkyr[t])

    # Al data
    t = .!isnan.(newx) .& .!isnan.(octopusdata.eal_mmkyr) .& (basin_srtm.avg_slope .< 300)
    x = append!(x, newx[t])
    y = append!(y, log10.(octopusdata.eal_mmkyr[t]))

    # Fit curve
    p = [0.5, 1/100]
    fobj = curve_fit(linear, x, y, p)

    function newemmkyr(slp)
        return 10^(slp * (fobj.param[2]) + (fobj.param[1]))
    end

    # Plot
    modellow = round(nanminimum(newx), sigdigits=2)
    modelhigh = round(nanmaximum(newx), sigdigits=2)
    modrng = range(start=modellow, stop=modelhigh, length=50)

    h = scatter(newx,octopusdata.ebe_mmkyr, label="OCTOPUS Be-10 data", msc=:auto, color=:blue)
    scatter!(h, newx,octopusdata.eal_mmkyr, label="OCTOPUS Al-26 data", msc=:auto, color=:orange)
    plot!(h, modrng, newemmkyr.(modrng), label = "E = 10^(0.056) + 1.00", width=3, color=:black, fg_color_legend=:white)
    plot!(xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, legend=:topleft, framestyle=:box)
    savefig(h, "modeltest.pdf")

## --- End of file