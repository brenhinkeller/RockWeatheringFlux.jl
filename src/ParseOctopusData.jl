## --- Load external packages
    using StatGeochem           # Used for slope
    using Plots; gr();          # Commented out
    using ProgressMeter         # Used for slope
    using Statistics            # Used for slope
    using DelimitedFiles
    using SpecialFunctions
    using JLD                   # Used for slope
    using NetCDF
    using LsqFit: curve_fit     # Used for slope
    using HDF5                  # Used for slope

    # Local utilities
    include("OCTOPUS_utilities.jl")     # TO DO: re-organize utilities files and give these better names


## --- Load and parse the OCTOPUS kml file
    @info "Loading OCTOPUS data"

    # Decompress, load, and recompress 
    run(`gunzip data/octopus/crn_basins_global.kml.gz`);
    (str, isfirstcoord, nbasins, subbasins) = load_octopus_kml("data/octopus/crn_basins_global.kml");
    run(`gzip data/octopus/crn_basins_global.kml`);                                           

    # Parse and export variables and data
    data = parse_octopus_kml_variables(str)
    exportdataset(data,"output/octopusdata.tsv",'\t')

    # Get basin polygon outlines
    (basin_polygon_n, basin_polygon_lat, basin_polygon_lon) = parse_octopus_polygon_outlines(str,isfirstcoord)


## --- Calculate slope for each basin
    @info "Calculating slope for each basin"

    srtm = h5read("data/srtm15plus_aveslope.h5","vars/")
    basin_srtm15plus_aveslope = get_basin_srtm15plus_aveslope(srtm, nbasins, subbasins, 
        basin_polygon_lat, basin_polygon_lon
    )

    # Slope can be loaded from here instead of calculating everything again
    # Renamed from OctopusSlopeRecalc.jld
    save("data/OCTOPUS_basin_aveslope_new.jld","basin_srtm15plus_aveslope",basin_srtm15plus_aveslope,
        "basin_polygon_n",basin_polygon_n,"basin_polygon_lat",basin_polygon_lat,
        "basin_polygon_lon",basin_polygon_lon, "subbasins",subbasins
    )

    # Alternatively, loading from this file is about 20x faster
    # save("data/basin_srtm15plus_aveslope.jld","basin_srtm15plus_aveslope",basin_srtm15plus_aveslope)


## --- Alternatively, load pregenerated slope data for each basin
    # @info "Loading basin slope data"
    # basin_srtm15plus_aveslope_2 =  load("data/OCTOPUS_basin_aveslope.jld")["basin_srtm15plus_aveslope"]

## --- Fit raw erosion rate as a function of slope (m/km)
    @info "Fitting erosion / slope curve"
        
    t = .!isnan.(basin_srtm15plus_aveslope) .& .!isnan.(data["ebe_mmkyr"]) .& (basin_srtm15plus_aveslope .< 300)
    x = basin_srtm15plus_aveslope[t]
    y = log10.(data["ebe_mmkyr"][t])

    t = .!isnan.(basin_srtm15plus_aveslope) .& .!isnan.(data["eal_mmkyr"]) .& (basin_srtm15plus_aveslope .< 300)
    x = append!(x, basin_srtm15plus_aveslope[t])
    y = append!(y, log10.(data["eal_mmkyr"][t]))

    p = [0.5, 1/100]
    fobj = curve_fit(linear, x, y, p)

    mse = mean(fobj.resid .^ 2)         # Mean-square error         0.30764514536299076
    ssr = sum((fobj.resid) .^ 2)        # Sum squared regression    984.4644651615704
    ybar = nanmean(y)                   # Mean                      1.7771592362257302
    sst = sum((y .- ybar) .^ 2)         # Total sum of squares      1626.2293532703075
    r2 = 1 - (ssr/sst)                  # r^2 value                 0.39463368854962777        

    # My parameters are different?
    # Not adding this to the utilities document because I need to figure out a way to get this to
    # talk to the flux code and if I want things to recalculate each time... TO DO

    # fobj.param[1]: 0.987237
    # fobj.param[2]: 0.00555159
    function emmkyr(slp)
        return 10^(slp*fobj.param[2] + fobj.param[1])
    end

    # In the final code, these should be put in as the variables for reproducability
    # I don't wanna run this every time I work on it though, so make a hardcoded version!

    # TO DO: error propagation through this calculation

    # function emmkyr(slp)
    #     return 10^(slp*(0.00567517 ± 0.001) + (0.971075 ± 0.1))
    # end






################################# PLOTS AND UNUSED CODE #################################

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


# Plots don't play nicely with an SSH connection.
# ## --- Plot raw Lat and Lon
#     h = plot(data["x_wgs84"],data["y_wgs84"],seriestype=:scatter,label="basin locations")
#     plot!(h, xlabel="Longitude", ylabel="Latitude")
#     savefig(h, "Basin_Locations.pdf")
#     display(h)

# ## --- Plot raw erosion rate vs basin area
#     h = plot(data["area"],data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
#     plot!(h, data["area"],data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
#     plot!(h, xlabel="Basin area (km^2)", ylabel="Erosion rate (mm/kyr)",yscale=:log10, xscale=:log10, legend=:topleft, fg_color_legend=:white)
#     savefig(h, "Area_vs_erosion.pdf")
#     display(h)

# ## --- Plot raw erosion rate vs latitude
#     h = plot(data["y_wgs84"],data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
#     plot!(h,data["y_wgs84"],data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
#     plot!(h, xlabel="Latitude", ylabel="Erosion rate (mm/kyr)", yscale=:log10,legend=:topleft, fg_color_legend=:white)
#     savefig(h, "Latitude_vs_erosion.pdf")
#     display(h)

# ## --- Plot raw erosion rate vs precipitation
#     h = plot(data["precip"]*31557600,data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
#     plot!(h,data["precip"]*31557600,data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
#     plot!(h, xlabel="Precipitation (kg/m^2/yr)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, fg_color_legend=:white)
#     savefig(h, "Precipitation_vs_erosion.pdf")
#     display(h)

# ## --- Plot raw erosion rate vs elevation
#     h = plot(data["elev_ave"],data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
#     plot!(h,data["elev_ave"],data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
#     plot!(h, xlabel="Elevation (m)", ylabel="Erosion rate (mm/kyr)", yscale=:log10,legend=:topleft, fg_color_legend=:white)
#     savefig(h, "Elevation_vs_erosion.pdf")
#     display(h)

# ## --- Plot new slope vs erosion rate
#     h = plot(basin_srtm15plus_aveslope,data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data", msc=:auto, alpha=:0.75);
#     plot!(h,basin_srtm15plus_aveslope,data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data", msc=:auto, alpha=:0.75);
#     plot!(h, xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10,legend=:topleft, fg_color_legend=:white,
#         framestyle=:box
#     );
#     savefig(h, "Slope_vs_Erosion.pdf");
#     display(h)

#     # log10(E) = slope/140 + 0.7
#     # E = 10^(slope/140 + 0.7)

# ## --- Plot new slope vs precipitation
#     h = plot(basin_srtm15plus_aveslope,data["precip"]*31557600,seriestype=:scatter,label="")
#     plot!(h, xlabel="SRTM15 Slope (m/km)", ylabel="Precipitation (kg/m^2/yr)", yscale=:log10,legend=:topleft, fg_color_legend=:white)
#     savefig(h, "Slope_vs_Precipitation.pdf")
#     display(h)

# ## --- Fit raw erosion rate as a function of slope (m/km)
#     @info "Fitting erosion / slope curve"
    
#     p = [0.5, 1/100]
#     t = .!isnan.(basin_srtm15plus_aveslope) .& .!isnan.(data["ebe_mmkyr"]) .& (basin_srtm15plus_aveslope .< 300)
#     x = basin_srtm15plus_aveslope[t]
#     y = log10.(data["ebe_mmkyr"][t])
#     t = .!isnan.(basin_srtm15plus_aveslope) .& .!isnan.(data["eal_mmkyr"]) .& (basin_srtm15plus_aveslope .< 300)
#     x = append!(x, basin_srtm15plus_aveslope[t])
#     y = append!(y, log10.(data["eal_mmkyr"][t]))
#     fobj = curve_fit(linear, x, y, p)

#     mse = mean(fobj.resid .^ 2)         # Mean-square error         0.30764514536299076
#     ssr = sum((fobj.resid) .^ 2)        # Sum squared regression    984.4644651615704
#     ybar = nanmean(y)                   # Mean                      1.7771592362257302
#     sst = sum((y .- ybar) .^ 2)         # Total sum of squares      1626.2293532703075
#     r2 = 1 - (ssr/sst)                  # r^2 value                 0.39463368854962777        

#     # My parameters are different?
#     # Not adding this to the utilities document because I need to figure out a way to get this to
#     # talk to the flux code and if I want things to recalculate each time... TO DO

#     # fobj.param[1]: 0.987237
#     # fobj.param[2]: 0.00555159
#     function emmkyr(slp)
#         return 10^(slp*fobj.param[2] + fobj.param[1])
#     end

#     # TO DO: error propagation through this calculation

#     # function emmkyr(slp)
#     #     return 10^(slp*(0.00567517 ± 0.001) + (0.971075 ± 0.1))
#     # end

    # h = Plots.plot(basin_srtm15plus_aveslope,data["ebe_mmkyr"], seriestype=:scatter, label="OCTOPUS Be-10 data", msc=:auto, alpha=:0.75);
    # Plots.plot!(h, basin_srtm15plus_aveslope,data["eal_mmkyr"], seriestype=:scatter, label="OCTOPUS Al-26 data", msc=:auto, alpha=:0.75);
    # Plots.plot!(h, xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, legend=:topleft, framestyle=:box);
    # #Plots.plot!(h, 1:500, emmkyr.(1:500), label = "E = 10^(slp*0.00556 + 0.987)", fg_color_legend=:white);
    # Plots.plot!(h, 1:500, emmkyr.(1:500), label = "", fg_color_legend=:white);
    # savefig(h, "Slope_vs_Erosion_Fitted.pdf")
    # display(h)

# ## --- Plot Be-10 vs Al-26 erosion rate
#     h = plot(data["ebe_mmkyr"],data["eal_mmkyr"],seriestype=:scatter,label="")
#     plot!(h, xlabel="Be-10 erosion rate (t1/2=1.387E6) (mm/kyr) ", ylabel="Al-26 erosion rate (t1/2=7.17E5) (mm/kyr)",xscale=:log,yscale=:log)
#     plot!(h,10.0.^(-1:4),10.0.^(-1:4),label="1:1",legend=:topleft,fg_color_legend=:white)

#     function linslp1(x,p)
#         y = p[1] .+ x
#     end
#     p = [0.0]
#     t = .~isnan.(data["ebe_mmkyr"]) .& .~isnan.(data["eal_mmkyr"])
#     x = log10.(data["ebe_mmkyr"][t])
#     y = log10.(data["eal_mmkyr"][t])
#     fobj = curve_fit(linslp1, x, y, p);

#     plot!(h,10.0.^(-1:4), 10.0.^((-1:4) .+ fobj.param[1]), label="1:$(10.0^fobj.param[1])")
#     savefig(h, "Be_vs_Al_rate.pdf")
#     display(h)

# ## --- Plot residual erosion rate as a function of latitude
#     resid = data["ebe_mmkyr"] - Emmkyr.(basin_srtm15plus_aveslope)
#     h = plot(abs.(data["y_wgs84"][resid.<0]), abs.(resid[resid.<0]), seriestype=:scatter, label="negative residuals", color=:red)
#     plot!(h, abs.(data["y_wgs84"][resid.>0]), resid[resid.>0], seriestype=:scatter, label="positive residuals", color=:blue)
#     plot!(h, xlabel="Latitude", ylabel="Residual erosion rate (mm/kyr)", yscale=:log10, legend=:topleft, fg_color_legend=:white)
#     savefig(h, "Latitude_vs_erosion_residual.pdf")
#     display(h)


# ## --- Plot residual erosion rate as a function of precipitation
#     resid = data["ebe_mmkyr"] - Emmkyr.(basin_srtm15plus_aveslope)
#     h = plot(data["precip"][resid.<0]*31557600, abs.(resid[resid.<0]), seriestype=:scatter, label="negative residuals", marker=(0.3,:red))
#     plot!(h, data["precip"][resid.>0]*31557600, resid[resid.>0], seriestype=:scatter, label="positive residuals", marker=(0.3,:blue))
#     plot!(h, xlabel="Precipitation (kg/m^2/yr)", ylabel="Residual erosion rate (mm/kyr)", yscale=:log10, fg_color_legend=:white)
#     savefig(h, "Precipitation_vs_erosion_residual.pdf")
#     display(h)


# ## --- Plot residual erosion rate as a function of elevation
#     resid = data["ebe_mmkyr"] - Emmkyr.(basin_srtm15plus_aveslope)
#     h = plot(abs.(data["elev_ave"][resid.<0]), abs.(resid[resid.<0]), seriestype=:scatter, label="negative residuals",  marker=(0.3,:red))
#     plot!(h, abs.(data["elev_ave"][resid.>0]), resid[resid.>0], seriestype=:scatter, label="positive residuals", marker=(0.3,:blue))
#     # h = plot(data["elev_ave"],abs.(resid),seriestype=:scatter,label="OCTOPUS 10Be data",legend=:topleft)
#     plot!(h, xlabel="Elevation", ylabel="Residual erosion rate (mm/kyr)",yscale=:log10,fg_color_legend=:white)
#     savefig(h, "Elevation_vs_erosion_residual.pdf")
#     display(h)

# ## --- Plot residual erosion ratio as a function of elevation
#     resid = data["ebe_mmkyr"] ./ Emmkyr.(basin_srtm15plus_aveslope)
#     h = plot(data["elev_ave"],abs.(resid),seriestype=:scatter,label="OCTOPUS 10Be data",legend=:topleft)
#     plot!(h, xlabel="Elevation", ylabel="Residual erosion ratio",yscale=:log10)
#     # savefig(h, "Latitude_vs_erosion_residual.pdf")
#     display(h)


## ---
#
#     # nbins = 20;
#     # nsims = 10000;
#     # (c, m, el, eu) = bin_bsr_means(data["slp_ave"], log10.(data["ebe_mmkyr"]), 0, 800, nbins, data["slp_std"]/10, nsims)
#     # plot(c, m, yerror=(el,eu), seriestype=:scatter, label="",xlabel="Slope", ylabel="Log10 Erosion rate (mm/kyr)")
#
#     (c,m,sigma) = binmeans(data["slp_ave"], log10.(data["ebe_mmkyr"]), 0, 800, 20)
#     plot(c,m,yerror=2*sigma,seriestype="scatter")
#     plot!(xlabel="Slope", ylabel="Log10 Erosion rate (mm/kyr)")
#
# ## ---
#     (c,m,sigma) = binmeans(data["elev_ave"], data["ebe_mmkyr"], 0, 6000, 12)
#     plot(c,m,yerror=2*sigma,seriestype="scatter")
#     plot!(xlabel="Elevation (m)", ylabel="Erosion rate (mm/kyr)")
#
# ## --- Plot raw erosion rate vs elevation scaled to snowline
#     plot(data["elev_ave"]./data["Snowline"],data["ebe_mmkyr"],seriestype=:scatter)
#     plot!(xlabel="Elevation/ELA", ylabel="Erosion rate (mm/kyr)",yscale=:log10)
#
# ## ---
#     (c,m,sigma) = binmeans(data["elev_ave"]./data["Snowline"], data["ebe_mmkyr"], 0, 1.2, 12)
#     plot(c,m,yerror=2*sigma,seriestype="scatter")
#     plot!(xlabel="Elevation/ELA", ylabel="Erosion rate (mm/kyr)")
#
# ## --- Plot resampled elevation vs erosion rate
#     nbins = 20;
#     nsims = 10000;
#     (c, m, el, eu) = bin_bsr_means(data["elev_ave"], data["ebe_mmkyr"], 0, 6000, nbins, data["elev_std"], nsims)
#     plot(c, m, yerror=(el,eu), seriestype=:scatter, label="",xlabel="Elevation (m)", ylabel="Erosion rate (mm/kyr)")
#
# ## --- Plot resampled elevation / snowline vs erosion rate
#     nbins = 20;
#     nsims = 10000;
#     (c, m, el, eu) = bin_bsr_means(data["elev_ave"]./data["Snowline"], data["ebe_mmkyr"], 0, 1.2, nbins, data["elev_std"]./data["Snowline"], nsims)
#     plot(c, m, yerror=(el,eu), seriestype=:scatter, label="")
#     plot!(xlabel="Elevation/ELA", ylabel="Erosion rate (mm/kyr)")

## --- End of File
