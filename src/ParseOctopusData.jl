    1+1

## --- Load external packages

    using StatGeochem
    using Plots; gr();
    using ProgressMeter: @showprogress

    if VERSION>=v"0.7"
        using Statistics, DelimitedFiles, SpecialFunctions
    end

## --- Load the OCTOPUS kml file
    # Unzip the kml file

    function load_octopus_kml(filename)

        # Check number and ordering of coordinate lists
        fid = open("data/octopus/crn_basins_global.kml")

        # Count the number of lines in the file
        i = 0
        for line in eachline(fid)
            i += 1
        end
        isfirstcoord = Array{Bool}(undef, i)

        # Parse the coords (basin polygon outlines) in each placemark
        i = 0
        ncoords = 0
        lastcoord = 0
        lastplacemark = 0
        seekstart(fid)
        for line in eachline(fid)
            i += 1
            if occursin("</ExtendedData>", line)
                lastplacemark = i
            elseif occursin("<coordinates>", line)
                ncoords += 1
                # If we've had a new placemark more recently than a cordinate, then
                # this is the first coord of the placemark
                if lastplacemark > lastcoord
                    isfirstcoord[ncoords] = true
                else
                    isfirstcoord[ncoords] = false
                end
                lastcoord = i
            end
        end
        isfirstcoord = isfirstcoord[1:ncoords]

        # Read entire string
        seekstart(fid)
        str = read(fid,String)

        # Close the file
        close(fid)

        nbasins = 0
        subbasins = Array{Any}(undef, sum(isfirstcoord))
        for i = 1:length(isfirstcoord)
            if isfirstcoord[i]
                nbasins += 1
                subbasins[nbasins] = i
            else
                subbasins[nbasins] = vcat(subbasins[nbasins],i)
            end
        end

        return (str, isfirstcoord, nbasins, subbasins)
    end

    run(`gunzip data/octopus/crn_basins_global.kml.gz`); # Decompress

    (str, isfirstcoord, nbasins, subbasins) = load_octopus_kml("data/octopus/crn_basins_global.kml")

    run(`gzip data/octopus/crn_basins_global.kml`); # Recompress

## --- Parse the file

    function parse_octopus_kml_variables(str)
        # Find a list of all the string variables
        lm = eachmatch(r"<SimpleField type=\"string\" name=\"(.*?)\"", str);
        stringvars = unique(map(x -> x[1], lm))

        # Find a list of all the numeric variables
        lm = eachmatch(r"<SimpleField type=\"double\" name=\"(.*?)\"", str);
        numvars = unique(map(x -> x[1], lm))

        # Make a dict to store all our parsed data
        data = Dict();

        # Parse all the string variables
        for i = 1:length(stringvars)
            lm = eachmatch(Regex("<SimpleData name=\"$(stringvars[i])\">(.*?)</SimpleData>"),str);
            data[stringvars[i]] = map(x -> x[1], lm);
        end

        # Parse all the numeric variables
        for i = 1:length(numvars)
            lm = eachmatch(Regex("<SimpleData name=\"$(numvars[i])\">(.*?)</SimpleData>"),str);
            data[numvars[i]] = map(x -> parse(Float64, x[1]), lm)
        end

        needsNaNs = ["alcorr", "be10ep", "beprod", "eal_err", "elev_std", "errbe_prod", "al26ep", "alprod", "be10ep_err", "beself", "eal_gcmyr", "erral_ams", "errbe_tot", "al26ep_err", "alself", "be10nc", "besnow", "eal_mmkyr", "erral_muon", "sizemax", "al26nc", "alsnow", "be10nc_err", "betopo", "ebe_err", "erral_prod", "sizemin", "al26nc_err", "altopo", "be10np", "betots", "ebe_gcmyr", "erral_tot", "slp_ave", "al26np", "altots", "be10np_err", "ebe_mmkyr", "errbe_ams", "slp_std", "al26np_err", "area", "becorr", "dbver", "elev_ave", "errbe_muon", ]

        for i = 1:length(needsNaNs)
            data[needsNaNs[i]][data[needsNaNs[i]] .< 0] .= NaN
        end

        return data
    end

    data = parse_octopus_kml_variables(str)

    # Parse the basin polygon outlines
    function parse_octopus_polygon_outlines(str,isfirstcoord)
        i = 0
        n = 0
        basin_polygon_lat = Array{Array}(undef, length(isfirstcoord))
        basin_polygon_lon = Array{Array}(undef, length(isfirstcoord))
        basin_polygon_n = Array{Int}(undef, length(isfirstcoord))
        for m in eachmatch(r"<coordinates>\n\t*(.*?)\n",str)
            i += 1 # Number of parsed values
            parsed = delim_string_function(x -> delim_string_parse(x, ',', Float64), m[1], ' ', Array{Float64,1})
            nparsed = length(parsed)

            basin_polygon_n[i] = nparsed;
            basin_polygon_lon[i] = parsed .|> x -> x[1]
            basin_polygon_lat[i] = parsed .|> x -> x[2]
        end
        return (basin_polygon_n, basin_polygon_lat, basin_polygon_lon)
    end

    (basin_polygon_n, basin_polygon_lat, basin_polygon_lon) = parse_octopus_polygon_outlines(str,isfirstcoord)

## --- Calculate precipitation for each basin

    using NetCDF
    prate = ncread("prate.sfc.mon.ltm.nc","prate")
    lat = ncread("prate.sfc.mon.ltm.nc","lat")
    lon = ncread("prate.sfc.mon.ltm.nc","lon")

    precip = mean(prate, dims=3)[:,:,1]
    imsc(collect(precip'), viridis)

    function find_precip(lat,lon)
        prate = ncread("prate.sfc.mon.ltm.nc","prate")
        precip = mean(prate, dims=3)[:,:,1]

        out = Array{Float64}(undef,size(lat))
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i] > 88.542 || lat[i] < -88.542
                out[i] = NaN
            else
                row = trunc(Int, mod(lon[i]*192/360 + 0.5, 192.0)) + 1
                col = trunc(Int, 46.99999 - lat[i]*94/177.084) + 1
                out[i] = precip[row, col]
            end
        end

        return out
    end

    data["precip"] = find_precip(data["y_wgs84"],data["x_wgs84"])

    # # Test find_precip to make sure it's not flipped
    # lat = repmat(-90:90,1,361)
    # lon = repmat((-180:180)',181,1)
    # imsc(find_precip(lat,lon),viridis)

## --- Recalculate slopes for each basin

    function get_basin_srtm15plus_aveslope(srtm::Dict,nbasins,subbasins,basin_polygon_lat,basin_polygon_lon)
        slope = srtm["slope"]
        x_lon_cntr = srtm["x_lon_cntr"]
        y_lat_cntr = srtm["y_lat_cntr"]

        basin_srtm15plus_aveslope = Array{Float64}(undef, nbasins)
        basin_N = Array{Int64}(undef, nbasins)
        for i = 1:nbasins
            rowsinbasin = Array{Int}(undef, 0)
            columnsinbasin = Array{Int}(undef, 0)
            for j = 1:length(subbasins[i])
                k = subbasins[i][j]
                subbasin_lon = basin_polygon_lon[k]
                subbasin_lat = basin_polygon_lat[k]
                (row, column) = find_grid_inpolygon(x_lon_cntr, y_lat_cntr, subbasin_lon, subbasin_lat)
                rowsinbasin = vcat(rowsinbasin,row)
                columnsinbasin = vcat(columnsinbasin,column)
            end

            pointsinbasin = unique(hcat(rowsinbasin,columnsinbasin),dims=1);
            basin_slopes = Array{UInt16}(undef, length(rowsinbasin))
            for j=1:size(pointsinbasin,1)
                basin_slopes[j] = slope[pointsinbasin[j,1],pointsinbasin[j,2]]
            end

            # Average all the slopes
            basin_srtm15plus_aveslope[i] = mean(basin_slopes)
            basin_N[i] = length(basin_slopes)

            print("Basin: $i, slope: $(round(basin_srtm15plus_aveslope[i],digits=3)), N: $(basin_N[i]) \n")
        end

        return basin_srtm15plus_aveslope
    end

    # using HDF5
    # srtm = h5read("data/srtm15plus_aveslope.h5","vars/")
    # basin_srtm15plus_aveslope = get_basin_srtm15plus_aveslope(srtm,nbasins,subbasins,basin_polygon_lat,basin_polygon_lon)


    using JLD
    # save("basin_srtm15plus_aveslope.jld","basin_srtm15plus_aveslope",basin_srtm15plus_aveslope)
    basin_srtm15plus_aveslope =  load("basin_srtm15plus_aveslope.jld")["basin_srtm15plus_aveslope"]

    # save("OctopusSlopeRecalc.jld","basin_srtm15plus_aveslope",basin_srtm15plus_aveslope,"basin_polygon_n",basin_polygon_n,"basin_polygon_lat",basin_polygon_lat,"basin_polygon_lon",basin_polygon_lon,"subbasins",subbasins)


## --- Plot raw Lat and Lon
    h = plot(data["x_wgs84"],data["y_wgs84"],seriestype=:scatter,label="basin locations")
    plot!(h, xlabel="Longitude", ylabel="Latitude")
    savefig(h, "Basin_Locations.pdf")
    display(h)

## --- Plot raw erosion rate vs basin area
    h = plot(data["area"],data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
    plot!(h, data["area"],data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
    plot!(h, xlabel="Basin area (km^2)", ylabel="Erosion rate (mm/kyr)",yscale=:log10, xscale=:log10, legend=:topleft, fg_color_legend=:white)
    savefig(h, "Area_vs_erosion.pdf")
    display(h)

## --- Plot raw erosion rate vs latitude
    h = plot(data["y_wgs84"],data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
    plot!(h,data["y_wgs84"],data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
    plot!(h, xlabel="Latitude", ylabel="Erosion rate (mm/kyr)", yscale=:log10,legend=:topleft, fg_color_legend=:white)
    savefig(h, "Latitude_vs_erosion.pdf")
    display(h)

## --- Plot raw erosion rate vs precipitation
    h = plot(data["precip"]*31557600,data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
    plot!(h,data["precip"]*31557600,data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
    plot!(h, xlabel="Precipitation (kg/m^2/yr)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, fg_color_legend=:white)
    savefig(h, "Precipitation_vs_erosion.pdf")
    display(h)

## --- Plot raw erosion rate vs elevation
    h = plot(data["elev_ave"],data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
    plot!(h,data["elev_ave"],data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
    plot!(h, xlabel="Elevation (m)", ylabel="Erosion rate (mm/kyr)", yscale=:log10,legend=:topleft, fg_color_legend=:white)
    savefig(h, "Elevation_vs_erosion.pdf")
    display(h)

## --- Plot new slope vs erosion rate
    h = plot(basin_srtm15plus_aveslope,data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS Be-10 data")
    plot!(h,basin_srtm15plus_aveslope,data["eal_mmkyr"],seriestype=:scatter,label="OCTOPUS Al-26 data")
    plot!(h, xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10,legend=:topleft, fg_color_legend=:white)
    savefig(h, "Slope_vs_Erosion.pdf")
    display(h)

    # log10(E) = slope/140 + 0.7
    # E = 10^(slope/140 + 0.7)

## --- Plot new slope vs precipitation
    h = plot(basin_srtm15plus_aveslope,data["precip"]*31557600,seriestype=:scatter,label="")
    plot!(h, xlabel="SRTM15 Slope (m/km)", ylabel="Precipitation (kg/m^2/yr)", yscale=:log10,legend=:topleft, fg_color_legend=:white)
    savefig(h, "Slope_vs_Precipitation.pdf")
    display(h)

## --- Fit raw erosion rate as a function of slope (m/km)

    using LsqFit: curve_fit

    function linear(x,p)
        y = p[1] .+ x * p[2]
    end
    p = [0.5, 1/100]
    t = .~isnan.(basin_srtm15plus_aveslope) .& .~isnan.(data["ebe_mmkyr"]) .& (basin_srtm15plus_aveslope .< 300)
    x = basin_srtm15plus_aveslope[t]
    y = log10.(data["ebe_mmkyr"][t])
    t = .~isnan.(basin_srtm15plus_aveslope) .& .~isnan.(data["eal_mmkyr"]) .& (basin_srtm15plus_aveslope .< 300)
    x = append!(x, basin_srtm15plus_aveslope[t])
    y = append!(y, log10.(data["eal_mmkyr"][t]))
    fobj = curve_fit(linear, x, y, p);

    # fobj.param[1]: 0.987237
    # fobj.param[2]: 0.00555159
    function Emmkyr(slp)
        return 10^(slp*0.00567517 + 0.971075)
    end

    h = plot(basin_srtm15plus_aveslope,data["ebe_mmkyr"], seriestype=:scatter, label="OCTOPUS Be-10 data")
    plot!(h, basin_srtm15plus_aveslope,data["eal_mmkyr"], seriestype=:scatter, label="OCTOPUS Al-26 data")
    plot!(h, xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)", yscale=:log10, legend=:topleft)
    plot!(h, 1:500, Emmkyr.(1:500), label = "E = 10^(slp/176.2 + 0.971)", fg_color_legend=:white)
    savefig(h, "Slope_vs_Erosion_Fitted.pdf")
    display(h)

## --- Plot Be-10 vs Al-26 erosion rate
    h = plot(data["ebe_mmkyr"],data["eal_mmkyr"],seriestype=:scatter,label="")
    plot!(h, xlabel="Be-10 erosion rate (t1/2=1.387E6) (mm/kyr) ", ylabel="Al-26 erosion rate (t1/2=7.17E5) (mm/kyr)",xscale=:log,yscale=:log)
    plot!(h,10.0.^(-1:4),10.0.^(-1:4),label="1:1",legend=:topleft,fg_color_legend=:white)

    function linslp1(x,p)
        y = p[1] .+ x
    end
    p = [0.0]
    t = .~isnan.(data["ebe_mmkyr"]) .& .~isnan.(data["eal_mmkyr"])
    x = log10.(data["ebe_mmkyr"][t])
    y = log10.(data["eal_mmkyr"][t])
    fobj = curve_fit(linslp1, x, y, p);

    plot!(h,10.0.^(-1:4), 10.0.^((-1:4) .+ fobj.param[1]), label="1:$(10.0^fobj.param[1])")
    savefig(h, "Be_vs_Al_rate.pdf")
    display(h)

## --- Plot residual erosion rate as a function of latitude
    resid = data["ebe_mmkyr"] - Emmkyr.(basin_srtm15plus_aveslope)
    h = plot(abs.(data["y_wgs84"][resid.<0]), abs.(resid[resid.<0]), seriestype=:scatter, label="negative residuals", color=:red)
    plot!(h, abs.(data["y_wgs84"][resid.>0]), resid[resid.>0], seriestype=:scatter, label="positive residuals", color=:blue)
    plot!(h, xlabel="Latitude", ylabel="Residual erosion rate (mm/kyr)", yscale=:log10, legend=:topleft, fg_color_legend=:white)
    savefig(h, "Latitude_vs_erosion_residual.pdf")
    display(h)


## --- Plot residual erosion rate as a function of precipitation
    resid = data["ebe_mmkyr"] - Emmkyr.(basin_srtm15plus_aveslope)
    h = plot(data["precip"][resid.<0]*31557600, abs.(resid[resid.<0]), seriestype=:scatter, label="negative residuals", marker=(0.3,:red))
    plot!(h, data["precip"][resid.>0]*31557600, resid[resid.>0], seriestype=:scatter, label="positive residuals", marker=(0.3,:blue))
    plot!(h, xlabel="Precipitation (kg/m^2/yr)", ylabel="Residual erosion rate (mm/kyr)", yscale=:log10, fg_color_legend=:white)
    savefig(h, "Precipitation_vs_erosion_residual.pdf")
    display(h)


## --- Plot residual erosion rate as a function of elevation
    resid = data["ebe_mmkyr"] - Emmkyr.(basin_srtm15plus_aveslope)
    h = plot(abs.(data["elev_ave"][resid.<0]), abs.(resid[resid.<0]), seriestype=:scatter, label="negative residuals",  marker=(0.3,:red))
    plot!(h, abs.(data["elev_ave"][resid.>0]), resid[resid.>0], seriestype=:scatter, label="positive residuals", marker=(0.3,:blue))
    # h = plot(data["elev_ave"],abs.(resid),seriestype=:scatter,label="OCTOPUS 10Be data",legend=:topleft)
    plot!(h, xlabel="Elevation", ylabel="Residual erosion rate (mm/kyr)",yscale=:log10,fg_color_legend=:white)
    savefig(h, "Elevation_vs_erosion_residual.pdf")
    display(h)

## --- Plot residual erosion ratio as a function of elevation
    resid = data["ebe_mmkyr"] ./ Emmkyr.(basin_srtm15plus_aveslope)
    h = plot(data["elev_ave"],abs.(resid),seriestype=:scatter,label="OCTOPUS 10Be data",legend=:topleft)
    plot!(h, xlabel="Elevation", ylabel="Residual erosion ratio",yscale=:log10)
    # savefig(h, "Latitude_vs_erosion_residual.pdf")
    display(h)


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
