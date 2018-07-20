## ---
    1+1
## --- Load external packages

    using StatGeochem
    using Plots
    gr();
    using ProgressMeter: @showprogress

## --- Load the OCTOPUS kml file

    # Unzip the kml file
    run(`gunzip data/octopus/crn_basins_global.kml.gz`) # Decompress

    # Check number and ordering of coordinate lists
    fid = open("data/octopus/crn_basins_global.kml")

    # Count the number of lines in the file
    i = 0
    for line in eachline(fid)
        i += 1
    end
    isfirstcoord = Array{Bool}(i)

    # Parse the coords (basin polygon outlines) in each placemark
    i = 0
    ncoords = 0
    lastcoord = 0
    lastplacemark = 0
    seekstart(fid)
    for line in eachline(fid)
        i += 1
        if ismatch(r"</ExtendedData>", line)
            lastplacemark = i
        elseif ismatch(r"<coordinates>", line)
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
    str = read(fid,String);

    # Close the file
    close(fid);

    # Recompress
    run(`gzip data/octopus/crn_basins_global.kml`) # Recompress

    nbasins = 0
    subbasins = Array{Any}(sum(isfirstcoord))
    for i = 1:length(isfirstcoord)
        if isfirstcoord[i]
            nbasins += 1
            subbasins[nbasins] = i
        else
            subbasins[nbasins] = vcat(subbasins[nbasins],i)
        end
    end
## --- Parse the file

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
        lm = eachmatch(Regex("<SimpleData name=\"$(stringvars[i])\">(.*?)</SimpleData>"),str)
        data[stringvars[i]] = map(x -> x[1], lm);
    end

    # Parse all the numeric variables
    for i = 1:length(numvars)
        lm = eachmatch(Regex("<SimpleData name=\"$(numvars[i])\">(.*?)</SimpleData>"),str)
        data[numvars[i]] = map(x -> parse(x[1]), lm);
    end

    # Parse the basin polygon outlines
    i = 0
    n = 0
    basin_polygon_lat = Array{Array}(length(isfirstcoord))
    basin_polygon_lon = Array{Array}(length(isfirstcoord))
    basin_polygon_n = Array{Int}(length(isfirstcoord))
    for m in eachmatch(r"<coordinates>\n\t*(.*?)\n",str)
        i += 1 # Number of parsed values
        parsed = delim_string_function(x -> delim_string_parse(x, ',', Float64), m[1], ' ', Array{Float64,1})
        nparsed = length(parsed)

        basin_polygon_n[i] = nparsed;
        basin_polygon_lon[i] = parsed .|> x -> x[1]
        basin_polygon_lat[i] = parsed .|> x -> x[2]
    end

    needsNaNs = ["alcorr", "be10ep", "beprod", "eal_err", "elev_std", "errbe_prod", "al26ep", "alprod", "be10ep_err", "beself", "eal_gcmyr", "erral_ams", "errbe_tot", "al26ep_err", "alself", "be10nc", "besnow", "eal_mmkyr", "erral_muon", "sizemax", "al26nc", "alsnow", "be10nc_err", "betopo", "ebe_err", "erral_prod", "sizemin", "al26nc_err", "altopo", "be10np", "betots", "ebe_gcmyr", "erral_tot", "slp_ave", "al26np", "altots", "be10np_err", "ebe_mmkyr", "errbe_ams", "slp_std", "al26np_err", "area", "becorr", "dbver", "elev_ave", "errbe_muon", ]

    for i = 1:length(needsNaNs)
        data[needsNaNs[i]][data[needsNaNs[i]] .< 0] = NaN;
    end


## --- Recalculate slopes

    using HDF5
    srtm = h5read("data/srtm15plus_aveslope.h5","vars/")
    slope = srtm["slope"]
    x_lon_cntr = srtm["x_lon_cntr"]
    y_lat_cntr = srtm["y_lat_cntr"]

    basin_srtm15plus_aveslope = Array{Float64}(nbasins)
    for i = 1:nbasins
        rowsinbasin = Array{Int}(0)
        columnsinbasin = Array{Int}(0)
        for j = 1:length(subbasins[i])
            k = subbasins[i][j]
            subbasin_lon = basin_polygon_lon[k]
            subbasin_lat = basin_polygon_lat[k]
            (row, column) = find_grid_inpolygon(x_lon_cntr, y_lat_cntr, subbasin_lon, subbasin_lat)
            rowsinbasin = vcat(rowsinbasin,row)
            columnsinbasin = vcat(columnsinbasin,column)
        end

        pointsinbasin = unique(hcat(rowsinbasin,columnsinbasin),1);
        basin_slopes = Array{UInt16}(length(rowsinbasin))
        for j=1:size(pointsinbasin,1)
            basin_slopes[j] = slope[pointsinbasin[j,1],pointsinbasin[j,2]]
        end

        # Average all the slopes
        basin_srtm15plus_aveslope[i] = mean(basin_slopes)

        print("Basin: $i, slope: $(basin_srtm15plus_aveslope[i]) \n")
    end

    # save("basin_srtm15plus_aveslope.jld","basin_srtm15plus_aveslope",basin_srtm15plus_aveslope)
    # basin_srtm15plus_aveslope =  load("basin_srtm15plus_aveslope.jld")["basin_srtm15plus_aveslope"]

    # save("OctopusSlopeRecalc.jld","basin_srtm15plus_aveslope",basin_srtm15plus_aveslope,"basin_polygon_n",basin_polygon_n,"basin_polygon_lat",basin_polygon_lat,"basin_polygon_lon",basin_polygon_lon,"subbasins",subbasins)


## --- Plot raw Lat and Lon
    plot(data["x_wgs84"],data["y_wgs84"],seriestype=:scatter)
    plot!(xlabel="Longitude", ylabel="Latitude")

## --- Plot raw erosion rate vs slope
    plot(data["slp_ave"],data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS 10Be", fg_color_legend=:white)
    plot!(xlabel="Slope (m/km)", ylabel="Erosion rate (mm/kyr)",yscale=:log10)
    # savefig("Slope_vs_Erosion.pdf")

## --- Plot new erosion rate vs slope

    plot(basin_srtm15plus_aveslope,data["ebe_mmkyr"],seriestype=:scatter,label="OCTOPUS 10Be", fg_color_legend=:white)
    plot!(xlabel="SRTM15 Slope (m/km)", ylabel="Erosion rate (mm/kyr)",yscale=:log10)
    # savefig("Slope_vs_Erosion.pdf")

    # log10(E) = slope/140 + 0.7
    # E = 10^(slope/140 + 0.7)

## --- Fit raw erosion rate as a function of slope (m/km)

    using LsqFit: curve_fit

    function linear(x,p)
        y = p[1] + x * p[2]
    end
    p = [0.5, 1/100]
    t = .~isnan.(basin_srtm15plus_aveslope) .& .~isnan.(data["ebe_mmkyr"]) .& (basin_srtm15plus_aveslope .< 300)
    x = basin_srtm15plus_aveslope[t]
    y = log10.(data["ebe_mmkyr"])[t]
    fobj = curve_fit(linear, x, y, p);

    # fobj.param[1]: 0.987237
    # fobj.param[2]: 0.00555159
    function Emmkyr(slp)
        return 10^(slp*0.00567517 + 0.971075)
    end

    plot!(1:500, Emmkyr.(1:500), label = "fit")

## ---

    # nbins = 20;
    # nsims = 10000;
    # (c, m, el, eu) = bin_bsr_means(data["slp_ave"], log10.(data["ebe_mmkyr"]), 0, 800, nbins, data["slp_std"]/10, nsims)
    # plot(c, m, yerror=(el,eu), seriestype=:scatter, label="",xlabel="Slope", ylabel="Log10 Erosion rate (mm/kyr)")

    (c,m,sigma) = binmeans(data["slp_ave"], log10.(data["ebe_mmkyr"]), 0, 800, 20)
    plot(c,m,yerror=2*sigma,seriestype="scatter")
    plot!(xlabel="Slope", ylabel="Log10 Erosion rate (mm/kyr)")

## --- Plot raw erosion rate vs elevation
    plot(data["elev_ave"],data["ebe_mmkyr"],seriestype=:scatter)
    plot!(xlabel="Elevation (m)", ylabel="Erosion rate (mm/kyr)",yscale=:log10)

## ---
    (c,m,sigma) = binmeans(data["elev_ave"], data["ebe_mmkyr"], 0, 6000, 12)
    plot(c,m,yerror=2*sigma,seriestype="scatter")
    plot!(xlabel="Elevation (m)", ylabel="Erosion rate (mm/kyr)")

## --- Plot raw erosion rate vs elevation scaled to snowline
    plot(data["elev_ave"]./data["Snowline"],data["ebe_mmkyr"],seriestype=:scatter)
    plot!(xlabel="Elevation/ELA", ylabel="Erosion rate (mm/kyr)",yscale=:log10)

## ---
    (c,m,sigma) = binmeans(data["elev_ave"]./data["Snowline"], data["ebe_mmkyr"], 0, 1.2, 12)
    plot(c,m,yerror=2*sigma,seriestype="scatter")
    plot!(xlabel="Elevation/ELA", ylabel="Erosion rate (mm/kyr)")

## --- Plot resampled elevation vs erosion rate
    nbins = 20;
    nsims = 10000;
    (c, m, el, eu) = bin_bsr_means(data["elev_ave"], data["ebe_mmkyr"], 0, 6000, nbins, data["elev_std"], nsims)
    plot(c, m, yerror=(el,eu), seriestype=:scatter, label="",xlabel="Elevation (m)", ylabel="Erosion rate (mm/kyr)")

## --- Plot resampled elevation / snowline vs erosion rate
    nbins = 20;
    nsims = 10000;
    (c, m, el, eu) = bin_bsr_means(data["elev_ave"]./data["Snowline"], data["ebe_mmkyr"], 0, 1.2, nbins, data["elev_std"]./data["Snowline"], nsims)
    plot(c, m, yerror=(el,eu), seriestype=:scatter, label="")
    plot!(xlabel="Elevation/ELA", ylabel="Erosion rate (mm/kyr)")

## --- End of File
