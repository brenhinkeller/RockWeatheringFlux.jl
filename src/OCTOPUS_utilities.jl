## --- Unzip the kml file
    """
    ```julia
    load_octopus_kml(filename)
    ```

    """
    function load_octopus_kml(filename)
        # Check number and ordering of coordinate lists
        fid = open("data/octopus/crn_basins_global.kml")

        # Count the number of lines in the file
        i = 0
        for line in eachline(fid)
            i += 1
        end
        isfirstcoord = Array{Bool}(undef, i)

        # Parse the coordinates (basin polygon outlines) in each placemark
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

## --- Parse the file
    """
    ```julia
    parse_octopus_kml_variables(str)
    ```
    """
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

        needsNaNs = ["alcorr", "be10ep", "beprod", "eal_err", "elev_std", "errbe_prod", "al26ep", "alprod", "be10ep_err", 
            "beself", "eal_gcmyr", "erral_ams", "errbe_tot", "al26ep_err", "alself", "be10nc", "besnow", "eal_mmkyr", 
            "erral_muon", "sizemax", "al26nc", "alsnow", "be10nc_err", "betopo", "ebe_err", "erral_prod", "sizemin", 
            "al26nc_err", "altopo", "be10np", "betots", "ebe_gcmyr", "erral_tot", "slp_ave", "al26np", "altots", "be10np_err", 
            "ebe_mmkyr", "errbe_ams", "slp_std", "al26np_err", "area", "becorr", "dbver", "elev_ave", "errbe_muon", ]

        for i = 1:length(needsNaNs)
            data[needsNaNs[i]][data[needsNaNs[i]] .< 0] .= NaN
        end

        return data
    end

## --- Parse the basin polygon outlines
    """
    ```julia
    parse_octopus_polygon_outlines(str,isfirstcoord)
    ```
    """
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

## --- Calculate precipitation
    """
    ```julia
    find_precip(lat,lon)
    ```
    for each basin? or each coordinate?
    """
    function find_precip(lat,lon)
        prate = ncread("data/prate.sfc.mon.ltm.nc","prate")
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

## --- (Re?)Calculate slope for each basin
    """
    ```julia
    get_basin_srtm15plus_aveslope(srtm::Dict,nbasins,subbasins,basin_polygon_lat,basin_polygon_lon)
    ```
    """
    # For the future, maybe I want this to be max slope?
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
                (column, row) = find_grid_inpolygon(x_lon_cntr, y_lat_cntr, subbasin_lon, subbasin_lat)
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

## --- Cutve fitting functions for slope / erosion rate
    """
    ```julia
    linear(x,p)
    ```
    """
    function linear(x,p)
        y = p[1] .+ x * p[2]
    end