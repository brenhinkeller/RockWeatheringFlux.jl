## --- Define function for picking random points

    function random_lat_lon(n)
        randlon = rand(n)*360-180
        randlat = 90 - acos.(rand(n)*2-1)*180/pi

        return (randlat, randlon)
    end


## --- Define functions for querying macrostrat

    function query_macrostrat(lat, lon, zoom)
        resp = HTTP.get("https://macrostrat.org/api/mobile/map_query?lat=$lat&lng=$lon&z=$zoom")
        str = String(resp.body)
        parsed = JSON.Parser.parse(str)
        try
            parsed["success"]["data"]["burwell"][1]["lith"]
        catch error
            resp = HTTP.get("https://macrostrat.org/api/mobile/map_query?lat=$lat&lng=$lon&z=1")
            str = String(resp.body)
            parsed = JSON.Parser.parse(str)
        end
        return parsed
    end

    function get_macrostrat_lith(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["lith"]
        catch error
            return "NA"
        end
    end

    function get_macrostrat_descrip(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["descrip"]
        catch error
            return "NA"
        end
    end

    function get_macrostrat_name(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["name"]
        catch error
            return "NA"
        end
    end

    function get_macrostrat_strat_name(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["strat_name"]
        catch error
            return "NA"
        end
    end

    function get_macrostrat_comments(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["comments"]
        catch error
            return "NA"
        end
    end

## --- Functions for dealing with SRTM15
    resourcepath = "data"

    # Read srtm15plus file from HDF5 storage, downloading from cloud if necessary
    function get_srtm15plus_aveslope(varname="")
        # Available variable names: "slope", "y_lat_cntr", "x_lon_cntr",
        # "nanval", "cellsize", "scalefactor", and "reference". Units are
        # meters of elevation and decimal degrees of latitude and longitude

        # Construct file path
        filepath = joinpath(resourcepath,"srtm15plus_aveslope.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            print("Downloading srtm15plus.h5 from google cloud storage\n")
            download("https://storage.googleapis.com/statgeochem/srtm15plus_aveslope.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath, "vars/"*varname)
    end

    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the SRTM15plus 15-arc-second elevation model.
    function find_srtm15plus_aveslope(srtm15plus,lat,lon)

        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        elseif isa(srtm15plus,Dict)
            data = srtm15plus["slope"]
        elseif isa(srtm15plus, Array)
            data = srtm15plus
        else
            error("wrong srtm15plus variable")
        end

        # Scale factor (cells per degree) = 60 * 4 = 240
        # (15 arc seconds goes into 1 arc degree 240 times)
        sf = 240

        # Create and fill output vector
        out=Array{Float64}(size(lat));
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                out[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                # Note that STRTM15 plus has N+1 columns where N = 360*sf
                row = 1 + round(Int,(90+lat[i])*sf)
                col = 1 + round(Int,(180+lon[i])*sf)
                # Find result by indexing
                res = data[row,col]
                if res > 1000
                    out[i] = NaN
                else
                    out[i] = res
                end
            end
        end

        return out
    end

    function find_srtm15plus_aveslope_around(srtm15plus,lat,lon; halfwidth_lat=10, halfwidth_lon=10)

        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        elseif isa(srtm15plus,Dict)
            data = srtm15plus["slope"]
        elseif isa(srtm15plus, Array)
            data = srtm15plus
        else
            error("wrong srtm15plus variable")
        end

        # Scale factor (cells per degree) = 60 * 4 = 240
        # (15 arc seconds goes into 1 arc degree 240 times)
        sf = 240

        # Create and fill output vector
        out=Array{Float64}(size(lat));
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                out[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                # Note that STRTM15 plus has N+1 columns where N = 360*sf
                row = 1 + round(Int,(90+lat[i])*sf)
                col = 1 + round(Int,(180+lon[i])*sf)

                # Find result by indexing
                k = 0;
                out[i] = 0;
                for r = row-halfwidth_lat:row+halfwidth_lat
                    for c = col-halfwidth_lon:col+halfwidth_lon
                        if r>0 && r<43202
                            k +=1
                            out[i] += data[r,mod(c-1,86400)+1]
                        end
                    end
                end
                out[i] /= k
            end
        end

        return out
    end
