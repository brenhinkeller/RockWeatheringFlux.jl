# Get and parse data from the Macrostrat API

## --- Generate random points on the continental crust
    """
    ```julia
    gen_continental_points(npoints, etopo)
    ```

    Generate an `npoints` element-long lists of latitudes, longitudes, and elevations for
    points randomly and uniformly distributed across the continental crust. 

    Requires `etopo` matrix of 1 arc minute resolution global relief model. Units are 
    meters of elevation and decimal degrees of latitude and longitude.

    See also: `get_etopo`.
    """
    function gen_continental_points(npoints, etopo)
        # Initialize
        rocklat = Array{Float64}(undef, npoints)
        rocklon = Array{Float64}(undef, npoints)
        elevations = Array{Float64}(undef, npoints)

        # Number of points added to lat / lon arrays
        currentpoints = 0

        while currentpoints < npoints
            # Generate some random latitudes and longitudes with uniform spatial density on the globe
            (randlat, randlon) = randlatlon(length(rocklat))

            # Find non-OIB points above sea level
            elev = find_etopoelev(etopo,randlat,randlon)
            isOIB = findOIBs(randlat, randlon)
            continental = (elev .> 0) .& .!isOIB
            newpoints = min(count(continental), length(rocklat) - currentpoints)

            # Concatenate together all the points that represent exposed crust
            rocklat[(currentpoints+1):(currentpoints+newpoints)] = randlat[continental][1:newpoints]
            rocklon[(currentpoints+1):(currentpoints+newpoints)] = randlon[continental][1:newpoints]
            elevations[(currentpoints+1):(currentpoints+newpoints)] = elev[continental][1:newpoints]

            currentpoints += newpoints
        end

        return rocklat, rocklon, elevations
    end
    export gen_continental_points


## --- Remove ocean islands
    """
    ```julia
    findOIBs(lats, lons)
    ```

    Return a BitVector that is `true` for any coordinates corresponding to known ocean 
    island hot-spots.

    Currently identifies: Azores, Balleny, Bowie, Canary, Cape Verde, Cobb, Discovery, 
    Easter, Galapagos, Hawaii, Iceland, Jan Mayan, Juan Fernandez, Kerguelen, Loisville,
    Pitcairn, Reunion, Saint Helena, Sao Tome, and Tristan) & (Gough.

    """
    function findOIBs(lats::AbstractArray, lons::AbstractArray)
        OIBs = @. (
            ((lats > 18) & (lats < 30) & (lons < -154) & (lons > -179.5)) |     # Hawaii
            ((lats > 36) & (lats < 41) & (lons < -22) & (lons > -29)) |         # Azores
            ((lats > -69) & (lats < -65) & (lons < 170) & (lons > 160)) |       # Balleny
            ((lats > 53) & (lats < 57) & (lons < -135) & (lons > -147)) |       # Bowie
            ((lats > 27) & (lats < 29) & (lons < -15) & (lons > -18.5)) |       # Canary
            ((lats > 14) & (lats < 18) & (lons < -22) & (lons > -26)) |         # Cape Verde
            ((lats > 46) & (lats < 50) & (lons < -130) & (lons > -135)) |       # Cobb
            ((lats > -44) & (lats < -41) & (lons < 3) & (lons > -5)) |          # Discovery
            ((lats > -30) & (lats < -20) & (lons < -80) & (lons > -110)) |      # Easter
            ((lats > -3) & (lats < 1) & (lons < -82) & (lons > -92)) |          # Galapagos
            ((lats > 62) & (lats < 67) & (lons < -13) & (lons > -26)) |         # Iceland
            ((lats > 68.5) & (lats < 71.5) & (lons < -7) & (lons > -10)) |      # Jan Mayen
            ((lats > -34.3) & (lats < -33) & (lons < -76) & (lons > -82)) |     # Juan Fernandez
            ((lats > -63) & (lats < -45) & (lons < 85) & (lons > 63)) |         # Kerguelen
            ((lats > -50) & (lats < -25) & (lons < -145) & (lons > -175)) |     # Louisville
            ((lats > -25) & (lats < -20) & (lons < -125) & (lons > -140)) |     # Pitcairn
            ((lats > -22) & (lats < -3) & (lons < 63) & (lons > 55)) |          # Reunion
            ((lats > -6) & (lats < -5.5) & (lons < -15.75) & (lons > -16.25)) | # Saint Helena
            ((lats > 0) & (lats < 4) & (lons < 9) & (lons > 6)) |               # Sao Tome
            ((lats > -40) & (lats < -38) & (lons < -9) & (lons > -13)) |        # Tristan & Gough
            ((lats > -35) & (lats < -19) & (lons < 11) & (lons > -5))           # Off the coast of Namibia?
        )
        return OIBs
    end
    export findOIBs


## --- Ping Macrostrat API
    """
    ```julia
    query_macrostrat(lat, lon)
    ```
    Get lithological data at highest available resolution for rocks at `lat`, `lon` 
    coordinate from the Macrostrat geologic map compilation v2 API.

    """
    function query_macrostrat(lat, lon)
        # Check highest possible resolution
        resp = HTTP.get("https://macrostrat.org/api/v2/geologic_units/map?lat=$lat&lng=$lon&scale=large")
        str = String(resp.body)
        parsed = JSON.Parser.parse(str)

        # If no data, check medium scale resolution 
        if isempty(parsed["success"]["data"])
            resp = HTTP.get("https://macrostrat.org/api/v2/geologic_units/map?lat=$lat&lng=$lon&scale=medium")
            str = String(resp.body)
            parsed = JSON.Parser.parse(str)
        else
            return parsed
        end

        # If still no data, check small scale resolution
        if isempty(parsed["success"]["data"])
            resp = HTTP.get("https://macrostrat.org/api/v2/geologic_units/map?lat=$lat&lng=$lon&scale=small")
            str = String(resp.body)
            parsed = JSON.Parser.parse(str)
        end

        # Fallback: global geologic map
        if isempty(parsed["success"]["data"])
            resp = HTTP.get("https://macrostrat.org/api/v2/geologic_units/map?lat=$lat&lng=$lon&scale=tiny")
            str = String(resp.body)
            parsed = JSON.Parser.parse(str)
        end

        return parsed 
    end
    export query_macrostrat


## --- Get data from the unparsed response dictionary
    function get_macrostrat_min_age(jobj)
        try
            return jobj["success"]["data"][1]["t_int_age"]::Number
        catch error
            return NaN
        end
    end
    export get_macrostrat_min_age

    function get_macrostrat_max_age(jobj)
        try
            return jobj["success"]["data"][1]["b_int_age"]::Number
        catch error
            return NaN
        end
    end
    export get_macrostrat_max_age

    function get_macrostrat_map_id(jobj)
        try
            return jobj["success"]["data"][1]["map_id"]::Number
        catch error
            return NaN
        end
    end
    export get_macrostrat_map_id

    function get_macrostrat_lith(jobj)
        try
            return jobj["success"]["data"][1]["lith"]
        catch error
            return "NA"
        end
    end
    export get_macrostrat_lith

    function get_macrostrat_descrip(jobj)
        try
            return jobj["success"]["data"][1]["descrip"]
        catch error
            return "NA"
        end
    end
    export get_macrostrat_descrip

    function get_macrostrat_name(jobj)
        try
            return jobj["success"]["data"][1]["name"]
        catch error
            return "NA"
        end
    end
    export get_macrostrat_name

    function get_macrostrat_strat_name(jobj)
        try
            return jobj["success"]["data"][1]["strat_name"]
        catch error
            return "NA"
        end
    end
    export get_macrostrat_strat_name

    function get_macrostrat_comments(jobj)
        try
            return jobj["success"]["data"][1]["comments"]
        catch error
            return "NA"
        end
    end
    export get_macrostrat_comments

    function get_macrostrat_refs(jobj)
        try
            return join(collect(values(jobj["success"]["refs"])), " | ")
        catch error
            return "NA"
        end
    end
    export get_macrostrat_refs


## --- Parse responses
    function parse_macrostrat_responses(responses, stop)
        # Preallocate
        rocktype = Array{String}(undef, stop, 1)
        rockdescrip = Array{String}(undef, stop, 1)
        rockname = Array{String}(undef, stop, 1)
        rockstratname = Array{String}(undef, stop, 1)
        rockcomments = Array{String}(undef, stop, 1)
        agemax = Array{Float64}(undef, stop, 1)
        agemin = Array{Float64}(undef, stop, 1)
        age = Array{Float64}(undef, stop, 1)
        refs = Array{String}(undef, stop, 1)

        # Parse responses into preallocated arrays
        for i in 1:stop
            # Catch undefined references
            rocktype[i] = get_macrostrat_lith(responses[i])
            rockdescrip[i] = get_macrostrat_descrip(responses[i])
            rockname[i] = get_macrostrat_name(responses[i])
            rockstratname[i] = get_macrostrat_strat_name(responses[i])
            rockcomments[i] = get_macrostrat_comments(responses[i])
            agemax[i] = get_macrostrat_max_age(responses[i])
            agemin[i] = get_macrostrat_min_age(responses[i])

            refs[i] = get_macrostrat_refs(responses[i])
        end

        # Filter ages younger than 0 and older than 4000 and compute age from min / max
        # Put age bounds in the correct order
        for i in 1:stop
            agemax[i] = ifelse(0<agemax[i]<4000, agemax[i], NaN)
            agemin[i] = ifelse(0<agemin[i]<4000, agemin[i], NaN)
            age[i] = nanmean([agemax[i], agemin[i]])

            if agemin[i] > agemax[i]
                agemin[i] = agemin[i] + agemax[i]
                agemax[i] = agemin[i] - agemax[i]
                agemin[i] = agemin[i] - agemax[i]
            end
        end

        # Convert strings to lowercase so they can be matched to known names of rock types
        rocktype = lowercase.(rocktype)
        rockdescrip = lowercase.(rockdescrip)
        rockname = lowercase.(rockname)
        rockstratname = lowercase.(rockstratname)
        rockcomments = lowercase.(rockcomments)

        # Replace tabs with spaces so they will not be confused with the \t delimitator
        rocktype = replace.(rocktype, "    " => " ")
        rockdescrip = replace.(rockdescrip, "    " => " ")
        rockname = replace.(rockname, "    " => " ")
        rockstratname = replace.(rockstratname, "    " => " ")
        rockcomments = replace.(rockcomments, "    " => " ")

        # Return as a tuple
        return (
            agemax = agemax,
            agemin = agemin,
            age = age,
            rocktype = rocktype,
            rockname = rockname,
            rockdescrip = rockdescrip,
            rockstratname = rockstratname,
            rockcomments = rockcomments,
            refstrings = refs
        )
    end
    export parse_macrostrat_responses


## --- Crop intermediate files
    """
    ```julia
    crop_intermediate!(pathin::String, [pathout::String])
    ```

    Intermediate Macrostrat save files contain all pre-generated latitudes, longitudes, 
    and elevations. Reduce the latitude / longitude / elevations saved in the file to
    only those with responses.

    File names must be relative paths, including file name extensions.
    """
    function crop_intermediate!(pathin::String, pathout::String=pathin)
        # Load the old file
        fid = h5open("$pathin", "r")
            rocklat = read(fid["rocklat"])
            rocklon = read(fid["rocklon"])
            elevation = read(fid["elevation"])
            agemax = read(fid["agemax"])
            agemin = read(fid["agemin"])
            age = read(fid["age"])
            rocktype = read(fid["rocktype"])
            rockname = read(fid["rockname"])
            rockdescrip = read(fid["rockdescrip"])
            rockstratname = read(fid["rockstratname"])
            rockcomments = read(fid["rockcomments"])
            refstrings = read(fid["reference"])
        close(fid)

        # Get number of samples
        npoints = length(rocktype)

        # Get the rock type of each sample
        macro_cats = match_rocktype(rocktype, rockname, rockdescrip)
        types = Array{String}(undef, npoints, 1)
        for i in eachindex(types)
            type = get_type(macro_cats, i)
            type==nothing && (type="")
            types[i] = string(type)
        end

        # Load the new file
        fid = h5open("$pathout", "w")
            fid["rocklat"] = rocklat[1:npoints]
            fid["rocklon"] = rocklon[1:npoints]
            fid["elevation"] = elevation[1:npoints]
            fid["agemax"] = agemax
            fid["agemin"] = agemin
            fid["age"] = age
            fid["rocktype"] = rocktype
            fid["rockname"] = rockname
            fid["rockdescrip"] = rockdescrip
            fid["rockstratname"] = rockstratname
            fid["rockcomments"] = rockcomments
            fid["reference"] = refstrings
            fid["npoints"] = npoints
            fid["typecategory"] = types
        close(fid)
    end
    export crop_intermediate!

## --- End of File