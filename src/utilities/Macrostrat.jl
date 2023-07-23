## --- Ping Macrostrat API
    """
    ```julia
    query_macrostrat(lat, lon, zoom::Number=11)
    ```
    Get lithological data for rocks at `lat`, `lon` coordinate from the Macrostrat API.

    Argument `zoom` controls precision; default is approximately 5km. Automatically retry with
    less precise window if initial query does not return data.
    """
    function query_macrostrat(lat, lon, zoom::Number=11)
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


## --- Parse Macrostrat responses
    function get_macrostrat_min_age(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["t_int_age"]::Number
        catch error
            return NaN
        end
    end

    function get_macrostrat_max_age(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["b_int_age"]::Number
        catch error
            return NaN
        end
    end

    function get_macrostrat_map_id(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["map_id"]::Number
        catch error
            return NaN
        end
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

    function get_macrostrat_ref_title(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["ref_title"]
        catch error
            return "NA"
        end
    end

    function get_macrostrat_ref_authors(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["authors"]
        catch error
            return "NA"
        end
    end

    function get_macrostrat_ref_year(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["ref_year"]
        catch error
            return "NA"
        end
    end

    function get_macrostrat_ref_doi(jobj)
        try
            return jobj["success"]["data"]["burwell"][1]["ref"]["isbn_doi"]
        catch error
            return "NA"
        end
    end


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
        authors = Array{String}(undef, stop, 1)
        years = Array{String}(undef, stop, 1)
        titles = Array{String}(undef, stop, 1)
        dois = Array{String}(undef, stop, 1)

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

            authors[i] = get_macrostrat_ref_authors(responses[i])
            years[i] =get_macrostrat_ref_year(responses[i])
            titles[i] = get_macrostrat_ref_title(responses[i])
            dois[i] = get_macrostrat_ref_doi(responses[i])
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

        # Condense references
        refstrings = @. authors * " | " * years * " | " * titles * " | " * dois

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
            refstrings = refstrings
        )
    end


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

## --- End of File