## --- Define function for picking random points

function random_lat_lon(n)
    randlon = rand(n)*360-180
    randlat = 90 - acos.(rand(n)*2-1)*180/pi

    return (randlat, randlon)
end

## --- Define functions for querying macrostrat

using HTTP, JSON

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
