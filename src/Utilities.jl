## --- Define function for finding elevation

function findelevation(etopoelev,lat,lon)
    # elev=FINDELEVATION(lat,lon,[etopoelev])
    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the ETOPO elevation model.

    # Scale factor used in map (=nrows/180=ncols/360)
    sf=60;
    maxrow = 180*sf;
    maxcol = 360*sf;

    # Create and fill output vector
    elev=Array{Float64}(length(lat));
    for i=1:length(lat)
        if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
            elev[i]=NaN; # Result is NaN if either input is NaN
        else
            # Convert latitude and longitude into indicies of the elevation map array
            row=round(Int,(90+lat[i])*sf+0.5);
            if row == (maxrow+1)
                row = maxrow;
            end

            col=round(Int,(180+lon[i])*sf+0.5);
            if col == (maxcol+1)
                col = maxcol;
            end

            elev[i]=etopoelev[row,col]; # Otherwise, find result
        end
    end

        return elev
end

## --- Define function for picking random points

function random_lat_lon(n)
    randlon = rand(n)*360-180
    randlat = 90 - acos.(rand(n)*2-1)*180/pi

    # randlat = Array{Float64}(n)
    # for i = 1:n
    #     while true
    #         proposed_lat = rand()*180-90
    #         if rand() < sin.((proposed_lat+90)*pi/180)
    #             randlat[i] = proposed_lat;
    #             break;
    #         end
    #     end
    # end

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
