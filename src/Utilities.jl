## --- Generate random points on the continental crust
"""
```julia
gen_continental_points(npoints, etopo)
```

Generate `npoints` element-long lists of latitudes, longitudes, and elevations for points 
uniformly distributed across the continental crust. Requires `etopo` matrix of 1 arc minute 
resolution global relief model. 

Units are meters of elevation and decimal degrees of latitude and longitude.

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

        # Find points above sea level
        elev = find_etopoelev(etopo,randlat,randlon)
        abovesea = elev .> 0
        newpoints = min(count(abovesea), length(rocklat) - currentpoints)

        # Concatenate together all the points that represent exposed crust
        rocklat[(currentpoints+1):(currentpoints+newpoints)] = randlat[abovesea][1:newpoints]
        rocklon[(currentpoints+1):(currentpoints+newpoints)] = randlon[abovesea][1:newpoints]
        elevations[(currentpoints+1):(currentpoints+newpoints)] = elev[abovesea][1:newpoints]

        currentpoints += newpoints
    end

    return rocklat, rocklon, elevations
end

## --- Slope and erosion rate relationship
"""
```julia
emmkyr(slp)
```

Find the erosion rate in mm/kyr given a slope `slp`.
"""
    function emmkyr(slp)
        return 10^(slp*0.00567517 + 0.971075)
    end


## --- Functions for querying macrostrat

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


## --- Functions for dealing with SRTM15










## --- Functions for dealing with SRTM15
# I think I can get rid of these and use the statgeochem ones, or at least a modification of them...
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

    function find_srtm15plus_aveslope_around(srtm15plus,lat,lon; halfwidth=1::Integer, max_allowed_slope=1000::Number)

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
                for r = (row-halfwidth):(row+halfwidth)
                    for c = (col-halfwidth):(col+halfwidth)
                        if r>0 && r<43202
                            res = data[r,mod(c-1,86400)+1]
                            if res < max_allowed_slope
                                k +=1
                                out[i] += res
                            end
                        end
                    end
                end
                out[i] /= k
            end
        end

        return out
    end


## --- Find the average value of slope over an area
"""
```julia
avg_over_area(data::Matrix, lat::Vector, lon::Vector, sf::Number=240; 
    halfwidth::Number=1, 
    maxpossible::Number=0xffff)
```
Find the average value of `data` over an area with radius `halfwidth` arc-seconds
at coordinates `lat` and `lon`.

This is distinct from `StatGeochem`'s `aveslope`. This function finds the average over an area, not the
average slope for a specific point. For example, this might be given a matrix of maximum slope at each point
on Earth, and would return the _average maximum slope_ at each point of interest.

Other optional arguments and defaults:

    sf::Number=240

Scale factor (cells per degree) for the SRTM15+ data. For 15 arc-second resolution, the scale factor is 240, 
because 15 arc-seconds go into 1 arc-degree 60 * 4 = 240 times.

    maxpossible::Number=0xffff

The maximum possible value for the variable of interest; variables with values greater than this are ignored.

"""
function avg_over_area(data::Matrix, lat::Vector, lon::Vector, sf::Number=240;
    halfwidth::Number=1, 
    maxpossible::Number=0xffff)
    # Make sure we will never index out of bounds
    @assert eachindex(lat) == eachindex(lon)

    # Make sure data has values that cover all of Earth
    (nrows, ncols) = size(data)
    @assert nrows == 180 * sf + 1   # Why 180 instead of 360?
    @assert ncols == 360 * sf + 1
    

    # Preallocate
    out = fill(NaN, length(lat))

    # Find result by indexing into the varname matrix
    for i in eachindex(lat)
        if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
            # Result is NaN if either input is NaN or out of bounds
            continue

        else
            # Convert latitude and longitude into indicies 
            row = 1 + round(Int,(90 + lat[i]) * sf)
            col = 1 + round(Int,(180 + lon[i]) * sf)

            # Index into the array - could I use @turbo here? It currently gets mad
            k = 0           # Generic counter
            out[i] = 0      # Starting value
            for r = (row-halfwidth):(row+halfwidth)
                for c = (col-halfwidth):(col+halfwidth)
                    # Only do the computation if we are in bounds
                    if 1 <= r <= nrows
                        res = data[r, mod(c-1,ncols-1)+1]

                        # Ignore values that are larger than the max possible value
                        if res < maxpossible
                            k +=1
                            out[i] += res
                        end
                    end
                end
            end

            # Save the average value
            out[i] /= k
        end
    end

    return out
end

"""
```julia
get_stats(data)
```
Ignoring `NaN`s, calculate the sum, mean, and standard deviation of `data`.

### Example
```
(data_s, data_m, data_e) = get_stats(data)
```
"""
function get_stats(data)
    return nansum(data), nanmean(data), nanstd(data)
end

## --- End of file