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
Find the average value of `data` over an area with radius `halfwidth` (units of gridcells) at 
coordinates `lat` and `lon`.

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


"""
```julia
match_rocktype(rocktype, rockname, rockdescrip; major=false)
```
Match samples to our rock types from the Burwell `rocktype`, `rockname`, and `rockdescrip`. If `major` is `true`, 
only returns matches for `sed`, `ign`, `met`, and `cover` rocks. If `major` is `false`, returns:
```
sed, ign, met, volc, plut, hypabyssal, metaign, metased, lowgrade, highgrade, cover
```
"""
function match_rocktype(rocktype, rockname, rockdescrip; major=false)
    npoints = length(rocktype)

    # Define our rock types based on names or partial names of different rock types from GetBurwellBulkAge.m
    # Sedimentary
    sedtypes = ["sediment", "fluv", " clast", "siliciclast", "conglomerat", "gravel", "sand", "psamm", "arenit", "arkos", "silt", 
    "mud", "marl", "clay", "shale", "wacke", "argillite", "argillaceous", "pelit", "pebble", "carbonate", "limestone", "dolo", 
    "caliche", "chalk", "travertine", "tavertine", "teravertine", "tufa", "evaporite", " salt", "salt flat", "gypsum", "boulder", 
    "diamict", "tillite", "stream", "beach", "terrace", "chert", "banded iron", "coal", "anthracite", "marine deposits", "turbidite", 
    "flysch", "paleosol"]

    # Igneous
    volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme", "pipe", "basalt", "andesit", "dacit", "rhyolit", 
    "pillow", "carbonatite", "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", "fonolito", 
    "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", "komatiite",]
    pluttypes = ["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", "monzonit", "syenit", "peridot", "dunit", 
    "harzburg", "anorthosite", "mangerite", "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", 
    "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",]
    hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", "dolerit", "diabase", "porphyry", 
    "microgranite"]
    
    igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",],volctypes,pluttypes,hypabyssaltypes)

    # Metamorphic
    metasedtypes = ["para", "metased", "meta-sed", "schist", "quartzite", "marble", "skarn", "slate", "phyllite",]
    metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone", "eclogite", "metabasite", "migma",]
    lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", "gossan", "alter", "hydrothermal", 
    "palagonite",]
    highgradetypes = ["crystalline", "basement", "marble", "skarn", "blueschist", "gneiss", "amphibolit", "eclogite", "granulit", 
    "hornfels", "granofels", "sanidinite", "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
    "harzburg", "high grade metamorphic"]
    
    mettypes = vcat(["gneiss", "granulit", "hornfels", "granofels", "mylonit", "meta", "cataclasite", "melange", "gouge", "tecton", 
    "calc silicate"], metasedtypes, metaigntypes, lowgradetypes, highgradetypes)

    # Other
    covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", "laterite", "surficial deposits", "talus", 
        "scree", "mass-wasting", "slide", "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", "debris"]
    cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]

    if major
        # Preallocate
        sed = falses(npoints)
        ign = falses(npoints)
        met = falses(npoints)
        cover = falses(npoints)

        # Check major lithology first
        for i = 1:length(sedtypes)
            sed .|= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], sedtypes[i]) : false)
        end
        for i = 1:length(igntypes)
            ign .|= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], igntypes[i]) : false)
        end
        for i = 1:length(mettypes)
            met .|= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], mettypes[i]) : false)
        end
        for i = 1:length(covertypes)
            cover .|= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], covertypes[i]) : false)
        end

        # Check the rest of rocktype
        not_matched = .~(sed .| ign .| met .| cover)
        for i = 1:length(sedtypes)
            sed[not_matched] .|= containsi.(rocktype[not_matched],sedtypes[i])
        end
        for i = 1:length(igntypes)
            ign[not_matched] .|= containsi.(rocktype[not_matched],igntypes[i])
        end
        for i = 1:length(mettypes)
            met[not_matched] .|= containsi.(rocktype[not_matched],mettypes[i])
        end
        for i = 1:length(covertypes)
            cover[not_matched] .|= containsi.(rocktype[not_matched],covertypes[i])
        end

        # Then rockname
        not_matched = .~(sed .| ign .| met .| cover)
        for i = 1:length(sedtypes)
            sed[not_matched] .|= containsi.(rockname[not_matched],sedtypes[i])
        end
        for i = 1:length(igntypes)
            ign[not_matched] .|= containsi.(rockname[not_matched],igntypes[i])
        end
        for i = 1:length(mettypes)
            met[not_matched] .|= containsi.(rockname[not_matched],mettypes[i])
        end
            for i = 1:length(covertypes)
            cover[not_matched] .|= containsi.(rockname[not_matched],covertypes[i])
        end

        # Then rockdescrip
        not_matched = .~(sed .| ign .| met .| cover)
        for i = 1:length(sedtypes)
            sed[not_matched] .|= sed[not_matched] .| containsi.(rockdescrip[not_matched],sedtypes[i])
        end
        for i = 1:length(igntypes)
            ign[not_matched] .|= containsi.(rockdescrip[not_matched],igntypes[i])
        end
        for i = 1:length(mettypes)
            met[not_matched] .|= containsi.(rockdescrip[not_matched],mettypes[i])
        end
        for i = 1:length(covertypes)
        cover[not_matched] .|= containsi.(rockdescrip[not_matched],covertypes[i])
        end

        return sed, ign, met, cover
    else
        # Preallocate arrays for each rock type
        sed = falses(npoints)
        ign = falses(npoints)
        met = falses(npoints)
        volc = falses(npoints)
        plut = falses(npoints)
        hypabyssal = falses(npoints)
        metaign = falses(npoints)
        metased = falses(npoints)
        lowgrade = falses(npoints)
        highgrade = falses(npoints)
        cover = falses(npoints)

        # Check all possible types, starting with major lithology
        for i = 1:length(sedtypes)
            sed .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], sedtypes[i]) : false)
        end
        for i = 1:length(igntypes)
            ign .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], igntypes[i]) : false)
        end
        for i = 1:length(mettypes)
            met .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], mettypes[i]) : false)
        end
        for i = 1:length(covertypes)
            cover .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], covertypes[i]) : false)
        end
        
        for i = 1:length(volctypes)
            volc .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], volctypes[i]) : false )
        end
        for i = 1:length(pluttypes)
            plut .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], pluttypes[i]) : false )
        end
        for i = 1:length(metaigntypes)
            metaign .|= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], metaigntypes[i]) : false)
        end
        for i = 1:length(metasedtypes)
            metased .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], metasedtypes[i]) : false)
        end
        for i = 1:length(lowgradetypes)
            lowgrade .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], lowgradetypes[i]) : false)
        end
        for i = 1:length(highgradetypes)
            highgrade .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], highgradetypes[i]) : false)
        end
        for i = 1:length(hypabyssaltypes)
            hypabyssal .|=  (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], hypabyssaltypes[i]) : false)
        end

        # Check the rest of rocktype
        not_matched = .~(sed .| ign .| met .| cover)
        for i = 1:length(sedtypes)
            sed[not_matched] .|= containsi.(rocktype[not_matched],sedtypes[i])
        end
        for i = 1:length(igntypes)
            ign[not_matched] .|= containsi.(rocktype[not_matched],igntypes[i])
        end
        for i = 1:length(mettypes)
            met[not_matched] .|= containsi.(rocktype[not_matched],mettypes[i])
        end
        for i = 1:length(covertypes)
            cover[not_matched] .|= containsi.(rocktype[not_matched],covertypes[i])
        end

        not_matched = .~(sed .| cover .| volc .| plut .| metaign .| metased .| lowgrade .| highgrade)
        for i = 1:length(volctypes)
        volc[not_matched] .|= containsi.(rocktype[not_matched],volctypes[i])
        end
        for i = 1:length(pluttypes)
        plut[not_matched] .|= containsi.(rocktype[not_matched],pluttypes[i])
        end
        for i = 1:length(metaigntypes)
        metaign[not_matched] .|= containsi.(rocktype[not_matched],metaigntypes[i])
        end
        for i = 1:length(metasedtypes)
        metased[not_matched] .|= containsi.(rocktype[not_matched],metasedtypes[i])
        end
        for i = 1:length(lowgradetypes)
        lowgrade[not_matched] .|= containsi.(rocktype[not_matched],lowgradetypes[i])
        end
        for i = 1:length(highgradetypes)
        highgrade[not_matched] .|= containsi.(rocktype[not_matched],highgradetypes[i])
        end
        for i = 1:length(hypabyssaltypes)
        hypabyssal[not_matched] .|= containsi.(rocktype[not_matched],hypabyssaltypes[i])
        end

        # Then rockname
        not_matched = .~(sed .| ign .| met .| cover)
        for i = 1:length(sedtypes)
            sed[not_matched] .|= containsi.(rockname[not_matched],sedtypes[i])
        end
        for i = 1:length(igntypes)
            ign[not_matched] .|= containsi.(rockname[not_matched],igntypes[i])
        end
        for i = 1:length(mettypes)
            met[not_matched] .|= containsi.(rockname[not_matched],mettypes[i])
        end
            for i = 1:length(covertypes)
            cover[not_matched] .|= containsi.(rockname[not_matched],covertypes[i])
        end

        not_matched = .~(sed .| cover .| volc .| plut .| metaign .| metased .| lowgrade .| highgrade)
        for i = 1:length(volctypes)
            volc[not_matched] .|= containsi.(rockname[not_matched],volctypes[i])
        end
        for i = 1:length(pluttypes)
            plut[not_matched] .|= containsi.(rockname[not_matched],pluttypes[i])
        end
        for i = 1:length(metaigntypes)
            metaign[not_matched] .|= containsi.(rockname[not_matched],metaigntypes[i])
        end
        for i = 1:length(metasedtypes)
            metased[not_matched] .|= containsi.(rockname[not_matched],metasedtypes[i])
        end
        for i = 1:length(lowgradetypes)
            lowgrade[not_matched] .|= containsi.(rockname[not_matched],lowgradetypes[i])
        end
        for i = 1:length(highgradetypes)
            highgrade[not_matched] .|= containsi.(rockname[not_matched],highgradetypes[i])
        end
        for i = 1:length(hypabyssaltypes)
            hypabyssal[not_matched] .|= containsi.(rockname[not_matched],hypabyssaltypes[i])
        end

        # Then rockdescrip
        not_matched = .~(sed .| ign .| met .| cover)
        for i = 1:length(sedtypes)
            sed[not_matched] .|= sed[not_matched] .| containsi.(rockdescrip[not_matched],sedtypes[i])
        end
        for i = 1:length(igntypes)
            ign[not_matched] .|= containsi.(rockdescrip[not_matched],igntypes[i])
        end
        for i = 1:length(mettypes)
            met[not_matched] .|= containsi.(rockdescrip[not_matched],mettypes[i])
        end
        for i = 1:length(covertypes)
        cover[not_matched] .|= containsi.(rockdescrip[not_matched],covertypes[i])
        end

        not_matched = .~(sed .| cover .| volc .| plut .| metaign .| metased .| lowgrade .| highgrade)
        for i = 1:length(volctypes)
            volc[not_matched] .|= containsi.(rockdescrip[not_matched],volctypes[i])
        end
        for i = 1:length(pluttypes)
            plut[not_matched] .|= containsi.(rockdescrip[not_matched],pluttypes[i])
        end
        for i = 1:length(metaigntypes)
            metaign[not_matched] .|= containsi.(rockdescrip[not_matched],metaigntypes[i])
        end
        for i = 1:length(metasedtypes)
            metased[not_matched] .|= containsi.(rockdescrip[not_matched],metasedtypes[i])
        end
        for i = 1:length(lowgradetypes)
            lowgrade[not_matched] .|= containsi.(rockdescrip[not_matched],lowgradetypes[i])
        end
        for i = 1:length(highgradetypes)
            highgrade[not_matched] .|= containsi.(rockdescrip[not_matched],highgradetypes[i])
        end
        for i = 1:length(hypabyssaltypes)
            hypabyssal[not_matched] .|= containsi.(rockdescrip[not_matched],hypabyssaltypes[i])
        end

        return sed, ign, met, volc, plut, hypabyssal, metaign, metased, lowgrade, highgrade, cover
    end
end


## --- End of file