## --- Generate random points on the continental crust
"""
```julia
gen_continental_points(npoints, etopo)
```

Generate an `npoints` element-long lists of latitudes, longitudes, and elevations for points 
randomly and uniformly distributed across the continental crust. 

Requires `etopo` matrix of 1 arc minute resolution global relief model. Units are meters of 
elevation and decimal degrees of latitude and longitude.

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

Find erosion rate in mm/kyr given slope `slp`.
"""
function emmkyr(slp)
    return 10^(slp*0.00567517 + 0.971075)
end


## --- Functions for querying macrostrat
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


## --- Find the average value of slope over an area
"""
```julia
avg_over_area(data::Matrix, lat::Vector, lon::Vector, sf::Number=240; 
    halfwidth::Number=1, 
    maxpossible::Number=0xffff)
```
Find the average value of `data` over an area with radius `halfwidth` (units of gridcells) at 
coordinates `lat` and `lon`.

This is distinct from `StatGeochem`'s `aveslope`. This function finds the average over an 
area, not the average slope for a specific point. For example, this might be given a matrix 
of maximum slope at each point on Earth, and would return the _average maximum slope_ at 
each point of interest.

# Arguments

    sf::Number=240

Scale factor (cells per degree) for the SRTM15+ data. For 15 arc-second resolution, the scale
factor is 240, because 15 arc-seconds go into 1 arc-degree 60 * 4 = 240 times.

    maxpossible::Number=0xffff

The maximum possible value for the variable of interest; variables with values greater than 
this are ignored.

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

## --- Simplify statistics calculations
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
Return the `NamedTuple` of `BitVector`s catagorizing Macrostrat `rocktype`, 
`rockname`, and `rockdescrip` as sedimentary, igneous, metamorphic, or cover. 

# Keyword argument `major`
If `major` is `true`, return only major rock types:
```
sed ign met cover
```

If `major` is `false`, return both major rock types and more granular subcatagories. Note 
that major rock types include more granular subcatagories; i.e. `ign` includes all rock 
catagorized as `volc` and `plut`, as well as rocks that do not fall into either of those 
subcatagories. To match this catagorization to the EarthChem system, `plut` includes
hypabyssal rocks, and `chert` includes banded iron formations.
```
siliciclast carb chert evaporite coal sed volc plut ign metased metaign cataclastic met cover
```

# Example
```julia-repl
cats = match_rocktype(rocktype, rockname, rockdescrip, major=true)
NamedTuple with 4 elements:
  sed    = BitVector(50000,)    [true ... true]
  ign    = BitVector(50000,)    [false ... false]
  met    = BitVector(50000,)    [false ... false]
  cover  = BitVector(50000,)    [false ... false]
```
"""
function match_rocktype(rocktype, rockname, rockdescrip; major=false)
    npoints = length(rocktype)

    # Rock types based on EarthChem system in RockNameInference.m
    # Sedimentary
    siliciclasttypes = ["siliciclast", "conglomerat", "sand", "psamm", "arenit", "arkos", "silt",]
    shaletypes = ["mud", "clay","shale", "wacke", "argillite", "argillaceous", "flysch", "pelit", 
        "turbidite",]
    carbtypes = ["carbonate", "limestone", "dolo", "marl", "chalk", "travertine", "tavertine", 
        "teravertine", "tufa",]
    cherttypes = ["chert", "banded iron",]
    evaporitetypes = ["evaporite", "gypsum", "salt flat",]
    coaltypes = ["coal", "anthracite",]
    
    sedtypes = vcat(["sediment", "fluv", "clast", "gravel", "pebble", "caliche", "boulder", 
        "diamict", "tillite", "stream", "beach", "terrace",  "marine deposits",  "paleosol"],
        siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, coaltypes)

    # Igneous
    volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme",
        "pipe", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", "tephra", 
        "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", 
        "fonolito", "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", 
        "komatiite",]
    hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", 
        "dolerit", "diabase", "porphyry", "microgranite"]
    pluttypes = vcat(["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", 
        "monzonit", "syenit", "peridot", "dunit", "harzburg", "anorthosite", "mangerite", 
        "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", 
        "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",], 
        hypabyssaltypes)
    
    igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",], volctypes, pluttypes)

    # Metamorphic
    metasedtypes = ["para", "metased", "meta-sed", "quartzite", "marble", "slate", "phyllite",]
    metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone",
        "eclogite", "metabasite",]
    lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", 
        "gossan", "alter", "hydrothermal", "palagonite",]
    highgradetypes = ["crystalline", "basement", "marble", "skarn", "schist", "blueschist", 
        "gneiss", "amphibolit", "eclogite", "granulit", "hornfels", "granofels", "sanidinite", 
        "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
        "harzburg", "high grade metamorphic"]
    cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]
    
    mettypes = vcat(["meta", "calc silicate",], metasedtypes, metaigntypes, lowgradetypes, 
        highgradetypes, cataclastictypes)

    # Cover
    covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", 
        "laterite", "surficial deposits", "talus", "scree", "mass-wasting", "slide", 
        "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", 
        "debris"]
    
    # Define what we're looking for
    if major
        # Preallocate
        cats = (
            sed = falses(npoints),
            ign = falses(npoints),
            met = falses(npoints),
            cover = falses(npoints)
        )
        typelist = [sedtypes, igntypes, mettypes, covertypes]
    else
        cats = (
            siliciclast = falses(npoints),
            shale = falses(npoints),
            carb = falses(npoints),
            chert = falses(npoints),
            evaporite = falses(npoints),
            coal = falses(npoints),
            sed = falses(npoints),

            volc = falses(npoints),
            plut = falses(npoints),
            ign = falses(npoints),

            metased = falses(npoints),
            metaign = falses(npoints),
            met = falses(npoints),

            cover = falses(npoints),
        )
        typelist = [siliciclasttypes, shaletypes, carbtypes, cherttypes, evaporitetypes, 
            coaltypes, sedtypes, volctypes, pluttypes, igntypes, metasedtypes, metaigntypes, 
            mettypes, covertypes
        ]
    end

    # Check major lithology first
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j] .|= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], typelist[j][i]) : false)
        end
    end

    println("finished major lithology")

    # Check the rest of rocktype
    not_matched = .~(cats.sed .| cats.ign .| cats.met .| cats.cover)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j][not_matched].|= containsi.(rocktype[not_matched], typelist[j][i])
        end
    end

    println("finished rest of rocktype")

    # Then rockname
    not_matched = .~(cats.sed .| cats.ign .| cats.met .| cats.cover)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j][not_matched].|= containsi.(rockname[not_matched], typelist[j][i])
        end
    end

    println("finished rockname")

    # Then rockdescrip
    not_matched = .~(cats.sed .| cats.ign .| cats.met .| cats.cover)
    for j in eachindex(typelist)
        for i = 1:length(typelist[j])
            cats[j][not_matched].|= containsi.(rockdescrip[not_matched], typelist[j][i])
        end
    end

    println("finished rock descrip")

    return cats
end


"""
```julia
match_earthchem(type; major=false)
```
Classify EarthChem bulk.mat `type` codes to rock types. Returns `BitVector`s. 
    
If  `major` is `true`, returns:
```
sed ign met
```
Metasedimentary and metaigneous rocks are respectively classified as sedimentary and igneous
rocks.

If  `major` is `false`, returns:
```
sed ign volc plut met metased metaign
```
Sedimentary rocks are not further catagorized. Igneous rocks are catagorized as volcanic (`volc`), 
plutonic (`plut`), _or_ undifferentiated igneous (`ign`). Metamorphic rocks are classified as 
metasedimentary (`metased`), metaigenous (`metaign`), _or_ undifferentiated metamorphic (`met`).

"""
function match_earthchem(type; major=false)
    npoints = length(type)

    # Preallocate
    cats = (
        alluvium = falses(npoints),
        siliciclast = falses(npoints),
        shale = falses(npoints),
        carb = falses(npoints),
        chert = falses(npoints),
        evaporite = falses(npoints),
        phosphorite = falses(npoints),
        coal = falses(npoints),
        volcaniclast = falses(npoints),
        sed = falses(npoints),

        volc = falses(npoints),
        plut = falses(npoints),
        ign = falses(npoints),

        metased = falses(npoints),
        metaign = falses(npoints),
        met = falses(npoints),
    )

    # Codes for types *in the order they appear in cats*
    codes = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.0, 
        3.1, 3.2, 3.0,
        2.1, 2.3, 2.0, 
    ]

    # Find all types
    for i in eachindex(codes)
        cats[i][findall(==(codes[i]), type)] .= true
    end
    
    # Assign all subtypes to parent type as appropriate. Alluvium is roughly equivilent to cover
    if major
    # To match Macrostrat, metaseds are seds, metaigns are igns
        cats.sed .= cats.sed .| cats.siliciclast .| cats.shale .| cats.carb .| cats.chert .| 
            cats.evaporite .| cats.phosphorite .| cats.coal .| cats.volcaniclast .| cats.metased
        cats.ign .= cats.ign .| cats.volc .| cats.plut .| cats.metaign

        majorcats = (sed=cats.sed, ign=cats.ign, met=cats.met)
        return majorcats
    else
        cats.sed .= cats.sed .| cats.siliciclast .| cats.shale .| cats.carb .| cats.chert .| 
            cats.evaporite .| cats.phosphorite .| cats.coal .| cats.volcaniclast
        cats.ign .= cats.ign .| cats.volc .| cats.plut
        cats.met .= cats.met .| cats.metased .| cats.metaign

        return cats
    end

end


## --- Functions for matching samples
    """
    ```julia
    major_elements(bulk, bulk_filter)
    ```
    Compute mean and standard deviation of major elements in `bulk` for each rock type.

    Major elements are defined based on Faye and Ødegård 1975 
    (https://www.ngu.no/filearchive/NGUPublikasjoner/NGUnr_322_Bulletin_35_Faye_35_53.pdf)
    """
    function major_elements(bulk, bulk_filter)
        vals = Array{NamedTuple}(undef, length(bulk_filter), 1)
        for i = 1:length(bulk_filter)
            elem = (
                SiO2 = (m = nanmean(bulk.SiO2[bulk_filter[i]]), e = nanstd(bulk.SiO2[bulk_filter[i]])),
                Al2O3 = (m = nanmean(bulk.Al2O3[bulk_filter[i]]), e = nanstd(bulk.Al2O3[bulk_filter[i]])),
                Fe2O3T = (m = nanmean(bulk.Fe2O3T[bulk_filter[i]]), e = nanstd(bulk.Fe2O3T[bulk_filter[i]])),
                TiO2 = (m = nanmean(bulk.TiO2[bulk_filter[i]]), e = nanstd(bulk.TiO2[bulk_filter[i]])),
                MgO = (m = nanmean(bulk.MgO[bulk_filter[i]]), e = nanstd(bulk.MgO[bulk_filter[i]])),
                CaO = (m = nanmean(bulk.CaO[bulk_filter[i]]), e = nanstd(bulk.CaO[bulk_filter[i]])),
                Na2O = (m = nanmean(bulk.Na2O[bulk_filter[i]]), e = nanstd(bulk.Na2O[bulk_filter[i]])),
                K2O = (m = nanmean(bulk.K2O[bulk_filter[i]]), e = nanstd(bulk.K2O[bulk_filter[i]])),
                MnO = (m = nanmean(bulk.MnO[bulk_filter[i]]), e = nanstd(bulk.MnO[bulk_filter[i]]))
            )
            vals[i] = elem
        end

        geochemkeys = (
            :alluvium, :siliciclast, :shale, :carb, :chert, :evaporite, :phosphorite, :coal, 
                :volcaniclast, :sed,
            :volc, :plut, :ign,
            :metased, :metaign, :met
        )

        # This is actually approx. 15ms slower than manually assigning things element by element...
        return NamedTuple{geochemkeys}(vals)
    end

    function major_elements(bulk, bulksamples::BitVector, bulk_filter::BitVector)
        elem = (
            SiO2 = (m = nanmean(bulk.SiO2[bulksamples][bulk_filter]), 
                e = nanstd(bulk.SiO2[bulksamples][bulk_filter])),
            Al2O3 = (m = nanmean(bulk.Al2O3[bulksamples][bulk_filter]), 
                e = nanstd(bulk.Al2O3[bulksamples][bulk_filter])),
            Fe2O3T = (m = nanmean(bulk.Fe2O3T[bulksamples][bulk_filter]), 
                e = nanstd(bulk.Fe2O3T[bulksamples][bulk_filter])),
            TiO2 = (m = nanmean(bulk.TiO2[bulksamples][bulk_filter]), 
                e = nanstd(bulk.TiO2[bulksamples][bulk_filter])),
            MgO = (m = nanmean(bulk.MgO[bulksamples][bulk_filter]), 
                e = nanstd(bulk.MgO[bulksamples][bulk_filter])),
            CaO = (m = nanmean(bulk.CaO[bulksamples][bulk_filter]), 
                e = nanstd(bulk.CaO[bulksamples][bulk_filter])),
            Na2O = (m = nanmean(bulk.Na2O[bulksamples][bulk_filter]), 
                e = nanstd(bulk.Na2O[bulksamples][bulk_filter])),
            K2O = (m = nanmean(bulk.K2O[bulksamples][bulk_filter]), 
                e = nanstd(bulk.K2O[bulksamples][bulk_filter])),
            MnO = (m = nanmean(bulk.MnO[bulksamples][bulk_filter]), 
                e = nanstd(bulk.MnO[bulksamples][bulk_filter]))
        )

        return elem
    end


## --- Find all EarthChem samples matching a Macrostrat sample
    """
    ```julia
    find_earthchem(rocktype::String, rockname::String, rockdescrip::String, 
        bulksamples::BitVector, bulk_rockname
    )
    ```

    For a single macrostrat sample, find all possible matching EarthChem samples.
    """
    function find_earthchem(rocktype::String, rockname::String, rockdescrip::String, bulkrockname, bulkrocktype, bulkmaterial)
        # Check for matches in rock name
        rocknamelist = unique(bulkrockname)
        found = check_matches(rocktype, rockname, rockdescrip, rocknamelist)

        if count(found)!=0
        # If samples were found, get all EarthChem samples that match the rock name
            earthchem_matches = falses(length(bulkrockname))
            to_find = rocknamelist[found]
            t = findall(!=(""), to_find)
            to_find = to_find[t]

            for i in eachindex(to_find)
                earthchem_matches .|= containsi.(bulkrockname, to_find[i])
            end

            return earthchem_matches
        
        else
        # If no samples were found, check material
            rocknamelist = unique(bulkmaterial)
            found = check_matches(rocktype, rockname, rockdescrip, rocknamelist)

            if count(found)!=0
            # Process and return
                earthchem_matches = falses(length(bulkmaterial))
                to_find = rocknamelist[found]
                t = findall(!=(""), to_find)
                to_find = to_find[t]

                for i in eachindex(to_find)
                    earthchem_matches .|= containsi.(bulkmaterial, to_find[i])
                end

                return earthchem_matches

            else
            # If still none, check type
                rocknamelist = unique(bulkrocktype)
                found = check_matches(rocktype, rockname, rockdescrip, rocknamelist)

                if count(found)!=0
                    earthchem_matches = falses(length(bulkrocktype))
                    to_find = rocknamelist[found]
                    t = findall(!=(""), to_find)
                    to_find = to_find[t]

                    for i in eachindex(to_find)
                        earthchem_matches .|= containsi.(bulkrocktype, to_find[i])
                    end

                    return earthchem_matches

                else
                    @warn "No matching EarthChem samples found for $rocktype, $rockname, $rockdescrip"
                    return falses(length(bulkrocktype))
                end
            end
            
        end
    end


## --- Check for matches
    function check_matches(rocktype::String, rockname::String, rockdescrip::String, rocknamelist)
        found = falses(length(rocknamelist))

        # Get list of matches from the sample. Check major lithology first
        for i in eachindex(rocknamelist)
            found[i] |= (match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], rocknamelist[i]) : false)
        end

        # If no matches, try the rest of rocktype
        if count(found)==0
            for i in eachindex(rocknamelist)
                found[i] |= containsi.(rocktype, rocknamelist[i])
            end

            # If still no matches, try rockname
            if count(found)==0
                for i in eachindex(rocknamelist)
                    found[i] |= containsi.(rockname, rocknamelist[i])
                end

                # Then rock descrip...
                if count(found)==0
                    for i in eachindex(rocknamelist)
                        found[i] |= containsi.(rockdescrip, rocknamelist[i])
                    end
                end
            end
        end

        return found
    end

## --- Find likelihoods
    function likelihood(bulkage::Vector, sampleage::Number,
            bulklat::Vector, bulklon::Vector, samplelat::Number, samplelon::Number,
            bulkgeochem::NamedTuple, samplegeochem::NamedTuple
        )
        # Preallocate
        npoints = length(bulkage)
        ll_age = Array{Float64}(undef, npoints, 1)
        ll_dist = Array{Float64}(undef, npoints, 1)
        ll_total = Array{Float64}(undef, npoints, 1)

        # Replace missing values: this will penalize but not exclude missing data
        @inbounds for i in 1:npoints
            if isnan(bulkage[i])
                bulkage[i] = -sampleage
            end

            # Assume if one coordinate is missing, so is the other one
            if isnan(bulklat[i])
                bulklat[i] = -samplelat
                bulklon[i] = -samplelon
            end
        end

        @turbo for i in 1:npoints
            # Age (σ = 38 Ma)
            ll_age[i] = -((bulkage[i] - sampleage)^2)/(38^2)

            # Distance (σ = 1.8 arc degrees)
            ll_dist[i] = -((haversine(samplelat, samplelon, bulklat[i], bulklon[i]))^2)/(1.8^2)
        end

        @. ll_total = ll_age + ll_dist
        
        # Geochemical log-likelihoods
        for elem in eachindex(bulkgeochem)
            @turbo for i in 1:npoints
                ll_total[i] += -((bulkgeochem[elem][i] - samplegeochem[elem].m)^2)/(samplegeochem[elem].e^2)
            end
        end

        matched_sample = rand_prop_liklihood(ll_total)
        return matched_sample
    end

## --- Randomly select a sample, weighted based on likelihood
    function rand_prop_liklihood(ll)
        sum_likelihoods = sum(exp, ll)
        r = rand()*sum_likelihoods
        s = zero(typeof(sum_likelihoods))
        @inbounds for i in eachindex(ll)
            s += exp(ll[i])
            if s > r
                return i
            end
        end
        return lastindex(ll)
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




## --- End of file
