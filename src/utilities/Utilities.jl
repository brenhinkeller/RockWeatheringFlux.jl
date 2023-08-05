## --- Run all sub-utilities
    include("Definitions.jl")
    include("Slope.jl")
    include("NaNMeasurements.jl")
    include("Macrostrat.jl")

    
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
    emmkyr(slp) = exp(slp * (0.0091 ± 0.0095) + (3.1 ± 1.9))

    # Previously 10^(slp*0.00567517 + 0.971075)


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

    # Defaults and Other Keyword Arguments

    - `sf::Number=240`: Scale factor (cells per degree) for the SRTM15+ data. For 15 arc-second 
    resolution, the scale factor is 240, because 15 arc-seconds go into 1 arc-degree 60 * 4 = 
    240 times.

    - `maxpossible::Number=0xffff`: The maximum possible value for the variable of interest; 
    variables with values greater than this are ignored.

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
    match_rocktype(rocktype, rockname, rockdescrip; major=false)
    match_rocktype(writtentype::AbstractArray{String})
    ```

    Return the `NamedTuple` of `BitVector`s catagorizing Macrostrat `rocktype`, 
    `rockname`, and `rockdescrip` as sedimentary, igneous, metamorphic, or cover. 
    Alternatively, return types already stored as strings in `writtentype`.

    ### Optional keyword argument `major`
    Note that this argument is only valid when parsing data directly from Macrostrat responses.

    `true` returns: `sed, ign, met`

    `false` returns: `siliciclast, shale, carb, chert, evaporite, coal, sed, volc, plut, 
    ign, metased, metaign, met, cover`

    Note that major rock types include more granular subcategories; i.e. `ign` includes 
    all rock catagorized as `volc` and `plut`, as well as rocks that do not fall into 
    either of those subcategories. To match the catagorization to the EarthChem system, 
    `plut` includes hypabyssal rocks, and `chert` includes banded iron formations.

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
    function match_rocktype(rocktype::AbstractArray, rockname::AbstractArray, 
        rockdescrip::AbstractArray; major=false)

        # Get rock type classifications and initialized BitVector
        typelist, cats = get_rock_class(major, length(rocktype))

        # If you can't improve the algorithm the least you can do is add a progress bar
        p = Progress(length(typelist)*4+1, desc="Finding Macrostrat rock types...")
        next!(p)

        # Check major lithology first
        for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    cats[j][k] = match(r"major.*?{(.*?)}", rocktype[k]) |> x -> isa(x,RegexMatch) ? containsi(x[1], typelist[j][i]) : false
                end
            end
            next!(p)
        end

        # Check the rest of rocktype
        not_matched = find_unmatched(cats, major=major)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rocktype[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # Then rockname
        not_matched = find_unmatched(cats, major=major)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rockname[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # Then rockdescrip
        not_matched = find_unmatched(cats, major=major)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rockdescrip[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        return un_multimatch!(cats, major)
    end

    function match_rocktype(writtentype::AbstractArray{String})
        typelist, cats = get_rock_class(false, length(writtentype))

        # Get all of the written types into their proper place
        for i in eachindex(writtentype)
            try
                cats[Symbol(writtentype[i])][i] = true
            catch
                continue
            end
        end

        # Make sure that sed / ign / met are true when their subtypes are true
        minortypes = (
            sed = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal),
            ign = (:volc, :plut),
            met = (:metased, :metaign)
        )
        for k in keys(minortypes)
            for i in minortypes[k]
                cats[k][cats[i]] .= true
            end
        end

        return cats
    end


    """
    ```julia
    find_unmatched(cats; [major])
    ```

    Given a `Tuple` of `BitVectors`, return a `BitVector` that is `true` at index `i` iff 
    all elements of the `Tuple` are `false` at index `i`.

    If `cats` is a `NamedTuple` of rock types defined by `get_rock_class`, specify `major` 
    as `true` or `false` to decrease runtime. `major` is `true` if `cats` contains only 
    `sed`, `ign`, `met`, and `cover`.

    # Example
    ```
    julia> find_unmatched(cats)
    500-element BitVector:
    0
    1
    1
    1
    ⋮
    1
    1
    1
    1
    ```
    """
    function find_unmatched(cats; major::Bool)
        if major
            return .!(cats.sed .| cats.ign .| cats.met .| cats.cover)
        else
            cats.sed .= cats.sed .| cats.siliciclast .| cats.shale .| cats.carb .| cats.chert .| 
                cats.evaporite .| cats.coal
            cats.ign .= cats.ign .| cats.volc .| cats.plut
            cats.met .= cats.met .| cats.metased .| cats.metaign

            return .!(cats.sed .| cats.ign .| cats.met .| cats.cover)
        end
    end

    function find_unmatched(cats)
        matched = falses(length(cats[1]))
        for i in eachindex(cats)
            matched .|= cats[i]
        end

        return .!matched
    end


## --- Return a list of matching rock names

    """
    ```
    match_rockname(rocktype::AbstractArray, rockname::AbstractArray, 
        rockdescrip::AbstractArray)
    ```

    Find samples in Macrostrat matching pre-defined rock names (e.g., basalt, sandstone,
    gneiss, etc.). Return a `NamedTuple` of `BitVectors` where the elements of the 
    `NamedTuple` are the rock names, and the `BitVectors` are true at element `i` iff
    sample `i` is that rock name.

    Does not match rock names defined as cover by `get_rock_class`.

    See also: `match_rocktype`.
    """
    function match_rockname(rocktype::AbstractArray, rockname::AbstractArray, 
        rockdescrip::AbstractArray)

        # Get rock names as one long list, excluding cover
        typelist = get_rock_class(true, length(rocktype))[1]
        typelist = unique((typelist.sed..., typelist.met..., typelist.ign...))

        # Initialize a NamedTuple with a BitVector for each rock name
        cats = NamedTuple{Symbol.(Tuple(typelist))}([falses(length(rocktype))
            for _ in 1:length(typelist)]
        )

        # If you can't improve the algorithm the least you can do is add a progress bar
        p = Progress(length(typelist)*4+1, desc="Finding Macrostrat rock names...")
        next!(p)

        # Check major lithology first
        for i in eachindex(typelist)
            for j in eachindex(rocktype)
                cats[i] .= match(r"major.*?{(.*?)}", rocktype[j]) |> x -> isa(x,RegexMatch) ? 
                containsi(x[1], typelist[i]) : false
            end
            next!(p)
        end

        # Check the rest of rocktype
        not_matched = find_unmatched(cats)
        @inbounds for i in eachindex(typelist)
            for j in eachindex(rocktype)
                !not_matched[j] && continue
                cats[i][j] |= containsi(rocktype[j], typelist[i])
            end
            next!(p)
        end

        # Then rockname
        not_matched = find_unmatched(cats)
        @inbounds for i in eachindex(typelist)
            for j in eachindex(rockname)
                !not_matched[j] && continue
                cats[i][j] |= containsi(rockname[j], typelist[i])
            end
            next!(p)
        end

        # Then rockdescrip
        not_matched = find_unmatched(cats)
        @inbounds for i in eachindex(typelist)
            for j in eachindex(rockdescrip)
                !not_matched[j] && continue
                cats[i][j] |= containsi(rockdescrip[j], typelist[i])
            end
            next!(p)
        end

        return cats
    end


## --- More Stuff. Maybe one day this file will be organized

    """
    ```julia
    get_type(cats::NamedTuple, i; [all_keys=false])
    ```

    Return the first key of `cats` where `i` is `true`. Optionally specify `allkeys`=`true`
    to return _all_ keys where `i` is true. Assumes equal length elements of `cats`.

    # Examples
    ```julia-repl
    julia> get_type(macro_cats, 254)
    :sed

    julia> get_type(name_cats, 3, all_keys=true)
    (:gravel, :sand, :silt, :clay)
    ```
    """
    get_type(cats::NamedTuple, i::Int64; all_keys::Bool=false) = _get_type(cats, i, static(all_keys))

    function _get_type(cats, i, all_keys::False)
        @assert 0 < i <= length(cats[1]) "Index $i out of bounds."

        @inbounds for k in keys(cats)
            cats[k][i] && return Symbol(k)
        end

        return nothing
    end

    function _get_type(cats, i, allkeys::True)
        @assert 0 < i <= length(cats[1]) "Index $i out of bounds."

        catkeys = keys(cats)
        keymatches = falses(length(catkeys))

        @inbounds for k in eachindex(catkeys)
            cats[k][i] && (keymatches[k]=true)
        end

        count(keymatches)==0 && return nothing
        return catkeys[keymatches]
    end

    """
    ```julia
    class_up(typelist, name::String)
    ```

    Find the rock type category for `name`. The `name` must be included in the rock name
    `typelist` returned by `get_rock_class`.

    ### Example
    ```julia-repl
    julia> class_up(typelist, "dacit")
    :volc
    ```
    """
    function class_up(typelist, name::String)
        for k in keys(typelist)
            for i in typelist[k]
                name == i && return k
            end
        end

        return nothing
    end

    """
    ```julia
    match_earthchem(type; major=false)
    ```
    Classify EarthChem bulk.mat `type` codes to rock types. Returns a `NamedTuple` of `BitVector`s. 

    # Optional Keyword Argument `major`
    `true` returns 

        sed, ign, met

    `false` returns

        alluvium, siliciclast, shale, carb, chert, evaporite, phosphorite, coal, volcaniclast, sed, volc, plut, ign, metased, metaign, met

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


    """
    ```julia
    major_elements(bulk::NamedTuple, [bulk_filter::BitVector])
    ```
    Compute mean and standard deviation of major elements (defined in `get_elements`) in 
    `bulk`. Optionally restrict samples by `bulk_filter`.
    """
    function major_elements(bulk::NamedTuple, bulk_filter::BitVector=trues(length(bulk[1])))
        major, = get_elements()

        element = [NamedTuple{(:m, :e)}(tuple.(
            nanmean(bulk[i][bulk_filter]),
            nanstd(bulk[i][bulk_filter])))
        for i in major]

        return NamedTuple{Tuple(major)}(element)
    end


## --- Find all EarthChem samples matching a Macrostrat sample
    """
    ```julia
    find_earthchem(name::String, rockname::AbstractArray, rocktype::AbstractArray, 
        rockmaterial::AbstractArray)
    ```

    Given a rock `name` (e.g. "basalt"), find all matching EarthChem samples from `rockname`,
    `rocktype`, and `rockmaterial`, where:
       * `rockname`s are rock names (e.g., sandstone, shale, chert)
       * `rocktype`s are "specialized" rock classes (e.g., volcanic, plutonic, breccia)
       * `rockmaterial`s are "generalized" rock classes (e.g., igneous, sedimentary, 
           metamorphic)
    Returns a `BitVector`.

    General rock names should return general results: i.e., "igneous" should return _all_
    igneous rocks. Therefore, EarthChem fields are searched from most to least general.
    All EarthChem fields are searched for all samples (unlike `match_rocktype` and 
    `match_rockname`, which stop searching when a sample is matched).

    """
    function find_earthchem(name::String, rockname::AbstractArray, rocktype::AbstractArray, 
        rockmaterial::AbstractArray)

        # Preallocate
        matches = falses(length(rockname))

        # Material, type, name for each sample
        for i in eachindex(rockmaterial)
            matches[i] = containsi(rockmaterial[i], name)
            matches[i] |= containsi(rocktype[i], name)
            matches[i] |= containsi(rockname[i], name)
        end

        return matches
    end


## --- Correct units to wt.%
    """
    ```julia
    standardize_units!(bulkdata::AbstractArray{Float64}, bulkunits::AbstractArray{String}, 
        density::NamedTuple;
        elem::Symbol
    )
    ```

    Modify data `bulkdata` with units `bulkunits` so all data is in units of wt.%. 
    
    Samples without listed units are assumed to be reported with the most commonly used 
    known unit for that data set. Requires the density of some element oxides where units 
    are sometimes stored as vol.%. The key `elem` is used to identify the correct element
    density.
    """
    function standardize_units!(bulkdata::AbstractArray{Float64}, bulkunits::AbstractArray{String}, 
        density::NamedTuple;
        elem::Symbol
    )
        # Preallocate
        unitsused = Dict(zip(unique(bulkunits), zeros(Int64, length(unique(bulkunits)))))
        missinginfo = falses(length(bulkdata))

        # If vol% is used to store units, ensure information is present to correct to wt.%
        if containsi(keys(unitsused), "VOL%")
            @assert containsi(keys(density), elem) "No density information present for $(string(elem))"
            density_i = density[elem]
        end

        # Correct units for data where unit is known
        for i in eachindex(bulkunits)
            # Record unit used
            unitsused[bulkunits[i]] += 1

            # If data exists, correct to wt.%
            if bulkunits[i]=="WT%" || bulkunits[i]=="NA"
                continue
            elseif bulkunits[i]=="DPM/G"
                bulkdata[i] = NaN
                unitsused[bulkunits[i]] -= 1
            elseif bulkunits[i]=="PPM" || bulkunits[i]=="MICROGRAM PER GRAM" || bulkunits[i]=="MILLIGRAM PER KILOGRAM"
                bulkdata[i] /= 10000
            elseif bulkunits[i]=="VOL%"
                bulkdata[i] *= density_i
            else
                missinginfo[i] = true
            end
        end

        # Find the most commonly used unit
        delete!(unitsused, "")
        assumedunit = maximum(unitsused)[1]

        # Correct units for data with missing unit information
        for i in eachindex(missinginfo)
            if missinginfo[i]
                if assumedunit=="WT%" || assumedunit=="NA"
                    break
                elseif assumedunit=="PPM" || assumedunit=="MICROGRAM PER GRAM" || assumedunit=="MILLIGRAM PER KILOGRAM"
                    bulkdata[i] /= 10000
                elseif assumedunit=="VOL%"
                    bulkdata[i] *= density_i
                end
            end
        end
    end


## --- Find likelihoods
    """
    ```julia
    likelihood(bulkage::AbstractArray, sampleage::Number,
        bulklat::AbstractArray, bulklon::AbstractArray, samplelat::Number, 
        samplelon::Number, bulkgeochem::NamedTuple, samplegeochem::NamedTuple,
        sampleidx::AbstractArray)
    ````

    For a Macrostrat sample with age `sampleage`, location (`samplelat`, `samplelon`), and
    estimated geochemistry `samplegeochem`, find the index of the EarthChem sample that 
    is most likely to represent the Macrostrat sample.
    """
    function likelihood(bulkage::AbstractArray, sampleage::Number,
            bulklat::AbstractArray, bulklon::AbstractArray, samplelat::Number, 
            samplelon::Number, bulkgeochem::NamedTuple, samplegeochem::NamedTuple,
            sampleidx::AbstractArray)

        # Preallocate
        npoints = length(bulkage)
        ll_age = Array{Float64}(undef, npoints, 1)
        ll_dist = Array{Float64}(undef, npoints, 1)
        ll_chem = zeros(npoints)

        # Replace missing values: this will penalize but not exclude missing data
        @inbounds for i in 1:npoints
            if isnan(bulkage[i])
                bulkage[i] = ifelse(sampleage < 1900, 3800, 0)
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

        # Geochemical log-likelihoods
        for name in eachindex(samplegeochem)
            for elem in eachindex(bulkgeochem)
                @turbo for i in 1:npoints
                    ll_chem[i] += -(((bulkgeochem[elem][i] - 
                        samplegeochem[name][elem].m)^2)/(samplegeochem[name][elem].e^2)
                    )
                end
            end
        end
        ll_chem ./= length(samplegeochem)

        matched_sample = rand_prop_liklihood(ll_age .+ ll_dist .+ ll_chem)
        return sampleidx[matched_sample]
    end


## --- Randomly select a sample, weighted based on likelihood
    """
    ```julia
    rand_prop_liklihood(ll)
    ```

    Weighted-random selection of an index based on log-likelihoods `ll`.
    """
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


## --- Calculate wt.% and flux by rock type
    # """
    # ```julia
    # function flux_source(bulk::AbstractArray, bulkidx::Vector{Int64}, erosion::NamedTuple, 
    #     macro_cats::NamedTuple, crustal_area::NamedTuple; 
    #     unitcodes::AbstractMatrix, unitdecoder::AbstractMatrix, crustal_density::Number=2750, 
    #     elem::String="")
    # ```

    # For a specified element in `bulk`, calculate the average wt.% and flux (kg/yr) by rocktype. 
    # Calculate the total global flux of that element (kg/yr). Return the number of samples `n`.

    # Note that `erosion`, `macro_cats`, and `crustal_area` _must_ contain at minimum the keys:
    # ```
    # :siliciclast, :shale, :carb, :chert, :evaporite, :coal, :sed, :volc, :plut, :ign, :metased, 
    # :metaign, :met
    # ```
    # Keys must be type `Symbol`.

    # ### Optional Keyword Arguments
    # - `crustal_density::Number=2750`: Average crustal density (kg/m³).
    # - `elem::String=""`: Element being analyzed, for terminal printout.

    # # Example
    # ```julia-repl
    # julia> wt, flux, global_flux, n = flux_source(bulk.P2O5, bulkidx, erosion, macro_cats, crustal_area, elem="phosphorus")
    # [ Info: 44307 of 50000 phosphorus samples (1%) are not NaN
    # ```
    # """
    # function flux_source(bulk::AbstractArray, bulkidx::Vector{Int64}, erosion::NamedTuple, 
    #     macro_cats::NamedTuple, crustal_area::NamedTuple; 
    #     crustal_density::Number=2750, elem::String="", printinfo=false)

    #     # Preallocate
    #     allkeys = collect(keys(macro_cats))
    #     deleteat!(allkeys, findall(x->x==:cover,allkeys))       # Do not compute cover

    #     allinitvals = fill(NaN ± NaN, length(allkeys))
    #     npoints = length(bulkidx)

    #     wt = Dict(zip(allkeys, allinitvals))
    #     flux = Dict(zip(allkeys, allinitvals))
    #     bulkdata = Array{Float64}(undef, npoints, 1)

    #     # Get EarthChem samples, if present
    #     for i in eachindex(bulkidx)
    #         (bulkidx[i] != 0) ? (bulkdata[i] = bulk[bulkidx[i]]) : (bulkdata[i] = NaN)
    #     end

    #     # Find how many samples have data for the element of interest
    #     n = length(findall(!isnan, bulkdata))
    #     if printinfo
    #         @info "$n of $npoints $elem samples ($(round(n/npoints*100, sigdigits=3))%) are not NaN"
    #     end

    #     # Calculate average wt.% for each rock type
    #     # TO DO: Maybe set no data to 0 instead of NaN? Would require re-writing a bit...
    #     for i in keys(wt)
    #         wt[i] = nanmean(bulkdata[macro_cats[i]]) ± nanstd(bulkdata[macro_cats[i]])
    #     end
    #     wt = NamedTuple{Tuple(keys(wt))}(values(wt))

    #     # Calculate provenance by rock type
    #     for i in keys(flux)
    #         flux[i] = erosion[i] * crustal_area[i] * wt[i] * crustal_density* 1e-8
    #     end
    #     flux = NamedTuple{Tuple(keys(flux))}(values(flux))

    #     # Compute global flux
    #     global_flux = nansum([flux.sed, flux.ign, flux.met])

    #     return wt, flux, global_flux, n
    # end
    

## --- Functions for dealing with measurements
    """
    ```julia
        unmeasurementify(A::AbstractArray{Measurement{Float64}})
        unmeasurementify(A::NamedTuple)
    ```
    
    Separate an Array or NamedTuple `A` of `measurements` into an array of values and an 
    array of errors.
    
    ## Example
    ```julia
    julia> A = [1 ± 0.1, 2 ± 0.1]
    2-element Vector{Measurement{Float64}}:
     1.0 ± 0.1
     2.0 ± 0.1
    
    julia> val, err = unmeasurementify(A)
    ([1.0, 2.0], [0.1, 0.1])
    ```
    """
    function unmeasurementify(A::AbstractArray{Measurement{Float64}})
        val = fill(NaN, length(A))
        err = fill(NaN, length(A))
        for i in eachindex(A)
            val[i] = A[i].val
            err[i] = A[i].err
        end
        return val, err
    end
    function unmeasurementify(A::NamedTuple)
        val = fill(NaN, length(A))
        err = fill(NaN, length(A))
        for i in eachindex(keys(A))
            val[i] = A[i].val
            err[i] = A[i].err
        end
        return val, err
    end

    
## --- Normalize compositions to 100%
    """
    ```julia
    normalize!(A::AbstractVector)
    ```

    Normalize the percentage values in `A` to 100%.
    """
    function normalize!(A::AbstractVector)
        sum_a = nansum(A)
        @turbo for i in eachindex(A)
            A[i] = A[i] / sum_a * 100
        end
        return A
    end


## --- Get points inside a polygon
    """
    ```julia
    points_in_shape(poly_x::AbstractArray, poly_y::AbstractArray,  
                    data_x::AbstractArray, data_y::AbstractArray) 
    ```

    Find all points from a set that lie within a polygon defined by the coordinates 
    (`poly_y`, `poly_x`). Points on an edge are included; returns a `BitVector` to filter
    points by position. Shapes are assumed to be self closing.

    Uses the even-odd rule to determine if points are inside the polygon
    (https://en.wikipedia.org/wiki/Even%E2%80%93odd_rule).

    To find latitude and longitude coordinates within a polygon, see `coords_in_shape`

    ### Example

    ```julia
    x, y, inpolygon = points_in_shape(poly_x, poly_y, data_x, data_y)
    ```
    """
    function points_in_shape(poly_x::AbstractArray{<:Number}, poly_y::AbstractArray{<:Number},  
        data_x::AbstractArray{<:Number}, data_y::AbstractArray{<:Number})

        # Define BitVectors for points in the polygon, and points that we know are in or outside
        inpolygon, known = falses(length(data_x)), falses(length(data_x))

        # Get the min / max of the polygon
        box_y = (min = minimum(poly_y), max = maximum(poly_y))
        box_x = (min = minimum(poly_x), max = maximum(poly_x))

        # Start by narrowing our list of unknowns. Any point that is outside the box defined
        # by the minimum and maximum is not in the polygon. Any point in a corner is in the
        # polygon
        for i in eachindex(data_y, data_x)
            # Check bounds. If outside the bounds, skip checking the corners
            # inpolygon[i] is false by default, so we don't have to change it
            if box_y.min > data_y[i] > box_y.max && box_x.min > data_x[i] > box_x.max
                known[i] = true
                continue
            end

            # Check corners
            for j in eachindex(poly_y, poly_x)
                if data_y[i]==poly_y[j] && data_x[i]==poly_x[j]
                    inpolygon[i] = known[i] = true
                end
            end
        end

        # The remaining points need to be determined using the even-odd algorithm. We're
        # going to work with the y coordinate, because this was originally designed for 
        # latitude / longitude coordinate points, and the distance between longitudes 
        # changes as you move poleward
        npoints = length(poly_y)
        for i in eachindex(data_x, data_y)
            # Don't bother if we already know
            known[i] && continue

            # Check each side
            pass = 0
            for j in eachindex(poly_y, poly_x)
                # Define the start and end points of the side
                A = (x = poly_x[j], y = poly_y[j])
                B = (x = poly_x[j%npoints+1], y = poly_y[j%npoints+1])

                # A ray extending infinitely in the positive direction along the latitude 
                # (x) axis will cross the side if and only if:
                    # 1) The latitude (y) coordinate of the ray falls between the latitude
                    #    (y) coordinates of the side
                    # 2) The starting longitude (x) coordinate of the ray is less than the
                    #    smallest longitude (x) coordinate of the side

                if (min(A.y,B.y) <= data_y[i] <= max(A.y,B.y)) && (data_x[i] <= min(A.x,B.x))
                    pass += 1
                end
            end

            # The point is in the polygon if it crosses through the sides an odd number of times
            inpolygon[i] = (pass % 2 != 0)
            known[i] = true
        end
        
        return data_x[inpolygon], data_y[inpolygon], inpolygon
    end


## --- Find coordinates in polygon
    """
    ```julia
    coords_in_shape(polylons::AbstractArray, polylats::AbstractArray, 
                    datalons::AbstractArray, datalats::AbstractArray)
    ```

    Return the latitude and longitude coordinates inside a polygon defined by `polylats`,
    `polylons`. Also returns a `BitVector` to index into the original data coordinate 
    vectors.

    ### Example

    ```julia
    lats, lons, points = coords_in_shape(polylon, polylat, datalon, datalat)
    ```
    """
    function coords_in_shape(polylons::AbstractArray{T1}, polylats::AbstractArray{T2}, 
        datalons::AbstractArray{T3}, datalats::AbstractArray{T4}) where {T1, T2, T3, T4}

        # Preallocate
        Tp = float(promote_type(T1, T2))
        Td = float(promote_type(T3, T4))

        xp, yp, zp = similar(polylats, Tp), similar(polylats, Tp), similar(polylats, Tp)
        xd, yd, zd = similar(datalats, Td), similar(datalats, Td), similar(datalats, Td)

        # Convert coordinates to cartesian
        @inbounds for i in eachindex(polylats, polylons)
            φ = deg2rad(90 - polylats[i])
            θ = deg2rad(polylons[i])
            xp[i], yp[i], zp[i] = cartesian(one(Tp), φ, θ)
        end
        @inbounds for i in eachindex(datalats, datalons)
            φ = deg2rad(90 - datalats[i])
            θ = deg2rad(datalons[i])
            xd[i], yd[i], zd[i] = cartesian(one(Td), φ, θ)
        end

        # Find points in the polygon. The z value is always constant, so we can ignore it
        xin, yin, inpolygon = points_in_shape(yp, xp, yd, xd)

        return datalons[inpolygon], datalats[inpolygon], inpolygon
    end


## --- End of file