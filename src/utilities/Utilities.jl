## --- Run all sub-utilities
    include("Definitions.jl")
    include("Slope.jl")
    include("NaNMeasurements.jl")
    include("Macrostrat.jl")
    include("Analysis.jl")


## --- Match rock names / types / descriptions to defined rock classes

    """
    ```julia
    match_rocktype(primary, secondary, tertiary; 
        source::Symbol, 
        [major::Bool], 
        [unmultimatch::Bool])
    ```

    Classify rock samples as sedimentary, igneous, or metamorphic (and associated subtypes)
    based on `primary`, `secondary`, and `tertiary` sample metadata. Use `primary` metadata 
    first; if no matches are made, attempt to classify sample using `secondary` metadata, 
    etc. 

    Specify the source of the samples as `:macrostrat` or `:earthchem`.

    ### Optional kwarg `major`
    `true` returns: `sed, ign, met`

    `false` returns: `siliciclast, shale, carb, chert, evaporite, coal, sed, volc, plut, 
    ign, metased, metaign, met, cover`

    Major rock types include subclasses; i.e. `ign` includes volcanic and plutonic samples.

    ### Optional kwarg `unmultimatch`
    Setting `unmultimatch=false` will not remove multiply-matched samples. Defaults to `true`.

    # Example
    ```julia
    cats = match_rocktype(rocktype, rockname, rockdescrip, source=:macrostrat, major=true)
    NamedTuple with 4 elements:
    sed    = BitVector(50000,)    [true ... true]
    ign    = BitVector(50000,)    [false ... false]
    met    = BitVector(50000,)    [false ... false]
    cover  = BitVector(50000,)    [false ... false]
    ```
    """
    match_rocktype(primary::AbstractArray, secondary::AbstractArray, tertiary::AbstractArray; 
        major::Bool=false, 
        unmultimatch::Bool=true,
        inclusive::Bool=true, 
        source::Symbol
    ) = _match_rocktype(primary, secondary, tertiary, major, unmultimatch, inclusive, static(source))

    """
    ```julia
    match_rocktype(rocktype, rockname, rockdescrip; source=:macrostrat, [major], [unmultimatch])
    ```

    Match Macrostrat rock names to defined rock classes.
    """
    function _match_rocktype(rocktype, rockname, rockdescrip, major, unmultimatch, inclusive,
            source::StaticSymbol{:macrostrat}
        )

        # Get rock type classifications and initialized BitVector
        typelist, cats = get_cats(major, length(rocktype))
        p = Progress(length(typelist)*4, desc="Finding Macrostrat rock types...")

        # Check major lithology 
        for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    cats[j][k] |= match(r"major.*?{(.*?)}", rocktype[k]) |> x -> 
                    isa(x,RegexMatch) ? containsi(x[1], typelist[j][i]) : false
                end
            end
            next!(p)
        end

        # Check the rest of rocktype
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rocktype[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # Then rockname
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rockname[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # Then rockdescrip
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= containsi(rockdescrip[k], typelist[j][i]))
                end
            end
            next!(p)
        end

        # If subtypes are true, major types must also be true
        if !major && inclusive
            minorsed, minorign, minormet = get_minor_types()
            for type in minorsed
                cats.sed .|= cats[type]
            end
            for type in minorign
                cats.ign .|= cats[type]
            end
            for type in minormet
                cats.met .|= cats[type]
            end
        end

        unmultimatch && return un_multimatch!(cats, major)
        return cats
    end

    """
    ```julia
    match_rocktype(Rock_Name, Type, Material; source=:earthchem, [major], [unmultimatch])
    ```

    Match Earthchem rock names to defined rock classes.
    """
    function _match_rocktype(Rock_Name, Type, Material, major, unmultimatch, inclusive,
            source::StaticSymbol{:earthchem}
        )

        # Get rock type classifications and initialized BitVector
        typelist, cats = get_cats(major, length(Rock_Name))
        p = Progress(length(typelist)*3, desc="Finding Earthchem rock types...")

        # Check rock name
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    cats[j][k] |= containsi(Rock_Name[k], typelist[j][i])
                end
            end
            next!(p)
        end

        # New typelist for Type
        # Omitted: vein
        typelist = (
            siliciclast = ("siliciclastic", "conglomerate&breccia",),
            volc = ("volcanic",),
            plut = ("plutonic", "pegmatitic"),
        )
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= (Type[k] == typelist[j][i]))
                end
            end
            next!(p)
        end

        # New typelist for Material
        # Omitted: vein, ore
        typelist = (
            sed = ("sedimentary",),
            plut = ("exotic",),
            ign = ("igneous", "xenolith",),
            met = ("metamorphic", "alteration",),
        )
        not_matched = find_unmatched(cats)
        @inbounds for j in eachindex(typelist)
            for i = eachindex(typelist[j])
                for k in eachindex(cats[j])
                    not_matched[k] && (cats[j][k] |= (Material[k] == typelist[j][i]))
                end
            end
            next!(p)
        end

        # If subtypes are true, major types must also be true
        if !major && inclusive
            minorsed, minorign, minormet = get_minor_types()
            for type in minorsed
                cats.sed .|= cats[type]
            end
            for type in minorign
                cats.ign .|= cats[type]
            end
            for type in minormet
                cats.met .|= cats[type]
            end
        end

        unmultimatch && return un_multimatch!(cats, major)
        return cats
    end


    """
    ```julia
    match_rocktype(writtentype::AbstractArray{String})
    ```

    Return a `NamedTuple` of `BitVector`s catagorizing Macrostrat samples as sedimentary, 
    igneous, metamorphic, and associated subtypes, or cover from types stored as strings 
    in `writtentype`.
    """
    function match_rocktype(writtentype::AbstractArray{String})
        cats = get_cats(false, length(writtentype))[2]

        # Parse all of the written types into cats
        for i in eachindex(writtentype)
            try
                cats[Symbol(writtentype[i])][i] = true
            catch
                continue
            end
        end

        # Define sed / ign / met as true when their subtypes are true
        minorsed, minorign, minormet = get_minor_types()
        minortypes = (
            sed = minorsed,
            ign = minorign,
            met = minormet
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
    find_unmatched(cats)
    ```

    Given a `Tuple` of `BitVectors`, return a `BitVector` that is `true` at index `i` iff 
    all elements of the `Tuple` are `false` at index `i`.

    If `cats` is a `NamedTuple` of rock types defined by `get_cats`, specify `major` 
    as `true` or `false` to decrease runtime. `major` is `true` if `cats` contains only 
    `sed`, `ign`, `met`, and `cover`.

    # Example
    ```julia-repl
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
    function find_unmatched(cats)
        matched = falses(length(cats[1]))
        @inbounds for i in eachindex(cats)
            matched .|= cats[i]
        end
        return .!matched
    end


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
        typelist = get_rock_class(major=true)
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
                cats[i][j] |= match(r"major.*?{(.*?)}", rocktype[j]) |> x -> isa(x,RegexMatch) ? 
                containsi(x[1], typelist[i]) : false
            end
            next!(p)
        end

        # Check the rest of rocktype
        not_matched = find_unmatched(cats)
        @inbounds for i in eachindex(typelist)
            for j in eachindex(rocktype)
                not_matched[j] && (cats[i][j] |= containsi(rocktype[j], typelist[i]))
            end
            next!(p)
        end

        # Then rockname
        not_matched = find_unmatched(cats)
        @inbounds for i in eachindex(typelist)
            for j in eachindex(rockname)
                not_matched[j] && (cats[i][j] |= containsi(rockname[j], typelist[i]))
            end
            next!(p)
        end

        # Then rockdescrip
        not_matched = find_unmatched(cats)
        @inbounds for i in eachindex(typelist)
            for j in eachindex(rockdescrip)
                not_matched[j] && (cats[i][j] |= containsi(rockdescrip[j], typelist[i]))
            end
            next!(p)
        end

        return cats
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


## --- Metadata for a specific sample

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

    function _get_type(cats, i, all_keys::True)
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
    class_up(typelist, name::String; 
        [all_types:Bool=false]
    )
    ```

    Get the rock type which describes the rock name `name`. The rock name must be an
    _exact_ match with the names in `typelist`. 
    
    Optionally specify `all_types` as `true` to return all types which match with `name`: 
    some names may be mapped to more than one rock type.

    See also: `get_cats` to return `typelist`.

    # Examples
    ```julia-repl
    julia> class_up(typelist, "dacit")
    :volc

    julia> class_up(typelist, "evaporite")
    :evaporite

    julia> class_up(typelist, "pyroxenite")
    :plut

    julia> class_up(typelist, "pyroxenite", all_types=true)
    (:plut, :met)
    ```
    """
    function class_up(typelist, name::String; all_types::Bool=false)
        if all_types
            rocktypes = keys(typelist)
            keymatches = falses(length(rocktypes))

            for k in eachindex(rocktypes)
                for i in eachindex(typelist[rocktypes[k]])
                    keymatches[k] |= (name == typelist[rocktypes[k]][i])
                end
            end

            return ifelse(count(keymatches)==0, nothing, rocktypes[keymatches])
        else
            @inbounds for k in keys(typelist)
                for i in typelist[k]
                    name == i && return k
                end
            end
        end

        return nothing
    end

    """
    ```julia
    class_up(name::Symbol, minorsed, minorign, minormet)
    ```

    Get the major rock type that describes the minor type `name`.

    # Examples
    ```julia-repl
    julia> class_up(:volc, minorsed, minorign, minormet)
    :ign

    julia> class_up(:ign, minorsed, minorign, minormet)
    :ign
    ```
    """
    function class_up(name::Symbol, minorsed, minorign, minormet)
        if name==:sed || name==:ign || name==:met
            return name
        end
    
        @inbounds for k in minorsed
            name==k && return :sed
        end
        @inbounds for k in minorign
            name==k && return :ign 
        end
        @inbounds for k in minormet
            name==k && return :met 
        end

        return nothing
    end

    """
    ```julia
    majorminor(type::Tuple, minorsed, minorign, minormet)
    ```

    Classify `type` as major or minor types. 

    ### Return
     *  `maj`: Vector of all types converted to major types.
     *  `ismaj`: BitVector which is `true` where `type` is a major type.
     *  `ctypes`: NamedTuple with the number of minor types present for each major type.

    # Example
    ```julia-repl
    julia> minorsed, minorign, minormet = get_minor_types();

    julia> type=(:ign, :volc);

    julia> maj, ismaj, ctypes = majorminor(type, minorsed, minorign, minormet)
    ([:ign, :ign], Bool[1, 0], (ign = 1,))
    ```

    Note that since `type` includes one major type (`:ign`) and one minor type (`:volc`),
    where `:volc` is a subtype of the `:ign` major type, `type` contains one minor type,
    so `ctypes` contains one element (`ctypes.ign`), equal to 1.
    """
    function majorminor(type::Tuple, minorsed, minorign, minormet)
        # Preallocate
        maj = similar(collect(type))
        ismaj = falses(length(maj))

        # Assign
        for j in eachindex(type)
            maj[j] = class_up(type[j], minorsed, minorign, minormet)
            ismaj[j] = maj[j] == type[j]
        end

        utype = unique(maj)
        ctypes = NamedTuple{Tuple(utype)}([count(==(i), maj) for i in utype] .- 1)

        return maj, ismaj, ctypes
    end


## --- Geochemistry

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


## --- Sample matching

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
        ll_total = Array{Float64}(undef, npoints, 1)

        # Replace missing values: this will penalize but not exclude missing data
        @inbounds for i in 1:npoints
            if isnan(bulkage[i])
                bulkage[i] = ifelse(sampleage < 1900, 3800, 0)
            end

            # Assume if one coordinate is missing, so is the other one
            if isnan(bulklat[i])
                bulklat[i] = -samplelat
                bulklon[i] = samplelon + 180
            end
        end

        @turbo for i in 1:npoints
            # Age (σ = 38 Ma)
            ll_age[i] = -((bulkage[i] - sampleage)^2)/(38^2)

            # Distance (σ = 1.8 arc degrees)
            ll_dist[i] = -((haversine(samplelat, samplelon, bulklat[i], bulklon[i]))^2)/(1.8^2)
        end

        @. ll_total = ll_age + ll_dist
        # ll_total .= 0
        
        # Geochemical log-likelihoods
        # for elem in eachindex(bulkgeochem)
        #     @turbo for i in 1:npoints
        #         ll_total[i] += -((bulkgeochem[elem][i] - samplegeochem[elem].m)^2)/(samplegeochem[elem].e^2)
        #     end
        # end

        matched_sample = rand_prop_liklihood(ll_total)
        return sampleidx[matched_sample]
    end


    """
    ```julia
    rand_prop_liklihood(ll)
    ```

    Weighted-random selection of an index based on log-likelihoods `ll`.

    See also: `weighted_rand` for weighted random selection based on likelihood / weight
    in normal space.
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


    """
    ```julia
    weighted_rand(p)
    ````

    Weighted random selection of an index based on weights `p`.

    See `rand_prop_liklihood` for selection based on log-likelihood.
    """
    function weighted_rand(p)
        sum_weights = nansum(p)
        r = rand()*sum_weights
        s = zero(typeof(sum_weights))
        @inbounds for i in eachindex(p)
            s += p[i]
            if s > r
                return i
            end
        end
        return lastindex(p)
    end


    """
    ```julia
    replace_major(type::Tuple, minortypes::NamedTuple, p::NamedTuple)
    ```

    Replace `:sed`, `:ign`, and `:met` types in `type` with a weighted-random selection 
    of a subtype from `minortypes` based on the weights in `p`. Remove duplicate types.

    Weights `p` must have the keys `sed`, `ign`, and `met`.

    ### Kwargs `ign`, `sed`, and `met`

    Minor rock types classified under the igneous, sedimentary, and metamorphic major 
    types.

    # Example
    ```julia-repl
    julia> minortypes = (
            sed=(:carb, :shale), 
            ign=(:volc, :plut), 
            met=(:metased, :metaign)
        );

    julia> p = (
            sed = [0.6, 0.4],
            ign = [0.5, 0.5],
            met = [0.3, 0.7],
        );

    julia> type = (:ign, :met, :carb, :plut,);

    julia> replace_major(type, p, minortypes)
    (:volc, :metaign, :carb, :plut)
    ```
    """
    function replace_major(type::Tuple, p::NamedTuple, minortypes::NamedTuple)
        type = collect(type)
        for i in eachindex(type)
            if type[i]==:ign 
                type[i] = minortypes.ign[weighted_rand(p.ign)]
            elseif type[i]==:sed
                type[i] = minortypes.sed[weighted_rand(p.sed)]
            elseif type[i]==:met 
                type[i] = minortypes.met[weighted_rand(p.met)]
            end
        end

        return Tuple(unique(type))
    end


    """
    ```julia
    get_descriptive_name(samplenames::Tuple, p_name::NamedTuple, sampletypes::Tuple,
            p_type::NamedTuple, typelist::NamedTuple, minortypes::NamedTuple
        )
    ```

    Given matched rock names `samplenames` and matched rock types `sampletypes`, return a 
    list of minor types and a randomly-selected descriptive rock name.

    # Example
    ```julia-repl
    name, types = get_descriptive_name(samplenames, p_name, sampletypes, p_type, typelist, minortypes)
    ```
    """
    function get_descriptive_name(samplenames::Tuple, p_name::NamedTuple, sampletypes::Tuple,
            p_type::NamedTuple, typelist::NamedTuple, minortypes::NamedTuple
        )

        # Randomly pick a matched type
        i = rand(1:length(sampletypes))
        t = sampletypes[i]

        if t==:sed || t==:ign || t==:met 
        # If the type is a major type, replace type and pick a rock name from mapped names
            t = minortypes[t][weighted_rand(p_type[t])]
            name = typelist[t][weighted_rand(p_name[t])]

            # Replace any other majors, but keep the already-replaced type
            sampletypes = collect(sampletypes)
            sampletypes[i] = t
            sampletypes = replace_major(Tuple(sampletypes), p_type, minortypes)

            return Symbol(name), sampletypes
            
        else
        # If the type is a minor type, randomly pick a matched name that maps to that type
            q = [string(n) in typelist[t] for n in samplenames]

            return rand(samplenames[q]), replace_major(sampletypes, p_type, minortypes)
        end 
    end
    

## --- Measurements

    """
    ```julia
        unmeasurementify(A::AbstractArray{Measurement{Float64}})
        unmeasurementify(A::NamedTuple)
    ```
    
    Separate an Array or NamedTuple `A` of `measurements` into an array of values and an 
    array of errors.
    
    # Example
    ```julia-repl
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


## --- (Coordinate) point in polygon

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

    # Example
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


    """
    ```julia
    coords_in_shape(polylons::AbstractArray, polylats::AbstractArray, 
                    datalons::AbstractArray, datalats::AbstractArray)
    ```

    Return the latitude and longitude coordinates inside a polygon defined by `polylats`,
    `polylons`. Also returns a `BitVector` to index into the original data coordinate 
    vectors.

    # Example
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