## --- Match rock type / class
    """
    ```julia
    match_rocktype(rocktype, rockname, rockdescrip; 
        [major], 
        [showprogress]
    )
    ```

    Match Macrostrat rock names to defined rock classes.

    Classify rock samples as sedimentary, igneous, or metamorphic (and associated subtypes)
    based on `rocktype`, `rockname`, and `rockdescrip` sample metadata. Use `rocktype` 
    first; if no matches are made, attempt to classify sample using `rockname` metadata, 
    etc. 

    Set `major=true` to return only matches for sedimentary, igneous, and metamorphic 
    rocks. See `get_rock_class` for discussion of rock types and subtypes. Set 
    `showprogress=false` to suppress the progress bar.

    # Example
    ```julia
    cats = match_rocktype(rocktype, rockname, rockdescrip, major=true)
    NamedTuple with 4 elements:
    sed    = BitVector(50000,)    [true ... true]
    ign    = BitVector(50000,)    [false ... false]
    met    = BitVector(50000,)    [false ... false]
    cover  = BitVector(50000,)    [false ... false]
    ```

    """
    function match_rocktype(rocktype::T, rockname::T, rockdescrip::T; major::Bool=false,
        showprogress::Bool=true) where T <: AbstractArray{<:String}
        
        # Get rock type classifications and initialized BitVector
        typelist, cats = get_cats(major, length(rocktype))
        set = keys(typelist)

        p = Progress(length(typelist)*4, desc="Finding Macrostrat rock types...", enabled=showprogress)

        # Check major lithology 
        for s in set
            for i in eachindex(typelist[s])
                cats[s] .|= (match.(r"major.*?{(.*?)}", rocktype) .|> 
                    x -> isa(x, RegexMatch) ? containsi.(x[1], typelist[s][i]) : false)
            end
            next!(p)
        end

        # Check the rest of rocktype
        not_matched = find_unmetamorphosed_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(rocktype[not_matched], i)
            end
            next!(p)
        end

        # Then rockname
        not_matched = find_unmetamorphosed_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(rockname[not_matched], i)
            end
            next!(p)
        end

        # Then rockdescrip
        not_matched = find_unmetamorphosed_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(rockdescrip[not_matched], i)
            end
            next!(p)
        end

        return rm_false_positives!(cats)
    end


    """
    ```julia
    match_rocktype(Rock_Name, Type, Material, sedrocks, ignrocks)
    ```

    Match Earthchem rock names to defined rock classes. Specify rock subtypes of sedimentary
    and igneous rocks as `sedrocks` and `ignrocks`. Recommended:

    ```julia
    sedrocks = (minorsed..., :sed,)
    ignrocks = (minorvolc..., minorplut..., minorign..., :ign)
    ```

    Rocks will be categorized as minor subtypes; there is no option to match only major rock 
    types.
    """
    function match_rocktype(Rock_Name::T, Type::T, Material::T, sedrocks, ignrocks) where T <: AbstractArray{<:String}
        # Get rock type classifications and initialized BitVector
        typelist, cats = get_cats(false, length(Rock_Name));

        # Broadly categorize samples as igneous and sedimentary based on material and type
        materials = (
            sed = ("sedimentary",),
            ign = ("exotic", "igneous", "xenolith",),
            met = ("metamorphic")
        )
        groups = NamedTuple{keys(materials)}([falses(length(Rock_Name)) for _ in 1:length(materials)])
        for k in keys(materials)
            for i in materials[k]
                groups[k] .|= (Material .== i)
            end
        end

        types = (
            sed = ("siliciclastic", "conglomerate&breccia",),
            ign = ("volcanic", "plutonic", "pegmatitic",),
        )
        for k in keys(types)
            for i in types[k]
                groups[k] .|= (Type .== i)
            end
        end

        # Match samples labeled as sedimentary / igneous with sedimentary / igneous rocks
        # This should stop "pyroclastic / clastic," etc. false positives
        match_subset!(cats, groups.sed, typelist, sedrocks, Rock_Name)
        match_subset!(cats, groups.ign, typelist, ignrocks, Rock_Name)

        # Siliciclastic type samples are named only "sedimentary" and need to be specifically
        # set as siliciclasts
        t = vec((Type .== "siliciclastic") .| (Type .== "conglomerate&breccia"));
        cats.siliciclast[t] .= true

        # Volcanic breccia isn't matched since the name "breccia" is sedimentary: it must 
        # be specifically set as volcanic
        t = vec((Rock_Name .== "breccia") .& (Type .== "volcanic"));
        cats.volc[t] .= true

        # Allow all unmatched samples to match against all rock types
        not_matched = find_unmatched(cats);
        match_subset!(cats, not_matched, typelist, keys(typelist), Rock_Name)

        # Match all unmatched samples with a defined type
        not_matched = find_unmatched(cats);

        t = vec(Type .== "volcanic");
        cats.volc[t .& not_matched] .= true

        t = vec((Type .== "plutonic") .| (Type .== "pegmatitic"));
        cats.plut[t .& not_matched] .= true

        t = vec(Type .== "siliceous")
        cats.chert[t .& not_matched] .= true

        # Match all unmatched samples with a defined sedimentary / igneous material 
        not_matched = find_unmatched(cats);
        cats.sed[groups.sed .& not_matched] .= true
        cats.ign[groups.ign .& not_matched] .= true

        # Match all unmatched samples with a defined metamorphic material
        # This means that protoliths should be mostly known, even though we didn't use
        # find_unmetamorphosed_unmatched: we don't have to use it here, because we know
        # the values of Material and Type
        not_matched = find_unmatched(cats);
        cats.met[groups.met .& not_matched] .= true

        return rm_false_positives!(cats)
    end


    """
    ```julia
    match_rocktype(rock_name, sample_description, qap_name, rgroup_id;
        rockgroup_id, rockgroup_name,
        [showprogress]
    )
    ```

    Match rock names from Gard et al., 2019 (10.5194/essd-11-1553-2019) to defined rock 
    classes. Attempts to match names based on information given by the original authors
    (`rock_name` and `sample_description`) before using computed or assigned rock classes
    (`qap_name` and `rgroup_id`).

    ## Required kwargs 

    `rockgroup_id`: rock group ID numbers from the rockgroup.csv file. Must be integers.

    `rockgroup_name`: A _single_ string of rock group descriptions (group, origin, and 
    facies) from the rockgroup.csv file.

    ## Optional kwargs
    `showprogress`: enable or disable the progress bar. Bool, defaults to `true`.

    """
    function match_rocktype(rock_name::T, sample_description::T, qap_name::T, rgroup_id::AbstractArray{Int};
            rockgroup_id::AbstractArray{Int}, rockgroup_name::T,
            showprogress::Bool=true
        ) where T <: AbstractArray{String}

        # Get rock type classifications and initialized BitVector
        typelist, cats = get_cats(false, length(rock_name))
        set = keys(typelist)

        p = Progress(length(typelist)*3,
            desc="Finding Gard et al., rock types...", enabled=showprogress
        )

        # Check rock name designated by original authors 
        for s in set
            for i in typelist[s]
                cats[s] .|= containsi.(rock_name, i)
            end
            next!(p)
        end

        # Check sample description inherited from previous databases
        not_matched = find_unmetamorphosed_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(sample_description[not_matched], i)
            end
            next!(p)
        end

        # Check QAP name
        not_matched = find_unmetamorphosed_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(qap_name[not_matched], i)
            end
            next!(p)
        end

        # Assign unmatched rocks by rock ID
        not_matched = find_unmetamorphosed_unmatched(cats)

        id_cats = get_cats(false, length(rockgroup_id))[2]
        for s in set
            for i in typelist[s]
                id_cats[s] .|= containsi.(rockgroup_name, i)
            end
        end
        
        for i in eachindex(rgroup_id)
            if not_matched[i]
                for s in set
                    cats[s][i] |= id_cats[s][rgroup_id[i]]
                end
            end
        end

        # For sedimentary rock group ID 88 (clastics), remove soils and cover
        cover = ("silt", "sand", "muddy", "mud", "gravel", "clay", "soil", "soils", 
            "andosol", "loam", "sediment", "ooze")
        for i in eachindex(rgroup_id)
            if rgroup_id[i] == 88
                # Split into groups of whole words
                samplenames = replace(rock_name[i], "/" => " ", "-" => " ", "." => " ", 
                    "," => " ", "(" => " ", ")"=> " ", 
                )
                samplenames = split(samplenames, " ")

                # Check if any of the words match any cover. If yes, remove all matches
                # Matches with clastics should have already been caught. Also, their 
                # clastics don't really match my clastics, because I differentiate 
                # siliciclastics and fine-grained siliciclastics
                for n in samplenames
                    for k in cover
                        if n == k
                            [cats[s][i] = false for s in set]
                            cats.cover[i] = true
                            break
                        end
                    end
                end
            end

        end
        
        # Remove false positives and return
        return rm_false_positives!(cats)
    end


    """
    ```julia
    match_rocktype(writtentype::AbstractArray{String})
    ```

    Return a `NamedTuple` of `BitVector`s catagorizing Macrostrat samples as sedimentary, 
    igneous, metamorphic, and associated subtypes, or cover from types stored as strings 
    in `writtentype`.

    Major types do not include minor subtypes.
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

        return cats
    end
    export match_rocktype


## --- Supporting, nonexported functions for match_rocktype
    """
    ```julia
    match_subset!(cats, filter, typelist, subset, rockname)
    ```

    Match a `subset` of rock type names in `typelist` to rock names in `rockname`. Filter
    the rock names to be matched by `filter`. The types in `subset` should be the same
    as the keys in `cats` and `typelist` which you want to search.

    To get `cats` and `typelist`, see `get_cats` or `get_rock_class`.

    # Example
    ```julia
    julia> typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    julia> typelist, cats = get_cats(false, length(bulkrockname));

    julia> sedlist = (minorsed..., :sed,);

    julia> match_subset!(cats, groups.sed, typelist, sedlist, bulkrockname)
    ```

    """
    function match_subset!(cats::NamedTuple, rfilter::BitVector, typelist::NamedTuple, 
            subset, rockname::AbstractArray{<:String}
        )

        for s in subset
            for i in typelist[s]
                cats[s][rfilter] .|= containsi.(rockname[rfilter], i)
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

    """
    function find_unmatched(cats)
        matched = falses(length(cats[1]))
        @inbounds for i in eachindex(cats)
            matched .|= cats[i]
        end
        return .!matched
    end

    """
    ```julia
    find_unmetamorphosed_unmatched(cats)
    ```

    As `find_unmatched`, but does not count uncategorized metamorphic rock as a matched
    type. This allows `match_rocktype` to continue searching rockname and rockdescrip for
    potential protoliths.
    """
    function find_unmetamorphosed_unmatched(cats)
        matched = falses(length(cats[1]))
        @inbounds for k in keys(cats)
            k == :met && continue
            matched .|= cats[k]
        end 
        return .!matched
    end

    """
    ```julia
    rm_false_positives!(cats)
    ```

    Remove false matches from `cats`, where a false match is when a rock name used to
    identify samples is present in another rock name. For example, all grano*diorite* 
    samples will match with *diorite*.

    Currently removes:
     *  Diorite from grano*diorite*.
     *  Lignite from ma*lignite*.
    """
    function rm_false_positives!(cats)
        cats.diorite .&= cats.granodiorite  # Diorite / granodiorite
        cats.coal .&= .!cats.alk_volc       # Lignite / malignite 

        return cats
    end


## --- Find metamorphic rocks 
    """
    ```julia
    find_metamorphics(rocktype, rockname, rockdescrip; 
        [major], 
        [showprogress]
    )
    ```

    As `match_rocktype`, but searches Macrostrat rock names for metamorphic rocks. 
    Matches are set to the protoliths of metamorphic rocks (e.g., a `true` value for 
    `:shale` refers to a metamorphosed shale), or to `:met` for undifferentiated 
    metamorphic rocks.

    """
    function find_metamorphics(rocktype::T, rockname::T, rockdescrip::T; major::Bool=false,
        showprogress::Bool=true) where T <: AbstractArray{<:String}
        
        # Get rock type classifications and initialized BitVector
        typelist = get_metamorphic_class()
        cats = NamedTuple{keys(typelist)}([falses(length(rocktype)) for _ in 1:length(typelist)]) 
        set = keys(typelist)[collect(.!isempty.(values(typelist)))]

        p = Progress(length(set)*4, desc="Finding metamorphic Macrostrat samples...", enabled=showprogress)

        # Check major lithology 
        for s in set
            for i in eachindex(typelist[s])
                cats[s] .|= (match.(r"major.*?{(.*?)}", rocktype) .|> 
                    x -> isa(x, RegexMatch) ? containsi.(x[1], typelist[s][i]) : false)
            end
            next!(p)
        end

        # Check the rest of rocktype
        not_matched = find_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(rocktype[not_matched], i)
            end
            next!(p)
        end

        # Then rockname
        not_matched = find_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(rockname[not_matched], i)
            end
            next!(p)
        end

        # Then rockdescrip
        not_matched = find_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(rockdescrip[not_matched], i)
            end
            next!(p)
        end

        return rm_false_positives!(cats)
    end


    """
    ```julia
    find_metamorphics(rock_name, sample_description, qap_name, rgroup_id;
        rockgroup_id, rockgroup_name,
        [showprogress]
    )
    ```

    As `match_rocktype`, but searches Gard et al., 2019 (10.5194/essd-11-1553-2019) for 
    metamorphic rocks.

    """
    function find_metamorphics(rock_name::T, sample_description::T, qap_name::T, rgroup_id::AbstractArray{Int};
        rockgroup_id::AbstractArray{Int}, rockgroup_name::T,
        showprogress::Bool=true
    ) where T <: AbstractArray{String}

        # Get rock type classifications and initialized BitVector
        typelist = get_metamorphic_class()
        cats = NamedTuple{keys(typelist)}([falses(length(rock_name)) for _ in 1:length(typelist)]) 
        set = keys(typelist)[collect(.!isempty.(values(typelist)))]

        p = Progress(length(set)*3,
            desc="Finding metamorphic rocks in Gard, et al...", enabled=showprogress
        )

        # Check rock name designated by original authors 
        for s in set
            for i in typelist[s]
                cats[s] .|= containsi.(rock_name, i)
            end
            next!(p)
        end

        # Check sample description inherited from previous databases
        not_matched = find_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(sample_description[not_matched], i)
            end
            next!(p)
        end

        # Check QAP name
        not_matched = find_unmatched(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(qap_name[not_matched], i)
            end
            next!(p)
        end

        # Assign unmatched rocks by rock ID
        not_matched = find_unmatched(cats)

        id_cats = get_cats(false, length(rockgroup_id))[2]
        for s in set
            for i in typelist[s]
                id_cats[s] .|= containsi.(rockgroup_name, i)
            end
        end
        
        for i in eachindex(rgroup_id)
            if not_matched[i]
                for s in set
                    cats[s][i] |= id_cats[s][rgroup_id[i]]
                end
            end
        end
        
        # Remove false positives and return
        return rm_false_positives!(cats)
    end
    export find_metamorphics


## --- End of file