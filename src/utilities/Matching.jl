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
        npoints = length(rocktype)
        typelist, cats = get_cats(major, npoints);
        set = keys(typelist);

        # Define descriptive and nondescriptive types
        is_nondescript = (:sed, :ign, :volc, :plut, :met, :cover)
        minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];
        minorset = (;
            sed = (minorsed..., :sed), 
            volc=(minorvolc..., :volc,), 
            plut=(minorplut..., :plut), 
            ign=(minorign..., :ign), 
            met=(:met,),
            cover=(:cover,),
        )

        # I fear progress
        p = Progress(length(set) + 3*length(is_nondescript), 
            desc="Finding Macrostrat rock types...", enabled=showprogress
        );

        # Start normally, checking major lithologies 
        for s in set
            for i in eachindex(typelist[s])
                cats[s] .|= (match.(r"major.*?{(.*?)}", rocktype) .|> 
                    x -> isa(x, RegexMatch) ? containsi.(x[1], typelist[s][i]) : false)
            end
            next!(p)
        end

        # Pause!! For samples matched only to a nondescriptive key (e.g. sed -> "sedimentary")
        # we want to keep looking, but only look in the "sedimentary" lists. We'll 
        # create catmatrix, a filter for which categories we're looking in. Start by defining 
        # it's columns...
        catmatrix_col = NamedTuple{is_nondescript}(i for i in eachindex(is_nondescript))
        catmatrix = find_unmatched(cats, catmatrix_col, is_nondescript, npoints)

        # Check the rest of rocktype
        for bigtype in eachindex(catmatrix_col)
            bigtype_filter = catmatrix[:,catmatrix_col[bigtype]]
            for s in minorset[bigtype]
                for i in typelist[s]
                    cats[s][bigtype_filter] .|= containsi.(rocktype[bigtype_filter], i)
                end
            end
            next!(p)
        end

        # Rockname 
        catmatrix = find_unmatched(cats, catmatrix_col, is_nondescript, npoints)
        for bigtype in eachindex(catmatrix_col)
            bigtype_filter = catmatrix[:,catmatrix_col[bigtype]]
            for s in minorset[bigtype]
                for i in typelist[s]
                    cats[s][bigtype_filter] .|= containsi.(rockname[bigtype_filter], i)
                end
            end
            next!(p)
        end

        # Rockdescrip
        catmatrix = find_unmatched(cats, catmatrix_col, is_nondescript, npoints)
        for bigtype in eachindex(catmatrix_col)
            bigtype_filter = catmatrix[:,catmatrix_col[bigtype]]
            for s in minorset[bigtype]
                for i in typelist[s]
                    cats[s][bigtype_filter] .|= containsi.(rockdescrip[bigtype_filter], i)
                end
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
        # find_unmatched_useful: we don't have to use it here, because we know
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
        not_matched = find_unmatched_useful(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(sample_description[not_matched], i)
            end
            next!(p)
        end

        # Check QAP name
        not_matched = find_unmatched_useful(cats)
        for s in set
            for i in typelist[s]
                cats[s][not_matched] .|= containsi.(qap_name[not_matched], i)
            end
            next!(p)
        end

        # Assign unmatched rocks by rock ID
        not_matched = find_unmatched_useful(cats)

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
    match_rocktype(Rock_Group, Rock_Subgroup, Rock_Composition, Rock_Name)
    ```
    """
    function match_rocktype(Rock_Group::T, Rock_Subgroup::T, Rock_Composition::T, 
        Rock_Facies::T, Rock_Name::T) where T <: AbstractArray{String}

        # Get rock type classifications, including metamorphic, and initialized BitVectors
        typelist, cats = get_cats(false, length(Rock_Name))
        metatypelist = get_metamorphic_class()
        metacats = deepcopy(cats)

        # Assign major lithologic class 
        metacats.met .= Rock_Group .== "metamorphic"
        metacats.sed .= ((Rock_Subgroup .== "metasedimentary") .| 
            (metacats.met .& (Rock_Subgroup .== "clastic")))
        metacats.ign .= Rock_Subgroup .== "metaigneous"
        metacats.volc .= Rock_Subgroup .== "metavolcanic"
        metacats.plut .= Rock_Subgroup .== "metaplutonic"

        cats.ign .= (Rock_Group .== "igneous") .| metacats.ign;
        cats.met .= (Rock_Group .== "metamorphic");
        cats.sed .= (Rock_Group .== "sedimentary") .| metacats.sed;
        cats.volc .= (Rock_Subgroup .== "volcanic") .| metacats.volc;
        cats.plut .= (Rock_Subgroup .== "plutonic") .| metacats.plut;

        # Match by major lithology 
        minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];

        match_subset!(cats, cats.sed, typelist, (minorsed..., :sed,), Rock_Name);
        match_subset!(cats, cats.ign, typelist, (minorign..., :ign,), Rock_Name);
        match_subset!(cats, cats.volc, typelist, (minorvolc..., :volc,), Rock_Name);
        match_subset!(cats, cats.plut, typelist, (minorplut..., :plut,), Rock_Name);

        match_subset!(metacats, metacats.sed, metatypelist, (minorsed..., :sed,), Rock_Name)
        match_subset!(metacats, metacats.ign, metatypelist, (minorign..., :ign,), Rock_Name)
        match_subset!(metacats, metacats.volc, metatypelist, (minorvolc..., :volc,), Rock_Name)
        match_subset!(metacats, metacats.plut, metatypelist, (minorplut..., :plut,), Rock_Name)

        # Ensure metamorphic samples are only those with no known protolith 
        include_minor!(cats)
        include_minor!(metacats)
        cats.met .&= .!(cats.ign .| cats.sed)
        metacats.met .&= .!(metacats.ign .| metacats.sed)

        # Let metamorphic samples match with whatever... Some misspellings or random rocks 
        # are named metased / metaign rock names but aren't tagged as such 
        match_subset!(cats, cats.met, typelist, keys(cats), Rock_Name);  
        match_subset!(metacats, metacats.met, typelist, keys(metacats), Rock_Name); 

        # Check rock facies (these will for sure be metamorphic)
        unmatched_cats = find_unmatched(cats)
        unmatched_metacats = find_unmatched(metacats)
        match_subset!(cats, unmatched_cats, typelist, keys(cats), Rock_Facies);  
        match_subset!(metacats, unmatched_metacats, typelist, keys(metacats), Rock_Facies);
        
        # And then re-exclude protoliths from metamorphic rocks
        include_minor!(cats)
        include_minor!(metacats)
        cats.met .&= .!(cats.ign .| cats.sed)
        metacats.met .&= .!(metacats.ign .| metacats.sed)

        # Remove soil samples 
        t = Rock_Composition .== "soil"
        for k in keys(cats)
            cats[k][t] .= false
            metacats[k][t] .= false
        end
        cats.cover[t] .= true

        # Remove false positives and return
        return rm_false_positives!(cats), metacats
    end

    
    """
    ```julia
    match_rocktype(writtentype::AbstractArray{String})
    ```

    Return a `NamedTuple` of `BitVector`s catagorizing samples as sedimentary, 
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
    rock_type_filter(set; [inclusiveign=false])
    ```

    Given a list of lithologic classes `set`, return a NamedTuple BitVector that filters 
    all the minor classes for a given major class.

    Set `inclusiveign=true` to make the BitVector filter for every single igneous class.

    # Example 
    ```julia-repl
    julia> typelist = get_rock_class()[1];

    julia> set = keys(typelist);

    julia> setfilter = rock_type_filter(set);

    julia> set[setfilter.sed]
    (:siliciclast, :shale, :carb, :evap, :chert, :phosphorite, :coal)

    julia> set[setfilter.sed] == minorsed
    true
    ```
    """
    function rock_type_filter(set; inclusiveign::Bool=false)
        # Preallocate
        majorlist = (:sed, :volc, :plut, :ign)
        setfilter = NamedTuple{majorlist}(falses(length(set)) for _ in majorlist)
        minorlist = zip(majorlist, get_rock_class()[2:end]);

        # Search and destroy
        for i in minorlist 
            for j in i[2]
                setfilter[i[1]] .|= (set .== j)
            end
        end

        if inclusiveign
            setfilter.ign .|= setfilter.volc
            setfilter.ign .|= setfilter.plut
            return setfilter
        else
            return setfilter
        end
    end
    export rock_type_filter


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
    find_unmatched_useful(cats)
    ```

    Imagine if `find_unmatched`, was more useful. 
        
    This finds all rocks that are unmatched to a minor / nondescriptive rock type. That means
    if a sample is matched with a... 
     * Undifferentiated metamorphic rock
     * Nondescriptive types (undifferentiated sedimentary, volcanic, plutonic, igneous)
    
    Then I count it as not matched and search the next category for useful information.

    So for example, something might have the rocktype "igneous" but specify "volcanic" in 
    the rockname, and "volcanic breccia and tuff" in the rock description. I don't want to 
    throw away the baby (tuff) in the bathwater (undifferentiated igneous).
    """
    function find_unmatched_useful(cats)
        matched = falses(length(cats[1]))

        @inbounds for k in keys(cats)
            k in (:sed, :volc, :plut, :ign, :met) && continue
            matched .|= cats[k]
        end 
        return .!matched
    end

    
    """
    ```julia
    find_unmatched(cats, catmatrix_col, is_nondescript, npoints)
    ```

    Return a `catmatrix` BitMatrix so for a given sample (row) each nondescriptive rock
    type (column) is true when that set of descriptive rocknames should be searched for 
    matches.

    # Arguments 
     * `cats`: A NamedTuple of BitVectors which is true when the corresponding sample 
       (index) has matched to a rock type.
     * `catmatrix_col`: Column indices of `catmatrix` which correspond to the nondescriptive 
       types.
     * `is_nondescript`: List of types considered nondescriptive. This should be equivalent 
       to the keys of `catmatrix_col`.
     * `npoints`: Number of samples in the dataset.

    # Output 
     * `catmatrix`: A BitMatrix with a row for each a sample, a column for each 
       nondescriptive type. A given (row, column) is true if that sample matched with 
       that nondescriptive type.

    # Method
     * If the sample matched to a descriptive type, all elements of the row are false. This 
       sample will not match with any more rock types.
     * If the sample did not match to *any* types, or *only* matched to metamorphic rocks, 
       all elements of the row are true. This sample may match with any defined rock type.
     * If the sample matched to *only* nondescriptive type(s), those types are true. In this 
       case, the metamorphic column is set to be false. Since the match to metamorphic samples 
       will not be undone, we don't need to search through the metamorphic names again. This 
       sample will only match with rocktypes that fall under the nondescriptive umbrella 
       (e.g., sedimentary rocks will only match with sedimentary subtypes)

    """
    function find_unmatched(cats::NamedTuple, catmatrix_col::NamedTuple, is_nondescript, npoints::Number)
        # Create catmatrix
        catmatrix = stack(cats[k] for k in is_nondescript)

        # Find matched and descriptively matched samples
        descriptive = falses(npoints);
        matched = falses(npoints);
        for k in keys(cats)
            matched .|= cats[k]
            k in is_nondescript && continue 
            descriptive .|= cats[k]
        end

        # If the sample matched with a descriptive type or with cover, we're DONE
        catmatrix[descriptive,:] .= 0
        catmatrix[cats.cover,:] .= 0

        # Samples that didn't match with anything, OR samples that ONLY matched with 
        # metamorphic rocks can match with anything
        only_metamorphic = vec(cats.met .& (sum(catmatrix, dims=2) .== 1) .& .!descriptive);
        catmatrix[only_metamorphic,:] .= 1;
        catmatrix[.!matched,:] .= 1;

        # If the sample matched with a nondescriptive type, we don't need to search through 
        # the metamorphic terms again. That match isn't gonna go away, so set that to zero.
        # Weird things happen with views, so we have to do it like this :(
        restricted = vec(0 .< sum(catmatrix, dims=2) .< 5);
        newmet = copy(catmatrix[:,catmatrix_col.met]);
        newmet[restricted] .= 0;
        catmatrix[:,catmatrix_col.met] .= newmet;

        return catmatrix
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
     *  Sedimentary (clastic) from volcani*clastic*.
    """
    function rm_false_positives!(cats)
        cats.diorite .&= .!cats.granodiorite  # granodiorite / diorite
        cats.coal .&= .!cats.alk_volc         # malignite / lignite 
        cats.sed .&= .!cats.volcaniclast      # volcaniclast / sed (clast)

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


## --- Convert NamedTuples to an array for file saving
    """
    ```julia
    cats_to_array(cats)
    ```

    Converts a NamedTuple of BitVectors `cats` to an array of 0s and 1s to save to a file. 
    Returns an array where each column corresponds to an element in `cats`, and a string of 
    header or row names.

    # Example 
    ```julia
    a, row_names_a = cats_to_array(cats)
    ``
    """
    function cats_to_array(cats::NamedTuple)
        a = Array{Int64}(undef, length(cats[1]), length(cats))
        for i in eachindex(keys(cats))
            for j in eachindex(cats[i])
                a[j,i] = ifelse(cats[i][j], 1, 0)
            end
        end

        return a, string.(collect(keys(cats))) 
    end
    export cats_to_array


## --- End of file