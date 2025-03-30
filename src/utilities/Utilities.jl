## --- Modify rock class BitVectors

    """
    ```julia
    delete_cover(cats)
    ```

    I don't like cover. It's coarse and rough and irritating, and it gets everywhere. 
    Not like `cats` after going through this function. Now the elements are only bedrock.
    """
    function delete_cover(cats)
        start = collect(keys(cats))
        notcover = start .!= :cover

        kittens = NamedTuple{Tuple(start[notcover])}(cats[k] for k in start[notcover])
        return kittens
    end
    export delete_cover


    """
    ```julia
    delete_volcaniclast(cats)
    ```

    As `delete_cover.`

    """
    function delete_volcaniclast(cats)
        start = collect(keys(cats))
        notvolcaniclast = start .!= :volcaniclast

        kittens = NamedTuple{Tuple(start[notvolcaniclast])}(cats[k] for k in start[notvolcaniclast])
        return kittens
    end
    export delete_volcaniclast


    """
    ```julia
    include_minor!(cats, [minorsed, minorvolc, minorplut, minorign])
    ```

    Include all minor rock classes in their major type. E.g., `cats.sed` is true when 
    `cats.shale` is true.

    Optionally specify `minorsed`, `minorvolc`, `minorplut`, and `minorign` for a slight
    speed-up. See `get_rock_class` for minor rock class lists.

    See also: `exclude_minor!`

    """
    function include_minor!(cats)
        minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];

        for type in minorsed
            cats.sed .|= cats[type]
        end
        for type in minorvolc
            cats.volc .|= cats[type]
        end
        for type in minorplut
            cats.plut .|= cats[type]
        end
        for type in (minorign..., minorvolc..., minorplut...)
            cats.ign .|= cats[type]
        end
        return cats
    end

    function include_minor!(cats, minorsed, minorvolc, minorplut, minorign)
        for type in minorsed
            cats.sed .|= cats[type]
        end
        for type in minorvolc
            cats.volc .|= cats[type]
        end
        for type in minorplut
            cats.plut .|= cats[type]
        end
        for type in (minorign..., minorvolc..., minorplut...)
            cats.ign .|= cats[type]
        end
        return cats
    end
    export include_minor!

    """
    ```julia
    exclude_minor!(cats, [minorsed, minorvolc, minorplut, minorign])
    ```

    Exclude all minor rock classes in their major type. E.g., `cats.sed` is false when 
    `cats.shale` is true. 

    Optionally specify `minorsed`, `minorvolc`, `minorplut`, and `minorign` for a slight
    speed-up. See `get_rock_class` for minor rock class lists.

    See also: `include_minor!`

    """
    function exclude_minor!(cats)
        minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:5];

        for type in minorsed
            cats.sed .&= .!cats[type]
        end
        for type in minorvolc
            cats.volc .&= .!cats[type]
        end
        for type in minorplut
            cats.plut .&= .!cats[type]
        end
        for type in (minorign..., minorvolc..., minorplut...)   # ALL igneous subclasses
            cats.ign .&= .!cats[type]
        end
        return cats
    end

    function exclude_minor!(cats, minorsed, minorvolc, minorplut, minorign)
        for type in minorsed
            cats.sed .&= .!cats[type]
        end
        for type in minorvolc
            cats.volc .&= .!cats[type]
        end
        for type in minorplut
            cats.plut .&= .!cats[type]
        end
        for type in (minorign..., minorvolc..., minorplut...)
            cats.ign .&= .!cats[type]
        end
        return cats
    end
    export exclude_minor!


## --- Metadata for a specific sample

    """
    ```julia
    get_type(cats::NamedTuple, i; [all_keys=false])
    ```

    Return the first key of `cats` where `i` is `true`. Optionally specify `allkeys`=`true`
    to return _all_ keys where `i` is true. Assumes equal length elements of `cats`.

    If there are no matches, return `:none`.

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

        return :none
    end

    function _get_type(cats, i, all_keys::True)
        @assert 0 < i <= length(cats[1]) "Index $i out of bounds."

        catkeys = keys(cats)
        keymatches = falses(length(catkeys))

        @inbounds for k in eachindex(catkeys)
            cats[k][i] && (keymatches[k]=true)
        end

        count(keymatches)==0 && return (:none,)
        return catkeys[keymatches]
    end
    export get_type


    # """
    # ```julia
    # class_up(typelist, name::String; 
    #     [all_types:Bool=false]
    # )
    # ```

    # Get the rock type which describes the rock name `name`. The rock name must be an
    # _exact_ match with the names in `typelist`. 
    
    # Optionally specify `all_types` as `true` to return all types which match with `name`: 
    # some names may be mapped to more than one rock type.

    # For `typelist` generated with `get_rock_class`, setting `inclusive=false` will not 
    # change the behavior of `alltypes=false`, as long as major types are listed after minor
    # types (because `class_up` will return the type for the first match it finds). If
    # `inclusive=true` and `alltypes=true`, the major type will always be returned.

    # See also: `get_rock_class` to return `typelist`.

    # # Examples
    # ```julia-repl
    # julia> class_up(typelist, "dacit")
    # :volc

    # julia> class_up(typelist, "evaporite")
    # :evaporite

    # julia> class_up(typelist, "pyroxenite")
    # :plut

    # julia> class_up(typelist, "pyroxenite", all_types=true)
    # (:plut, :met)
    # ```
    # """
    # function class_up(typelist, name::String; all_types::Bool=false)
    #     if all_types
    #         rocktypes = keys(typelist)
    #         keymatches = falses(length(rocktypes))

    #         for k in eachindex(rocktypes)
    #             for i in eachindex(typelist[rocktypes[k]])
    #                 keymatches[k] |= (name == typelist[rocktypes[k]][i])
    #             end
    #         end

    #         return ifelse(count(keymatches)==0, nothing, rocktypes[keymatches])
    #     else
    #         @inbounds for k in keys(typelist)
    #             for i in typelist[k]
    #                 name == i && return k
    #             end
    #         end
    #     end

    #     @warn "No upwards class found for $name"
    #     return nothing
    # end


    """
    ```julia
    class_up(name::Symbol, minorsed, minorign)
    ```

    Get the major rock type that describes the minor type `name`.

    _**Important:**_ If `minorign` contains only `(:volc, :plut, :carbonatite)`, `class_up`
    cannot recognize minor types within those subtypes.

    # Examples
    ```julia-repl
    julia> class_up(:volc, minorsed, minorign)
    :ign

    julia> class_up(:ign, minorsed, minorign)
    :ign

    julia> class_up(:dacite, minorsed, (minorign..., minorvolc..., minorplut...))
    ```
    """
    function class_up(name::Symbol, minorsed, minorign)
        if name==:sed || name==:ign || name==:met
            return name
        end
    
        @inbounds for k in minorsed
            name==k && return :sed
        end
        @inbounds for k in minorign
            name==k && return :ign 
        end

        @warn "No upwards class found for $name"
        return nothing
    end

    """
    ```julia
    class_up(name::Symbol, minorsed, minorvolc, minorplut, minorign)
    ```

    Get the major rock type that describes the minor type `name`.

    # Examples
    ```julia-repl
    julia> class_up(:volc, minorsed, minorign)
    :ign

    julia> class_up(:ign, minorsed, minorign)
    :ign
    ```
    """
    function class_up(name::Symbol, minorsed, minorvolc, minorplut, minorign)
        if name==:sed || name==:ign || name==:met || name==:volc || name==:plut
            return name
        end
    
        @inbounds for k in minorsed
            name==k && return :sed
        end
        @inbounds for k in minorvolc
            name==k && return :volc
        end
        @inbounds for k in minorplut
            name==k && return :plut
        end
        @inbounds for k in minorign
            name==k && return :ign 
        end

        @warn "No upwards class found for $name"
        return nothing
    end
    export class_up

    # """
    # ```julia
    # majorminor(type::Tuple, minorsed, minorign, minormet)
    # ```

    # Classify `type` as major or minor types. 

    # ### Return
    #  *  `maj`: Vector of all types converted to major types.
    #  *  `ismaj`: BitVector which is `true` where `type` is a major type.
    #  *  `ctypes`: NamedTuple with the number of minor types present for each major type.

    # # Example
    # ```julia-repl
    # julia> minorsed, minorign, minormet = get_minor_types();

    # julia> type=(:ign, :volc);

    # julia> maj, ismaj, ctypes = majorminor(type, minorsed, minorign, minormet)
    # ([:ign, :ign], Bool[1, 0], (ign = 1,))
    # ```

    # Note that since `type` includes one major type (`:ign`) and one minor type (`:volc`),
    # where `:volc` is a subtype of the `:ign` major type, `type` contains one minor type,
    # so `ctypes` contains one element (`ctypes.ign`), equal to 1.
    # """
    # function majorminor(type::Tuple, minorsed, minorign, minormet)
    #     # Preallocate
    #     maj = similar(collect(type))
    #     ismaj = falses(length(maj))

    #     # Assign
    #     for j in eachindex(type)
    #         maj[j] = class_up(type[j], minorsed, minorign, minormet)
    #         ismaj[j] = maj[j] == type[j]
    #     end

    #     utype = unique(maj)
    #     ctypes = NamedTuple{Tuple(utype)}([count(==(i), maj) for i in utype] .- 1)

    #     return maj, ismaj, ctypes
    # end


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
    export major_elements


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
    export standardize_units!

    dol = (12.01+2*16)/((24.869+40.08)/2+12.01+16*3)*100          # Dolomite
    gyp = (32.07+16*3+2*(18))/(40.08+32.07+16*4+2*(18))*100       # Gypsum
    bas = (32.07+16*3+0.5*(18))/(40.08+32.07+16*4+0.5*(18))*100   # Bassanite (2CaSO₄⋅H₂O)

    """
    ```julia
    screenvolatiles!(t, volatiles; 
        isevap, isgypsum, [general], [evaporite], [gypsum]
    )
    ```

    Restrict the filtered samples `t` based on the total calculated and assumed wt.% 
    `volatiles`. 

    Test all samples against `general`, evaporites against `evaporite`, and gypsum 
    against `gypsum`. Values respectively default to the wt.% volatiles in dolomite (47.6%), 
    bassanite (61.4%), and gypsum (67.4%).

    """
    function screenvolatiles!(t::BitVector, volatiles::AbstractArray; 
            isevap::BitVector, isgypsum::BitVector, 
            general::Number=dol, evaporite::Number=bas, gypsum::Number=gyp
        )

        t .&= volatiles .<= general;                        # All samples
        t[isevap] .= volatiles[isevap] .<= evaporite;       # All evaporites
        t[isgypsum] .= volatiles[isgypsum] .<= gypsum;      # Gypsum           

        return vec(t)
    end
    export screenvolatiles!


    """
    ```julia
    ca_plagioclase(K2O, Al2O3, CaO, K2O_to_Al2O3, Na2O_to_Al2O3)
    ```

    Assuming all `K₂O`, `Al₂O₃`, and `Na₂O` comes from feldspars, calculate the expected 
    CaO content of a sample. Give all values in wt.%, and specify conversion factors from 
    K₂O and Na₂O to Al₂O₃:

    ```julia
    using StatGeochem
    K2O_to_Al2O3 = (2*molarmass["Al"] + 3*molarmass["O"])/(molarmass["K"]*2 + molarmass["O"]);
    Na2O_to_Al2O3 = (2*molarmass["Al"] + 3*molarmass["O"])/(molarmass["Na"]*2 + molarmass["O"]);
    Al2O3_to_CaO = (molarmass["Ca"]+molarmass["O"])/(2*molarmass["Al"] + 3*molarmass["O"]);
    ```

    """
    function ca_plagioclase(K2O, Al2O3, Na2O, K2O_to_Al2O3, Na2O_to_Al2O3, Al2O3_to_CaO)
        # Calculate Al2O3 not from Na-plag or K-spar
        Al2O3_Ca_plag = Al2O3 - (K2O*K2O_to_Al2O3 + Na2O*Na2O_to_Al2O3)

        # Calculate CaO expected from known amount of Al2O3
        return Al2O3_Ca_plag*Al2O3_to_CaO
    end
    export ca_plagioclase


## --- Sample matching

    """
    ```julia
    likelihood(bulkage::AbstractArray, sampleage::Number,
        bulklat::AbstractArray, bulklon::AbstractArray, samplelat::Number, 
        samplelon::Number, bulkgeochem::NamedTuple, samplegeochem::NamedTuple,
        sampleidx::AbstractArray
    )
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
        ll_total = zeros(npoints, 1)

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
            ll_dist[i] = -((haversine(samplelat, samplelon, bulklat[i], bulklon[i]))^2)/(18.0^2)
        end
        
        # Geochemical log-likelihoods
        for elem in eachindex(bulkgeochem)
            @turbo for i in 1:npoints
                ll_total[i] += -((bulkgeochem[elem][i] - samplegeochem[elem].m)^2)/(samplegeochem[elem].e^2)
            end
        end

        # # Men love to scale their log-likelihoods so geochemistry and spatiotemporal 
        # # similarity have the same weight
        # n = length(bulkgeochem)
        # ll_total ./= n

        # Everyone loves addition
        ll_total .+= (ll_dist .+ ll_age)

        matched_sample = rand_prop_liklihood(ll_total)
        return sampleidx[matched_sample]
    end
    export likelihood

    """
    ```julia
    rand_prop_liklihood(ll)
    ```

    Weighted-random selection of an index based on log-likelihoods `ll`.

    See also: `weighted_rand` for weighted random selection based on likelihood / weight
    in normal space.
    """
    function rand_prop_liklihood(ll)
        log_sum_likelihoods = logsumexp(ll)
        r = rand()*exp(log_sum_likelihoods)
        s = zero(typeof(log_sum_likelihoods))
        @inbounds for i in eachindex(ll)
            s += exp(ll[i])
            if s > r
                return i
            end
        end
        return lastindex(ll)
    end
    export rand_prop_liklihood

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
    export weighted_rand

    # """
    # ```julia
    # replace_major(type::Tuple, minortypes::NamedTuple, p::NamedTuple)
    # ```

    # Replace `:sed`, `:ign`, and `:met` types in `type` with a weighted-random selection 
    # of a subtype from `minortypes` based on the weights in `p`. Remove duplicate types.

    # Weights `p` must have the keys `sed`, `ign`, and `met`.

    # ### Kwargs `ign`, `sed`, and `met`

    # Minor rock types classified under the igneous, sedimentary, and metamorphic major 
    # types.

    # # Example
    # ```julia-repl
    # julia> minortypes = (
    #         sed=(:carb, :shale), 
    #         ign=(:volc, :plut), 
    #         met=(:metased, :metaign)
    #     );

    # julia> p = (
    #         sed = [0.6, 0.4],
    #         ign = [0.5, 0.5],
    #         met = [0.3, 0.7],
    #     );

    # julia> type = (:ign, :met, :carb, :plut,);

    # julia> replace_major(type, p, minortypes)
    # (:volc, :metaign, :carb, :plut)
    # ```
    # """
    # function replace_major(type::Tuple, p::NamedTuple, minortypes::NamedTuple)
    #     type = collect(type)
    #     for i in eachindex(type)
    #         if type[i]==:ign 
    #             type[i] = minortypes.ign[weighted_rand(p.ign)]
    #         elseif type[i]==:sed
    #             type[i] = minortypes.sed[weighted_rand(p.sed)]
    #         elseif type[i]==:met 
    #             type[i] = minortypes.met[weighted_rand(p.met)]
    #         end
    #     end

    #     return Tuple(unique(type))
    # end


    # """
    # ```julia
    # get_descriptive_name(samplenames::Tuple, p_name::NamedTuple, sampletypes::Tuple,
    #         p_type::NamedTuple, typelist::NamedTuple, minortypes::NamedTuple
    #     )
    # ```

    # Given matched rock names `samplenames` and matched rock types `sampletypes`, return a 
    # list of minor types and a randomly-selected descriptive rock name.

    # # Example
    # ```julia-repl
    # name, types = get_descriptive_name(samplenames, p_name, sampletypes, p_type, typelist, minortypes)
    # ```
    # """
    # function get_descriptive_name(samplenames::Tuple, p_name::NamedTuple, sampletypes::Tuple,
    #         p_type::NamedTuple, typelist::NamedTuple, minortypes::NamedTuple
    #     )

    #     # Randomly pick a matched type
    #     i = rand(1:length(sampletypes))
    #     t = sampletypes[i]

    #     if t==:sed || t==:ign || t==:met 
    #     # If the type is a major type, replace type and pick a rock name from mapped names
    #         t = minortypes[t][weighted_rand(p_type[t])]
    #         name = typelist[t][weighted_rand(p_name[t])]

    #         # Replace any other majors, but keep the already-replaced type
    #         sampletypes = collect(sampletypes)
    #         sampletypes[i] = t
    #         sampletypes = replace_major(Tuple(sampletypes), p_type, minortypes)

    #         return Symbol(name), sampletypes
            
    #     else
    #     # If the type is a minor type, randomly pick a matched name that maps to that type
    #         q = [string(n) in typelist[t] for n in samplenames]

    #         return rand(samplenames[q]), replace_major(sampletypes, p_type, minortypes)
    #     end 
    # end
    

## --- Measurements and tuples

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
    export unmeasurementify

    """
    ```julia
    untupleify(A)
    ```

    Separate an array of `Tuple`s into two arrays.

    # Example
    ```julia-repl
    julia> A = [(1,2),(3,4),(5,6)];

    julia> untupleify(A)
    ([1, 3, 5], [2, 4, 6])
    ```
    
    """
    function untupleify(A::AbstractArray{T}) where T <: Tuple{Any, Any}
        Ta = typeof(A[1][1])
        Tb = typeof(A[1][2])
        a = Array{Ta}(undef, length(A))
        b = Array{Tb}(undef, length(A))
        for i in eachindex(A)
            a[i] = A[i][1]
            b[i] = A[i][2]
        end
        return a, b
    end
    export untupleify


## --- Geography (PIP, geolprov)

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
    export points_in_shape

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
    export coords_in_shape


    """
    ```julia
    decode_find_geolprov(geolprov::AbstractArray)
    ```

    Translate numeric codes from `find_geolprov` to province names.

    # Example
    ```julia-repl
    julia> lats = [-24.18,53.72,27.27,37.80];

    julia> lons = [-68.92,-90.90,86.45,-115.88];

    julia> provs = find_geolprov(lats, lons)
    4-element Vector{Int64}:
     12
     31
     13
     21

    julia> decode_find_geolprov(provs)
     4×1 Vector{String}:
     "Continental_Arc"
     "Shield"
     "Collisional_Orogen"
     "Rift"
    ```
    """
    function decode_find_geolprov(geolprov::AbstractVector)
        # Preallocate / define
        out = Array{String}(undef, length(geolprov), 1)
        decoder = (
            Accreted_Arc = 10,
            Island_Arc = 11,
            Continental_Arc = 12,
            Collisional_Orogen = 13,
            Extensional = 20,
            Rift = 21,
            Plume = 22,
            Shield = 31,
            Platform = 32,
            Basin = 33,
            No_Data = 00,
        )

        provnames = collect(keys(decoder))
        for k in provnames
            t = @. geolprov == decoder[k]
            out[t] .= string.(k)
        end

        return vec(out)
    end
    export decode_find_geolprov


    """
    ```julia
    decode_find_geolcont(geolcont::AbstractArray)
    ```

    Translate numeric codes from `find_geolcont` to province names.
    """
    function decode_find_geolcont(geolcont::AbstractVector)
        # Preallocate / define
        out = Array{String}(undef, length(geolcont), 1)
        decoder = (
            Africa = 1,
            Eurasia = 2,
            North_America = 3,
            South_America = 4,
            Australia = 5,
            Antarctica = 6,
            NA = 7,
        )

        contnames = collect(keys(decoder))
        for k in contnames
            t = @. geolcont == decoder[k]
            out[t] .= string.(k)
        end

        return vec(out)
    end
    export decode_find_geolcont


## --- Strings 
    """
    ```
    replace_malformed_char(str::String)
    ```

    Replace invalid or malformed characters in `str` with spaces.

    # Example 
    ```julia-repl
    julia> str = "volcanic basementb (b\xec\x81lsamo fmt.) "

    julia> replace_malformed_char(str)
    "volcanic basementb (b lsamo fmt.) "
    ```
    """
    function replace_malformed_char(str::String)
        str = collect(str)
        for i in eachindex(str)
            str[i] = ifelse(isvalid(str[i]), str[i], ' ')
        end
        return String(str)
    end
    export replace_malformed_char


## --- Resampling 

    """
    ```julia 
    resampling_age(geochem_age, geochem_age_min, geochem_age_max, 
        map_age, map_age_min, map_age_max;
        [error_default], [uncert_rel], [uncert_abs]
    ) 
    ```

    Calculate sample age and age uncertainties for a set of matched samples. Prefer ages 
    associated with geochemical age. If no ages are listed, use mapped ages. 
    
    Optionally specify a default uncertainty (`error_default` _percent_), relative 
    (`error_rel` _percent_) and / or absolute (`error_abs`, _Myr._) minimum age uncertainties. 
    If both relative and absolute uncertainties are specified, the larger uncertainty will 
    be assigned.

    # Example
    ```julia
    sampleage, ageuncert = resampling_age(mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
        macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
        uncert_rel=5, uncert_abs=50
    )
    ```

    """
    function resampling_age(geochem::T, geochem_min::T, geochem_max::T, 
            map::T, map_min::T, map_max::T;
            error_default::Number=NaN, uncert_rel::Number=NaN, uncert_abs::Number=NaN
        ) where T <: AbstractArray{<:Number}

        # Try geochemical age
        sampleage = copy(geochem);
        ageuncert = nanadd.(geochem_max, .- geochem_min) ./ 2;

        # Map age for missing values
        t = isnan.(sampleage);
        sampleage[t] .= map[t]
        ageuncert[t] .= nanadd.(map_max[t], .- map_min[t]) ./ 2;
        
        # Add minimum uncertainties, if specified
        for i in eachindex(ageuncert)
            ageuncert[i] = nanmaximum([sampleage[i]*uncert_rel/100, ageuncert[i], uncert_abs])
        end

        # Add default uncertainty, if specified 
        if !isnan(error_default)
            for i in eachindex(ageuncert)
                ageuncert[i] = ifelse(isnan(ageuncert[i]), sampleage[i]*error_default/100, ageuncert[i])
            end
        end

        return sampleage, ageuncert
    end
    export resampling_age


## --- Merge two arrays, alternating columns
    """
    mesh(A::AbstractArray{<: Number}, B::AbstractArray{<: Number})


    ```julia
    mesh(A, B)
    ```

    Combine two equally sized arrays `A` and `B` into a single array, alternating columns 
    from A and B.

    # Example
    ```julia-repl
    julia> A = ones(Int,2,2);

    julia> B = zeros(Int,2,2);

    julia> mesh(A,B)
    2×4 Matrix{Int64}:
    1  0  1  0
    1  0  1  0
    ```
    
    """
    function mesh(A::AbstractArray, B::AbstractArray)
        @assert size(A) == size(B)
        ncols, nrows = size(A[:,:])
        AB = Array{Base.promote_type(eltype(A), eltype(B))}(undef, ncols, nrows*2)
        for i in 1:nrows
            AB[:,2i - 1] .= A[:,i]
            AB[:,2i] .= B[:,i]
        end

        return AB
    end
    export mesh


## --- Numbers and statistics
    """
    ```julia 
    rescale_in_range(x, min1, max1, min2, max2)
    ```
    Rescale `x` from from a range (`min1`, `max1`) to a range (`min2`, `max2`)

    # Examples
    ```
    190.0 = rescale_in_range(1900, 0, 3800, 0, 380)

    120.0 = rescale_in_range(60, 40, 80, 0, 240)
    ```
    
    """
    rescale_in_range(x, min1, max1, min2, max2) = (max2 - min2)*(x - min1)/(max1 - min1) + min2
    export rescale_in_range 


    """
    ```julia
    normalize!(A::AbstractVector; [total])
    ```

    Normalize the values in `A` to `total` (default: 100).
    """
    function normalize!(A::AbstractVector, total::Number=100)
        sum_a = nansum(A)
        @turbo for i in eachindex(A)
            A[i] = A[i] / sum_a * total
        end
        return A
    end
    export normalize!


## --- End of file