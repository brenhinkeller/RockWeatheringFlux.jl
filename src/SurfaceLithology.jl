## --- Set up 
    # FIXME: Carbonatites coming out to NaN and not zero for continent values?
    # also check if we have at least one carbonatite sample... maybe we need more 
    # sig figs?

    # Figure out the mapped lithologies on Earth's surface and compare this to the 
    # abundance in the combined geochemical datasets  
    using RockWeatheringFlux
    using HDF5, DelimitedFiles


## --- Load data 
    # List of samples matched to a lithologic class 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Sample locations and mapped lithologic classes
    fid = h5open(macrostrat_io, "r")
    rocklat = read(fid["vars"]["rocklat"])[t]
    rocklon = read(fid["vars"]["rocklon"])[t]

    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    macro_cats = delete_cover(macro_cats)
    close(fid)

    # Bulk geochemical data 
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)

    include_minor!(bulk_cats);
    bulk_cats = delete_cover(bulk_cats)


## --- Definitions
    npoints = count(t)

    # Lithologic class of matched samples 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

    # Lithologic class definitions 
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:end];
    

## --- Fraction of land represented by each continent 
    # We're going to use this to 
    # (a) filter samples by location so we can get mapped lithology on each continent 
    # (b) convert our values from % area of each continent to % area of the world

    # Find the location (continent) of each sample and the relative contribution of each 
    # continent to total samples. Exclude NA
    continent = find_geolcont(rocklat, rocklon)
    continent = NamedTuple{Tuple([Symbol.(continents[1:end-1]); :Global])}([
        [continent .== i for i in eachindex(continents[1:end-1])]; [trues(npoints)]
    ])
    dist_cont = merge(NamedTuple{keys(continent)[1:end-1]}(normalize!(
            [count(continent[k])/npoints for k in keys(continent)[1:end-1]])
        ), (Global=100.0,)
    )


## --- Calculate the fraction of crust covered by each lithology 
    # We're saying 100% of rocks can be described as sed, ign, or undifferentiated met
    # Normalize global values so sed + ign + undiff met = 100%
    # Normalize continent values so sed + ign + undiff met = 100% continental area, which 
    # we calculated in the cell above this one. 

    # Preallocate 
    region = collect(keys(continent))
    ncols = (length(minorsed)+2) + (length(minorplut)+2) + (length(minorvolc)+2) + 3 + 4;
    result = Array{Float64}(undef, ncols, length(region));

    # Calculate surficial abundance by region
    for i in eachindex(region)
        # Filter for the region and number of samples
        s = continent[region[i]]
        nregion = count(s)

        # Calculate abundances of constituent classes in this region
        # Normalize them to 100%
        include_minor!(macro_cats);
        ign_undiff = .!(macro_cats.volc .| macro_cats.plut .| macro_cats.carbonatite);
        abundance = NamedTuple{(:sed, :volc, :plut, :carbonatite, :ign_undiff, :met_undiff)}(
            normalize!([
                count(macro_cats.sed .& s),                 # Sedimentary
                count(macro_cats.volc .& s),                # Volcanic
                count(macro_cats.plut .& s),                # Plutonic
                count(macro_cats.carbonatite .& s),         # Carbonatite
                count(macro_cats.ign .& ign_undiff .& s),   # Undifferentiated ign
                count(megaclass.met_undiff .& s)            # Undifferentiated met
            ]./nregion.*100
        ))
        abundance_ign = (abundance.volc +                   # All igneous rocks 
            abundance.plut + abundance.carbonatite + 
            abundance.ign_undiff
        )

        # The above calculated the percent abundance of each major lithologic class, 
        # normalized to 100% (i.e., 100% of the continent we're looking at)
        
        # We want to calculate the abundance of all the minor (and undifferentiated) 
        # lithologic classes, and normalize that abundance to the percentage of major 
        # classes (i.e., all sed classes should add to the total sed abundance)

        # Count the number of occurances of each minor class
        exclude_minor!(macro_cats);
        sed = float.([[count(macro_cats[k] .& s) for k in minorsed]; count(macro_cats.sed .& s)])
        volc = float.([[count(macro_cats[k] .& s) for k in minorvolc]; count(macro_cats.volc .& s)])
        plut = float.([[count(macro_cats[k] .& s) for k in minorplut]; count(macro_cats.plut .& s)])

        # Convert our counts to percentage of the total, where total is total number of
        # points in the region / continent 
        sed .= sed ./ nregion * 100
        volc .= volc ./ nregion * 100
        plut .= plut ./ nregion * 100
        
        # Normalize those minor class percentages to the abundance of the constituent class
        # This just means the percent abundances of e.g. sedimentary rocks will add to the 
        # total percent abundance of all sedimentary rocks that we calcualted earlier 
        sed .= sed ./ sum(sed) .* abundance.sed
        volc .= volc ./ sum(volc) .* abundance.volc
        plut .= plut ./ sum(plut) .* abundance.plut
        
        # Calculate the percentage of metased and metaign rocks as the fraction of total 
        # sed / ign rocks tagged as metamorphic 
        include_minor!(macro_cats);
        metased = count(megaclass.metased .& s) / count(macro_cats.sed .& s)
        metaign = count(megaclass.metaign .& s) / count(macro_cats.ign .& s)
        metased *= abundance.sed 
        metaign *= abundance_ign

        met_total = metased + metaign + abundance.met_undiff

        # Convert data from 100% of region to % of total surface area 
        area_frac = dist_cont[region[i]] / 100

        ign_out = [abundance_ign, abundance.carbonatite, abundance.ign_undiff] .* area_frac
        sed_out = [sum(sed); sed] .* area_frac
        volc_out = [sum(volc); volc] .* area_frac
        plut_out = [sum(plut); plut] .* area_frac
        met_out = [met_total, metased, metaign, abundance.met_undiff] .* area_frac

        # Save all data to the results array in the order [total; subtypes] 
        # Be SURE to check that this is in the same order as labeled in the output 🤦
        result[:,i] .= [sed_out; volc_out; plut_out; ign_out; met_out]
    end
    
    # Normalize each set of results to the global total 
    for i = 1:size(result)[1]
        normresult = result[i, 1:end-1] ./ sum(result[i, 1:end-1]) * result[i,end]
        result[i, 1:end-1] .= normresult
    end


## --- Save to file 
    # Define labels (down here so it's easier to compare to the arrays above)
    sed_label = ["Total Sedimentary"; string.(collect(minorsed)); "Undifferentiated Sedimentary"]
    volc_label = ["Total Volcanic"; string.(collect(minorvolc)); "Undifferentiated Volcanic"]
    plut_label = ["Total Plutonic"; string.(collect(minorplut)); "Undifferentiated Plutonic"]
    ign_label = ["Total Igneous", "Carbonatite", "Undifferentiated Igneous"]
    met_label = ["Total Metamorphic", "Metasedimentary", "Metaigneous", "Undifferentiated Metamorphic"]
    
    cols = ["Lithology" reshape(collect(string.(keys(continent))), 1, :)]

    # Export file, values in percentages
    labels = [sed_label; volc_label; plut_label; ign_label; met_label];
    writedlm(mapped_surface_lith, vcat(cols, hcat(labels, result)), ',')

    # Check to make sure the sums work out 
    sums_match = isapprox.(nansum(result[:, 1:end-1], dims=2), result[:, end])
    if count(sums_match) != size(result)[1]
        @warn "Mapped lithology: discrepancy in values for: $(labels[.!vec(sums_match)])" 
    end

    # Print % undifferentiated lithologies to terminal 
    target = containsi.(labels, "Undifferentiated")
    undiff_label = labels[target]
    undiff_value = round.(result[:,end][target], sigdigits=3)
    @info """ Abundance of:
    $(undiff_label[1]): $(undiff_value[1])%
    $(undiff_label[2]): \t$(undiff_value[2])%
    $(undiff_label[3]): \t$(undiff_value[3])%
    $(undiff_label[4]): \t$(undiff_value[4])%
    $(undiff_label[5]): $(undiff_value[5])%
    """


# ## --- Store major lithologic class abundances 
#     target = ("Total Sedimentary", "Total Volcanic", "Total Plutonic",
#         "Carbonatite", "Undifferentiated Igneous", "Undifferentiated Metamorphic",
#     )
#     target_i = [findfirst(x->x == target[i], labels) for i in eachindex(target)]
#     major_global_abundance = result[target_i, end]


## --- Check lithologic class abundance in geochemical dataset 
    # Undifferentiated igneous 
    ign_undiff = .!(bulk_cats.volc .| bulk_cats.plut .| bulk_cats.carbonatite);

    # Undifferentiated metamorphic can't be matched with anything else
    met_undiff = bulk_cats.met
    for k in keys(bulk_cats)    
        k == :met && continue
        met_undiff .&= .!bulk_cats[k]
    end

    # Sed, volc, plut, carbonatite, undiff ign, undiff met = 100
    n = length(bulk.Sample_ID)
    abundance = NamedTuple{(:sed, :volc, :plut, :carbonatite, :ign_undiff, :met_undiff)}(
        normalize!([
            count(bulk_cats.sed),                 # Sedimentary
            count(bulk_cats.volc),                # Volcanic
            count(bulk_cats.plut),                # Plutonic
            count(bulk_cats.carbonatite),         # Carbonatite
            count(bulk_cats.ign .& ign_undiff),   # Undifferentiated ign
            count(met_undiff)                     # Undifferentiated met
        ]./n.*100
    ))
    ign_total = (abundance.volc +                 # All igneous rocks 
        abundance.plut + abundance.carbonatite + 
        abundance.ign_undiff
    )
    abundance = merge(abundance, (; ign_total=ign_total))
    

    # Count the number of occurances of each minor class
    exclude_minor!(macro_cats);
    sed = float.([[count(macro_cats[k]) for k in minorsed]; count(macro_cats.sed)])
    volc = float.([[count(macro_cats[k]) for k in minorvolc]; count(macro_cats.volc)])
    plut = float.([[count(macro_cats[k]) for k in minorplut]; count(macro_cats.plut)])

    # Convert above counts to percentage of the total number of samples
    sed .= sed ./ n * 100
    volc .= volc ./ n * 100
    plut .= plut ./ n * 100
    
    # Calculate the percentage of metased and metaign rocks as the fraction of total 
    # sed / ign rocks tagged as metamorphic 
    include_minor!(macro_cats);
    metased = count(megaclass.metased) / count(macro_cats.sed)
    metaign = count(megaclass.metaign) / count(macro_cats.ign)
    metased *= abundance.sed 
    metaign *= abundance.ign_total

    met_total = metased + metaign + abundance.met_undiff

    # Normalize those minor class percentages to the abundance of the constituent class
    # This just means the percent abundances of e.g. sedimentary rocks will add to the 
    # total percent abundance of all sedimentary rocks that we calcualted earlier 
    sed .= sed ./ sum(sed) .* abundance.sed
    volc .= volc ./ sum(volc) .* abundance.volc
    plut .= plut ./ sum(plut) .* abundance.plut

    # Slam that bad boy together
    # We don't need to do any of the area shit because we're not separating by continent here
    ign_out = [abundance.ign_total, abundance.carbonatite, abundance.ign_undiff]
    sed_out = [sum(sed); sed]
    volc_out = [sum(volc); volc]
    plut_out = [sum(plut); plut]
    met_out = [met_total, metased, metaign, abundance.met_undiff]

    # Save all data to the results array in the order [total; subtypes] 
    # Be SURE to check that this is in the same order as labeled in the output 🤦
    result = [sed_out; volc_out; plut_out; ign_out; met_out]

    # Wheeeeeeeeeeeeeeeee! We could export this if we wanted 
    bulkabundance = vcat(["Lithology" "Abundance"], hcat(labels, result))


## --- Compare the two abundances 
    # Load the global file
    surfabundance = readdlm(mapped_surface_lith, ',')[:,end]

    # See how far off we are 
    difference = bulkabundance[2:end,2] .- surfabundance[2:end]

    smashed = (hcat(bulkabundance[2:end,2], surfabundance[2:end], difference))
    smashed = round.(smashed, sigdigits=3)
    smashed = (hcat(bulkabundance[:,1], vcat(["Geochem" "Maps" "Difference"], smashed)))
    display(smashed)
    for i in eachindex(smashed[:,1])
        println(join(smashed[i,:], ";"))
    end



## --- Spatially resampled geochemical dataset 
#     # Convert lithology to a resample-able data array
#     a, head = cats_to_array(bulk_cats)
#     a_err = zeros(size(a))

#     # Calculate spatial weights and resample
#     nsims = 10_000
#     k = Array{Float64}(undef, count(matched))
#     try 
#         k .= readdlm("output/bulk_k.csv", ',')
#     catch
#         k .= invweight_location(bulk.Latitude[matched], bulk.Longitude[matched])
#         writedlm("output/bulk_k.csv", k, ',')
#     end
#     p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
#     simout = bsresample(a, a_err, nsims, vec(p))

#     # Parse data back into a Tuple
#     data = @. simout > 0
#     sim_cats = NamedTuple{Tuple(Symbol.(head))}([data[:,i] for i in eachindex(head)])
#     include_minor!(sim_cats)

#     # Do as above 
#     ign_undiff = .!(sim_cats.volc .| sim_cats.plut .| sim_cats.carbonatite);
#     abundance = NamedTuple{(:sed, :volc, :plut, :carbonatite, :ign_undiff, :met_undiff)}(
#         normalize!([
#             count(sim_cats.sed),                 # Sedimentary
#             count(sim_cats.volc),                # Volcanic
#             count(sim_cats.plut),                # Plutonic
#             count(sim_cats.carbonatite),         # Carbonatite
#             count(sim_cats.ign .& ign_undiff),   # Undifferentiated ign
#             count(sim_cats.met)                  # Undifferentiated met
#         ]./npoints.*100
#     ))
#     abundance_ign = (abundance.volc +                   # All igneous rocks 
#         abundance.plut + abundance.carbonatite + 
#         abundance.ign_undiff
#     )

#     # Store major lithologic class abundances 
#     major_resampled_abundance = [abundance.sed, abundance.volc, abundance.plut, 
#         abundance.carbonatite, abundance.ign_undiff, abundance.met_undiff]


# ## --- Print to terminal for LaTeX 
#     cols = ["" "Geochemical Data" "Resampled" "Mapped"]
#     table = round.(
#         hcat(major_dataset_abundance, major_resampled_abundance, major_global_abundance),
#         digits=1
#     )
#     table = vcat(cols, hcat(collect(target), table))
#     for i in 1:size(table)[1]
#         println("$(join(table[i,:], " & ")) \\ \b\\")
#     end


## --- End of file 