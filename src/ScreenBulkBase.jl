## --- Load and parse data file
    bulk = matread("data/bulk.mat")["bulk"];
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk));


## --- Load and parse metadata
    bulktext = matread("data/bulktext.mat")["bulktext"]

    # Type stabilize sample metadata
    elements = unique([String.(bulktext["elements"]); ["Units", "Methods"]])
    for n in elements
        bulktext[n] = String.(bulktext[n])
    end

    # Densify sparse arrays, correct for zero-indexing and Float64 type
    for k in keys(bulktext["unit"])
        bulktext["unit"][k] = collect(bulktext["unit"][k])
        bulktext["unit"][k] .+= 1.0
        bulktext["unit"][k] = Int.(bulktext["unit"][k])
    end
    for k in keys(bulktext["method"])
        bulktext["method"][k] = collect(bulktext["method"][k])
        bulktext["method"][k] .+= 1.0
        bulktext["method"][k] = Int.(bulktext["method"][k])
    end
    for k in keys(bulktext["index"])
        bulktext["index"][k] = collect(bulktext["index"][k])
        bulktext["index"][k] .+= 1.0
        bulktext["index"][k] = Int.(bulktext["index"][k])
    end

    # Information by element
    bulktext_units = (
        units = bulktext["Units"],
        index = NamedTuple{Tuple(Symbol.(collect(keys(bulktext["unit"]))))}(
            [bulktext["unit"][i] for i in keys(bulktext["unit"])]
        )
    )
    bulktext_methods = (
        methods = bulktext["Methods"],
        index = NamedTuple{Tuple(Symbol.(collect(keys(bulktext["method"]))))}(
            [bulktext["method"][i] for i in keys(bulktext["method"])]
        )
    )

    # Information by sample
    bulktext = (
        elements = NamedTuple{Tuple(Symbol.(bulktext["elements"]))}(
            [bulktext[i] for i in bulktext["elements"]]
        ),
        index = NamedTuple{Tuple(Symbol.(collect(keys((bulktext["index"])))))}(
            [bulktext["index"][i] for i in keys(bulktext["index"])]
        )
    )


## --- Get rock type matches for all samples: DIY and save to a file
    # # Rock names / types / materials for all EarthChem data
    # bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name])
    # bulkrocktype = lowercase.(bulktext.elements.Type[bulktext.index.Type])
    # bulkmaterial = lowercase.(bulktext.elements.Material[bulktext.index.Material])

    # typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();

    # # Match rock types
    # bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial, 
    #     (minorsed..., :sed,), (minorvolc..., minorplut..., minorign..., :ign)
    # )

    # # Save to file
    # fid = h5open("output/bulk_unrestricted_types.h5", "w")
    # a = Array{Int64}(undef, length(bulk_cats[1]), length(bulk_cats))
    # for i in eachindex(keys(bulk_cats))
    #     for j in eachindex(bulk_cats[i])
    #         a[j,i] = ifelse(bulk_cats[i][j], 1, 0)
    #     end
    # end
    # fid["bulk_cats"] = a
    # fid["bulk_cats_head"] = string.(collect(keys(bulk_cats)))
    # close(fid)


## --- Or save yourself the time and load from a file
    fid = h5open("output/bulk_unrestricted_types.h5", "r")
    header = read(fid["bulk_cats_head"])
    data = read(fid["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # We want major types to include minor types
    # Remember there are no minor metamorphic rocks anymore! Goodbye! I'll miss you :(
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorsed
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorvolc
        bulk_cats.volc .|= bulk_cats[type]
    end
    for type in minorplut
        bulk_cats.plut .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end


## --- Define and organize elements to convert units
    # Define keys for the output file. Some already exist, volatile information needs added
    majors, minors = get_elements()
    allelements = [majors; minors]
    deleteat!(allelements, findall(x->x==:Volatiles,allelements))

    existing = [:Latitude, :Longitude, :Loc_Prec, :Age, :Age_Max, :Age_Min]
    new = [:Volatiles, :Volatiles_Reported, :Volatiles_Calc, :Volatiles_Assumed]
    allkeys = [allelements; new; existing]

    # LOI, CaCO3, H2O, and CO2 don't get reported in the output, but we'll want them with 
    # the correct units for conversions 
    source_elements = [allelements; [:Loi, :CaCO3, :H2O, :CO2]]

    # Get all units present in bulktext
    presentunits = join(collect(keys(bulktext_units.index)), " ")

    # Define densities for relevant major element oxides (g/cm³)
    density = (Al2O3=3.95, K2O=2.35, MgO=3.58, Na2O=2.27, P2O5=2.39, SiO2=2.65, TiO2=4.23, 
        CaO=3.34, FeOT=5.74
    )


## --- Convert all units to wt.%
    # This method is resistant to changes in element order
    for i in source_elements
        # Check unit is present, and convert indices to a list of units
        unit_i = string(i) * "_Unit"
        @assert contains(presentunits, unit_i) "$unit_i not present"
        bulkunits = bulktext_units.units[bulktext_units.index[Symbol(unit_i)]]

        # Convert units of all sample data
        standardize_units!(bulk[i], bulkunits, density, elem=i)
    end


## --- Prepare to calculate volatiles and note all steps 
    # Preallocate
    volatiles_reported = Array{Float64}(undef, length(bulk.SiO2), 1)    # [CO₂ + H₂O] OR LOI
    volatiles_calc = Array{Float64}(undef, length(bulk.SiO2), 1)        # E.g. carbonate CO2
    volatiles_assumed = Array{Float64}(undef, length(bulk.SiO2), 1)     # 100 - bulk weight
    volatiles = Array{Float64}(undef, length(bulk.SiO2), 1)             # Calc + assumed

    # Compute some conversion factors
    CaCO3_to_CaO = (40.078+15.999)/(40.078+15.999*3+12.01)
    CaCO3_to_CO2 = (15.999*2+12.01)/(40.078+15.999*3+12.01)

    CaO_to_CO2 = (15.999*2+12.01)/(40.078+15.999)               # Assumes CaCO3
    CaO_to_H2O = 2*(2*1.00784+15.999)/(40.078+15.999)           # Gypsum
    CaO_to_SO3 = (32.065+3*15.999)/(40.078+15.999)              # Gypsum

    K2O_to_Al2O3 = (26.98*2+15.999*3)/(39.098*2+15.999)       # K-spar
    Na2O_to_Al2O3 = (26.98*2+15.999*3)/(22.990*2+15.999)      # Na-plag
    Al2O3_to_CaO = (40.078+15.999)/(26.98*2+15.999*3)         # Ca-plag


## --- Calculate reported volatiles as the larger of [CO₂ + H₂O] OR LOI
    @turbo volatiles_reported .= nanadd.(bulk.CO2, bulk.H2O)

    @turbo for i in eachindex(volatiles)
        volatiles_reported[i] = ifelse(bulk.Loi[i] > volatiles_reported[i], 
            bulk.Loi[i], volatiles_reported[i]
        )
    end


## --- Calculate CO2 and SO3 from CaO for carbonates and evaporites 
    # There are a lot of limestones with only CaO :(

    # All samples: convert CaCO₃ to CaO + CO₂, and get CO₂ from CaO
    converted = falses(length(bulk.CaCO3));
    @turbo for i in eachindex(bulk.CaCO3)
        # Replace the existing CaO value with the computed CaO value
        bulk.CaO[i] = ifelse(bulk.CaCO3[i] > 0, CaCO3_to_CaO * bulk.CaCO3[i], bulk.CaO[i])

        # Add calculated CO2 to the reported values
        bulk.CO2[i] = ifelse(bulk.CaCO3[i] > 0, 
            nanadd(CaCO3_to_CO2 * bulk.CaCO3[i], bulk.CO2[i]), 
            bulk.CO2[i]
        )
        converted[i] = (bulk.CaCO3[i] > 0)
    end

    # Carbonates: calculate CO₂ from CaO
    # Replace existing values, in case any samples already did this
    target = bulk_cats.carb .| bulk_cats.siliciclast .| bulk_cats.shale;
    @turbo for i in eachindex(bulk.CaO)
        bulk.CO2[i] = ifelse(target[i] & !converted[i],     # Unless we did that before...
            bulk.CaO[i]*CaO_to_CO2,                         # Replace value with calculated CO₂
            bulk.CO2[i]                                     # Otherwise, use reported CO₂
        )
    end

    # Siliciclasts and shales: calculate CO₂ from CaO not associated with Ca-plag
    CaO_excess = similar(bulk.CO2)
    @turbo for i in eachindex(CaO_excess)
        CaO_plag = ca_plagioclase(bulk.K2O[i], bulk.Al2O3[i], bulk.Na2O[i], 
            K2O_to_Al2O3, Na2O_to_Al2O3, Al2O3_to_CaO
        )
        CaO_excess[i] = bulk.CaO[i] - CaO_plag
    end
    
    target = @. (bulk_cats.siliciclast .| bulk_cats.shale) & !converted
    @turbo for i in eachindex(CO2)
        bulk.CO2[i] = ifelse(target[i] & (CaO_excess[i] > 0), CaO_excess[i]*CaO_to_CO2, bulk.CO2[i])
    end


    # Gypsum: calculate H₂O and SO₃ from CaO 
    SO3 = Array{Float64}(undef, length(bulk.SiO2), 1)
    isgypsum = bulktext.elements.Rock_Name[bulktext.index.Rock_Name] .== "GYPSUM"
    @turbo for i in eachindex(bulk.H2O)
        bulk.H2O[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_H2O, bulk.H2O[i])
        SO3[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_SO3, NaN)
    end

    # Compute volatiles
    # Initialize array with the larger of (calculated!) [CO₂ + H₂O + SO₃] OR LOI
    novalue = @. isnan(bulk.H2O) & isnan(bulk.CO2) & isnan(SO3);
    zeronan!(bulk.CO2)
    zeronan!(bulk.H2O)
    zeronan!(SO3)

    @turbo volatiles .= bulk.CO2 .+ bulk.H2O .+ SO3
    @turbo for i in eachindex(volatiles)
        volatiles[i] = ifelse(bulk.Loi[i] > volatiles[i], bulk.Loi[i], volatiles[i])
    end

    # Save the calculated volatiles
    volatiles_calc .= volatiles
    volatiles_calc[novalue] .= NaN


## --- Compute total wt.% analyzed for all samples
    # Exclude OIBs and rocks below the shelf break
    isOIB = findOIBs(bulk.Latitude, bulk.Longitude);

    # etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    # elev = find_etopoelev(etopo, bulk.Latitude, bulk.Longitude)
    # abovesea = elev .> -140;            # Shelf break at 140 m
    abovesea = bulk.Elevation .> -140   # Equivalent results, but faster

    # Compute bulk analyzed weight for rocks above sea level
    bulkweight = Array{Float64}(undef, length(bulk.SiO2), 1)

    p = Progress(length(bulkweight) ÷ 10, desc="Calculating wt.% ...", enabled=show_progress)
    @inbounds for i in eachindex(bulkweight)
        if abovesea[i] & !isOIB[i]
            bulkweight[i] = nansum([bulk[j][i] for j in allelements]) + volatiles[i]
        else
            bulkweight[i] = NaN
        end
        i%10==0 && next!(p)
    end


## --- Correct for assumed unmeasured volatile loss in sedimentary rocks
    # If the sample is a sedimentary rock with a total analyzed wt.% below 100%, assume 
    # the "missing" data are volatiles that were not included in the analysis
    # additional = zeros(length(bulkweight))
    for i in eachindex(bulkweight)
        volatiles_assumed[i] = ifelse(bulk_cats.sed[i] && bulkweight[i] < 100, 
            (100 - bulkweight[i]), 0.0)
    end

    # How many samples before assuming additional volatiles?
    t = @. 84 <= bulkweight <= 104
    tᵢ = count(t)

    # Get all samples that are within acceptable range after assuming volatiles
    t = @. 84 <= bulkweight .+ volatiles_assumed <= 104;

    # Add calculated and assumed volatiles for the total volatiles value
    @turbo volatiles .+= volatiles_assumed

    # Restrict!
    t₂ = vec(copy(t))
    t = vec(t)
    screenvolatiles!(t, volatiles, isevap=bulk_cats.evap, isgypsum=vec(isgypsum),
        general=dol, evaporite=bas, gypsum=gyp    
    )

    # temp test to make sure we're actually modifying the value
    @assert t != t₂

    # Print to terminal
    nsamples = round(count(t)/length(t)*100, digits=2)
    up = count(t) - tᵢ
    @info """Saving $(count(t)) samples ($nsamples%)
    Assuming volatiles increased count from $tᵢ to $(count(t))
    Total increase = $up
    """


## --- Restrict to in-bounds only and normalize compositions
    # This is inefficient, but is not sensitive to changes in the order of any elements
    newbulk = merge(bulk, (
        Volatiles=volatiles, 
        Volatiles_Reported=volatiles_reported, 
        Volatiles_Calc=volatiles_calc, 
        Volatiles_Assumed=volatiles_assumed,
        Sample_ID = collect(1:length(bulk.SiO2))
    ))
    bulk = NamedTuple{Tuple(allkeys)}([newbulk[i][t] for i in allkeys])

    # Normalize compositions to 100%
    contributing = [allelements; :Volatiles]             # Need to re-include volatiles!
    for i in eachindex(bulk.SiO2)
        sample = [bulk[j][i] for j in contributing]      # Get it
        normalize!(sample)                               # Normalize it
        for j in eachindex(contributing)                 # Put it back
            bulk[contributing[j]][i] = sample[j]
        end
    end


## --- Restrict bulktext to data we're keeping in the final output file
    # Methods
    restrictmethods = Symbol.(string.(allelements) .* "_Meth")
    bulktext_methods = (
        methods = bulktext_methods.methods,
        index = NamedTuple{Tuple(restrictmethods)}(
            [bulktext_methods.index[i][t] for i in restrictmethods]
        )
    )

    # Other / general elements
    bulktext = (
        elements = bulktext.elements,
        index = NamedTuple{Tuple(keys(bulktext.index))}(
            [bulktext.index[i][t] for i in keys(bulktext.index)]
        )
    )


## -- Restrict rock types
    # Restrict
    bulk_kittens = NamedTuple{keys(bulk_cats)}(bulk_cats[k][t] for k in keys(bulk_cats))

    # Check that all rock types have samples
    for k in keys(bulk_kittens)
        if count(bulk_kittens[k]) < 10
            @warn "After filtering, type \"$k\" has $(count(bulk_kittens[k])) samples."
        end
    end


## --- Write data to an HDF5 file
    fid = h5open(fileout, "w")
    data = create_group(fid, "bulk")
    text = create_group(fid, "bulktext")

    # Bulk
    write(data, "header", string.(allkeys))
    writebulk = create_dataset(data, "data", Float64, (count(t), length(allkeys)))
    for i in eachindex(allkeys)
        writebulk[:,i] = bulk[allkeys[i]]
    end

    # Bulktext methods
    methods = create_group(text, "methods")
    write(methods, "methods", bulktext_methods.methods)
    write(methods, "header", string.(restrictmethods))
    index = create_dataset(methods, "index", Int64, (count(t), length(restrictmethods)))
    for i in eachindex(restrictmethods)
        index[:,i] = bulktext_methods.index[restrictmethods[i]]
    end

    # Bulk text sample data
    sampledata = create_group(text, "sampledata")
    elements = create_group(sampledata, "elements")
    for i in keys(bulktext.elements)
        write(elements, string(i), bulktext.elements[i])
    end

    textkeys = collect(keys(bulktext.index))
    write(sampledata, "header", string.(textkeys))
    index = create_dataset(sampledata, "index", Int64, (count(t), length(textkeys)))
    for i in eachindex(textkeys)
        index[:,i] = bulktext.index[textkeys[i]]
    end

    # Rock types and rock names
    bulktypes = create_group(fid, "bulktypes")

    # Rock types
    a = Array{Int64}(undef, length(bulk_kittens[1]), length(bulk_kittens))
    for i in eachindex(keys(bulk_kittens))
        for j in eachindex(bulk_kittens[i])
            a[j,i] = ifelse(bulk_kittens[i][j], 1, 0)
        end
    end
    bulktypes["bulk_cats"] = a
    bulktypes["bulk_cats_head"] = string.(collect(keys(bulk_kittens))) 

    close(fid)

    @info "Saved new bulk.h5 file to $fileout."


## --- End of File