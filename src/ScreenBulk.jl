# This file will read bulk.mat and bulktext.mat and standardize all units as wt.%, save
# data as a new file.

# TO DO: 
    # Re-read in old files, and figure out what to do for parsing / correcting element 
        # oxide weights. For example, some samples may have different values for CaO and 
        # CaCO3, but Mg and MgO may represent the same data.
    # Read in raw data files. Useful code:
        # # Filter ages younger than 0 or greater than the age of the earth
        # invalid_age = vcat(findall(>(4000), bulk.Age), findall(<(0), bulk.Age))
        # bulk.Age[invalid_age] .= NaN

        # # Fill in any missing ages from bounds
        # for i in eachindex(bulk.Age)
        #     if isnan(bulk.Age[i])
        #         bulk.Age[i] = nanmean([bulk.Age_Max[i], bulk.Age_Min[i]])
        #     end
        # end

## --- Set up
    # Packages
    using MAT
    using StatGeochem
    using LoopVectorization
    using HDF5
    using Static
    using Measurements
    using ProgressMeter

    # Local utilities
    include("utilities/Utilities.jl")


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

    # # List of rock names
    # rocknames = get_rock_class(major=true, inclusive=false)
    # rocknames = unique((rocknames.sed..., rocknames.met..., rocknames.ign...))

    # # Rock subtypes, major types do not include minors (for class_up function call)
    # typelist = get_rock_class(major=false, inclusive=false)

    # # Match rock types
    # bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial; unmultimatch=false, 
    #     inclusive=false, source=:earthchem
    # )

    # # Match rock names
    # p = Progress(length(rocknames), desc="Finding EarthChem samples for each rock name")
    # bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([falses(length(bulkrockname)) 
    #     for _ in eachindex(rocknames)]
    # )
    # for i in eachindex(rocknames)
    #     bulk_lookup[i] .= find_earthchem(rocknames[i], bulkrockname, bulkrocktype, 
    #         bulkmaterial
    #     )

    #     # If no matches, jump up a class. Find everything within that class
    #     if count(bulk_lookup[i]) == 0
    #         searchlist = typelist[class_up(typelist, rocknames[i])]

    #         # Search all of those names; each class should at least have something
    #         for j in eachindex(searchlist)
    #             bulk_lookup[i] .|= find_earthchem(searchlist[j], bulkrockname, bulkrocktype, 
    #                 bulkmaterial
    #             )
    #         end
    #     end
    #     next!(p)
    # end

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

    # a = Array{Int64}(undef, length(bulk_lookup[1]), length(bulk_lookup))
    # for i in eachindex(keys(bulk_lookup))
    #     for j in eachindex(bulk_lookup[i])
    #         a[j,i] = ifelse(bulk_lookup[i][j], 1, 0)
    #     end
    # end
    # fid["bulk_lookup"] = a
    # fid["bulk_lookup_head"] = string.(collect(keys(bulk_lookup)))
    # close(fid)


## --- Or save yourself the time and load from a file
    fid = h5open("output/bulk_unrestricted_types.h5", "r")
    header = read(fid["bulk_cats_head"])
    data = read(fid["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    
    # We want major elements to include major elements
    minorsed, minorign, minormet = get_minor_types()
    for type in minorsed
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end
    for type in minormet
        bulk_cats.met .|= bulk_cats[type]
    end

    # And load rock names
    rocknames = read(fid["bulk_lookup_head"])
    data = read(fid["bulk_lookup"])
    data = @. data > 0
    bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}(
        [data[:,i] for i in eachindex(rocknames)])
    close(fid)


## --- Define and organize elements to convert units
    # Get major and minor elements, but remove :Volatiles from before it causes problems
    majors, minors = get_elements()
    allelements = [majors; minors]
    allkeys = [allelements; [:Latitude, :Longitude, :Age]]
    deleteat!(allelements, findall(x->x==:Volatiles,allelements))

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


## --- Fix some missing data and weird conversions
    # Pre-compute some conversion factors
    CaCO3_to_CaO = (40.078+15.999)/(40.078+15.999*3+12.01)
    CaCO3_to_CO2 = (15.999*2+12.01)/(40.078+15.999*3+12.01)

    CaO_to_CO2 = (15.999*2+12.01)/(40.078+15.999)
    CaO_to_H2O = 2*(2*1.00784+15.999)/(40.078+15.999)
    CaO_to_SO3 = (32.065+3*15.999)/(40.078+15.999)

    # CaCO₃ → CaO + CO₂, and get CO₂ from CaO
    # Replace the existing CaO value with the computed CaO value
    # Add the existing CO2 value to the computed CO2 values 
    @turbo for i in eachindex(bulk.CaCO3)
        bulk.CaO[i] = ifelse(bulk.CaCO3[i] > 0, CaCO3_to_CaO * bulk.CaCO3[i], bulk.CaO[i])
        bulk.CO2[i] = ifelse(bulk.CaCO3[i] > 0, CaCO3_to_CO2 * bulk.CaCO3[i] + bulk.CO2[i],
            bulk.CO2[i]
        )
    end

    # CARBONATES ONLY: Calculate CO₂ from CaO
    # There are a lot of limestones with only CaO :(

    # What if we also do this for siliciclasts and shales?
    target = bulk_cats.carb .| bulk_cats.siliciclast .| bulk_cats.shale;
    @turbo for i in eachindex(bulk.CaO)
        bulk.CO2[i] += ifelse(target[i], bulk.CaO[i]*CaO_to_CO2, NaN)
    end

    # GYPSUM ONLY: CaO → H₂O, SO₄ (technically SO₃ because already have an oxygen in CaO)
    SO3 = Array{Float64}(undef, length(bulk.SiO2), 1)
    isgypsum = bulktext.elements.Rock_Name[bulktext.index.Rock_Name] .== "GYPSUM"
    @turbo for i in eachindex(bulk.H2O)
        bulk.H2O[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_H2O, NaN)
        SO3[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_SO3, 0)
    end


## --- Create new volatiles measurement
    # Initialize as the larger of (CO₂ + H₂O) OR LOI, and add calculated SO₃
    volatiles = Array{Float64}(undef, length(bulk.SiO2), 1)

    zeronan!(bulk.CO2)
    zeronan!(bulk.H2O)
    @turbo volatiles .= bulk.CO2 .+ bulk.H2O

    @turbo for i in eachindex(volatiles)
        volatiles[i] = ifelse(bulk.Loi[i] > volatiles[i], bulk.Loi[i], volatiles[i])
    end

    @turbo volatiles .+= SO3


## --- Compute total wt.% analyzed for all samples
    # Exclude rocks below the shelf break
    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, bulk.Latitude, bulk.Longitude)
    abovesea = elev .> -140;    # Shelf break at 140 m

    # Compute bulk analyzed weight for rocks above sea level
    bulkweight = Array{Float64}(undef, length(bulk.SiO2), 1)

    p = Progress(length(bulkweight) ÷ 10, desc="Calculating wt.% ...")
    @inbounds for i in eachindex(bulkweight)
        if abovesea[i]
            bulkweight[i] = nansum([bulk[j][i] for j in allelements]) + volatiles[i]
        else
            bulkweight[i] = NaN
        end
        i%10==0 && next!(p)
    end


## --- Correct for assumed unmeasured volatile loss in sedimentary rocks
    # If the sample is a sedimentary rock with a total analyzed wt.% below 84%, assume 
    # the "missing" data is volatiles that were not included in the analysis
    additional = zeros(length(bulkweight))
    for i in eachindex(bulkweight)
        additional[i] = ifelse(bulk_cats.sed[i] && bulkweight[i] < 100, 
            (100 - bulkweight[i]), 0.0)
    end

    # How many samples before assuming additional volatiles?
    t = @. 84 <= bulkweight <= 104
    tᵢ = count(t)

    # Don't let samples through if the assumed volatile is more than 50 wt.%
    t = @. 84 <= bulkweight .+ additional <= 104
    t .&= additional .< 50
    t = vec(t)

    volatiles .+= additional

    # Print to terminal
    nsamples = round(count(t)/length(t)*100, digits=2)
    up = count(t) - tᵢ
    @info """Saving $(count(t)) samples ($nsamples%)
    Assuming volatiles increased count from $tᵢ to $(count(t))
    Total increase = $up
    """

## --- Save an intermediate file for analysis
    fid = h5open("output/itermediate_screen.h5", "w")
    g = create_group(fid, "vars")
        g["SiO2"] = bulk.SiO2
        g["bulkweight"] = bulkweight
        g["additional"] = additional
    close(fid)


## --- Restrict to in-bounds only and normalize compositions
    # This is inefficient, but is not sensitive to changes in the order of any elements
    newbulk = merge(bulk, (Volatiles=volatiles,))
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
    bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name])

    # Other / general elements
    bulktext = (
        elements = bulktext.elements,
        index = NamedTuple{Tuple(keys(bulktext.index))}(
            [bulktext.index[i][t] for i in keys(bulktext.index)]
        )
    )


## --- Classify rock names (getting from file is too tricky to debug)
    # Rock names / types / materials for all EarthChem data
    bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name])
    bulkrocktype = lowercase.(bulktext.elements.Type[bulktext.index.Type])
    bulkmaterial = lowercase.(bulktext.elements.Material[bulktext.index.Material])

    # List of rock names
    rocknames = get_rock_class(major=true, inclusive=false)
    rocknames = unique((rocknames.sed..., rocknames.met..., rocknames.ign...))

    # Rock subtypes, major types do not include minors (for class_up function call)
    typelist = get_rock_class(major=false, inclusive=false)

    # Match rock types
    bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial, 
        unmultimatch=false, inclusive=false, source=:earthchem
    )

    # Match rock names
    p = Progress(length(rocknames), desc="Finding EarthChem samples for each rock name")
    bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([falses(length(bulkrockname)) 
        for _ in eachindex(rocknames)]
    )
    for i in eachindex(rocknames)
        bulk_lookup[i] .= find_earthchem(rocknames[i], bulkrockname, bulkrocktype, 
            bulkmaterial
        )

        # If no matches, jump up a class. Find everything within that class
        if count(bulk_lookup[i]) == 0
            searchlist = typelist[class_up(typelist, rocknames[i])]

            # Search all of those names; each class should at least have something
            for j in eachindex(searchlist)
                bulk_lookup[i] .|= find_earthchem(searchlist[j], bulkrockname, bulkrocktype, 
                    bulkmaterial
                )
            end
        end
        next!(p)
    end


## --- Write data to an HDF5 file
    @info "Saving new bulk.h5 file."

    fid = h5open("output/bulk.h5", "w")
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
    a = Array{Int64}(undef, length(bulk_cats[1]), length(bulk_cats))
    for i in eachindex(keys(bulk_cats))
        for j in eachindex(bulk_cats[i])
            a[j,i] = ifelse(bulk_cats[i][j], 1, 0)
        end
    end
    bulktypes["bulk_cats"] = a
    bulktypes["bulk_cats_head"] = string.(collect(keys(bulk_cats))) 

    # Rock names
    a = Array{Int64}(undef, length(bulk_lookup[1]), length(bulk_lookup))
    for i in eachindex(keys(bulk_lookup))
        for j in eachindex(bulk_lookup[i])
            a[j,i] = ifelse(bulk_lookup[i][j], 1, 0)
        end
    end
    bulktypes["bulk_lookup"] = a
    bulktypes["bulk_lookup_head"] = string.(collect(keys(bulk_lookup)))

    close(fid)


## --- End of File