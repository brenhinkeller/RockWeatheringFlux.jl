# This file will read bulk.mat and bulktext.mat and standardize all units as wt.%, save
# data as a new file.

# TO DO: 
    # Re-read in old files, and figure out what to do for parsing / correcting element 
        # oxide weights. For example, some samples may have different values for CaO and 
        # CaCO3, but Mg and MgO represent the same data.
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
    
    # Load and parse data files
    bulk = matread("data/bulk.mat")["bulk"]
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))


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


## --- Define and organize elements to correct
    majors, minors = get_elements()

    # Collect major and minor elements together. Make sure to keep rock type!
    allelements = [majors; minors; :H2O; :CO2; :Loi]
    allkeys = [allelements; [:Type, :Latitude, :Longitude, :Age]]
    ndata = length(allelements)

    # Get all units present in bulktext
    presentunits = join(collect(keys(bulktext_units.index)), " ")

    # Define densities for relevant major element oxides (g/cm³)
    density = (Al2O3=3.95, K2O=2.35, MgO=3.58, Na2O=2.27, P2O5=2.39, SiO2=2.65, TiO2=4.23, 
        CaO=3.34, FeOT=5.74
    )


## --- Convert all units to wt.%
    # This method is resistant to changes in element order
    for i in allelements
        # Check unit is present, and convert indices to a list of units
        unit_i = string(i) * "_Unit"
        @assert contains(presentunits, unit_i) "$unit_i not present"
        bulkunits = bulktext_units.units[bulktext_units.index[Symbol(unit_i)]]

        # Convert units of all sample data
        standardize_units!(bulk[i], bulkunits, density, elem=i)
    end


## --- Restrict bulk to whole rock analysis of continental crust
    # Samples are 84-104 wt.% analyzed
    bulkweight = zeros(length(bulk.SiO2))
    t = falses(length(bulkweight))
    bounds = (84, 104)
    @time for i in eachindex(bulkweight)
        bulkweight[i] = nansum([bulk[j][i] for j in allelements])
        t[i] = ifelse(bounds[1] < bulkweight[i] < bounds[2], true, false)
    end
    
    # Samples are on land
    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, bulk.Latitude, bulk.Longitude)
    abovesea = elev .> 0
    t .&= abovesea


## --- Terminal print, normalization, restriction
    nsamples = round(count(t)/length(t)*100, digits=2)
    @info "Saving $nsamples% samples between $(bounds[1])% and $(bounds[2])% analyzed wt.%"

    # Restrict bulk to in-bounds only
    bulk = NamedTuple{Tuple(allkeys)}([bulk[i][t] for i in allkeys])

    # Normalize compositions to 100%
    for i in eachindex(bulk.SiO2)
        sample = [bulk[j][i] for j in allelements]      # Get it
        normalize!(sample)                              # Normalize it
        for j in eachindex(allelements)                 # Put it back
            bulk[allelements[j]][i] = sample[j]
        end
    end

## --- Restrict bulktext to in-bounds only data for elements of interest
    restrictunits = Symbol.(string.(allelements) .* "_Unit")
    restrictmethods = Symbol.(string.(allelements) .* "_Meth")

    bulktext_units = (
        units = bulktext_units.units,
        index = NamedTuple{Tuple(restrictunits)}(
            [bulktext_units.index[i][t] for i in restrictunits]
        )
    )
    bulktext_methods = (
        methods = bulktext_methods.methods,
        index = NamedTuple{Tuple(restrictmethods)}(
            [bulktext_methods.index[i][t] for i in restrictmethods]
        )
    )
    bulktext = (
        elements = bulktext.elements,
        index = NamedTuple{Tuple(keys(bulktext.index))}(
            [bulktext.index[i][t] for i in keys(bulktext.index)]
        )
    )


## --- Write data to an HDF5 file
    fid = h5open("output/bulk.h5", "w")
    data = create_group(fid, "bulk")
    text = create_group(fid, "bulktext")

    # Bulk
    write(data, "header", string.(allkeys))
    writebulk = create_dataset(data, "data", Float64, (count(t), length(allkeys)))
    for i in eachindex(allkeys)
        writebulk[:,i] = bulk[allkeys[i]]
    end

    # Bulktext units
    units = create_group(text, "units")
    write(units, "units", bulktext_units.units)
    write(units, "header", string.(restrictunits))
    index = create_dataset(units, "index", Int64, (count(t), length(restrictunits)))
    for i in eachindex(restrictunits)
        index[:,i] = bulktext_units.index[restrictunits[i]]
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


## --- Find rock type and rock name matches
    # Interpret index / element formation of bulktext
    bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name])
    bulkrocktype = lowercase.(bulktext.elements.Type[bulktext.index.Type])
    bulkmaterial = lowercase.(bulktext.elements.Material[bulktext.index.Material])

    # Rock types
    bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial; unmultimatch=false, 
        inclusive=false, source=:earthchem
    )
    
    # Rock names
    rocknames = get_rock_class(true)
    rocknames = unique((rocknames.sed..., rocknames.met..., rocknames.ign...))

    typelist = get_rock_class(false, true)      # Subtypes, major types include minors
    nbulk = length(bulkrockname)
    bulk_lookup = NamedTuple{Symbol.(Tuple(rocknames))}([falses(nbulk) for _ in eachindex(rocknames)])
    
    p = Progress(length(rocknames)+1, desc="Finding EarthChem samples for each rock name")
    next!(p)
    for i in eachindex(rocknames)
        bulk_lookup[i] .= find_earthchem(rocknames[i], bulkrockname, bulkrocktype, bulkmaterial)
    
        # If no matches, jump up a class
        if count(bulk_lookup[i]) == 0
            newsearch = class_up(typelist, rocknames[i])
            newsearch==:carb && (newsearch=="carbonate")    # No carbonatites!
            bulk_lookup[i] .= find_earthchem(string(newsearch), bulkrockname, bulkrocktype,
                bulkmaterial
            )
    
            # If still no matches, jump up a class again
            if count(bulk_lookup[i]) == 0
                newsearch = class_up(typelist, string(newsearch))
                bulk_lookup[i] .= find_earthchem(string(newsearch), bulkrockname, bulkrocktype,
                bulkmaterial
            )
            end
        end
        next!(p)
    end


## --- Write matches to file 
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