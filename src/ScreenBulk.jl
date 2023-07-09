# This file will read bulk.mat and bulktext.mat and standardize all units as wt.%
# Data will be saved as a new file
    # Currently a new .mat file

# TO DO: 
    # save as an .h5 file. Requires re-writing other files so putting off for now.
    # Re-read in old files, and figure out what to do for parsing / correcting element 
        # oxide weights. For example, some samples may have different values for CaO and 
        # CaCO3, but Mg and MgO represent the same data.
    # Consider making multiple utilities files...
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

    # Local utilities
    include("Utilities.jl")
    
    # Load and parse data files
    bulk = matread("data/bulk.mat")["bulk"]
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))


## --- Load and parse metadata
    bulktext = matread("data/bulktext.mat")["bulktext"]

    # Type stabilize sample metadata
    firstshell = unique([String.(bulktext["elements"]); ["Units", "Methods"]])
    for n in firstshell
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
    allelements = [majors; minors]
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


## --- Restrict bulk to 84-104 wt.% analyzed
    bulkweight = zeros(length(bulk.SiO2));
    @time for i in eachindex(bulkweight)
        bulkweight[i] = nansum([bulk[j][i] for j in allelements])
    end

    bounds = (84, 104)
    t = @. bounds[1] < bulkweight < bounds[2]

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


## --- Write to an HDF5 file
    fid = h5open("data/bulk.h5", "w")
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

    close(fid)


## --- End of File