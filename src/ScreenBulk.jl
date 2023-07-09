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
        units = bulktext["Methods"],
        index = NamedTuple{Tuple(Symbol.(collect(keys(bulktext["method"]))))}(
            [bulktext["method"][i] for i in keys(bulktext["method"])]
        )
    )

    # Information by sample
    bulktext = (
        elements = NamedTuple{Tuple(Symbol.(collect(keys(bulktext["elements"]))))}(
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


## --- Restrict data to 84-104 wt.% analyzed
    bulkweight = zeros(length(bulk.SiO2));
    @time for i in eachindex(bulkweight)
        bulkweight[i] = nansum([bulk[j][i] for j in allelements])
    end

    bounds = (84, 104)
    t = @. bounds[1] < bulkweight < bounds[2]

    nsamples = round(count(t)/length(t)*100, digits=2)
    @info "Saving $nsamples% samples between $(bounds[1])% and $(bounds[2])% analyzed wt.%"

    bulknew = NamedTuple{Tuple(allkeys)}([fill(NaN, count(t)) for i in allkeys])
    for i in allkeys
        bulknew[i] .= bulk[i][t]
    end

    # Normalize compositions to 100%
    for i in eachindex(bulknew.SiO2)
        sample = [bulknew[j][i] for j in allelements]   # Get it
        normalize!(sample)                              # Normalize it
        for j in eachindex(allelements)                 # Put it back
            bulknew[allelements[j]][i] = sample[j]
        end
    end

    
## --- Restrict 


## --- Write to an HDF5 file
    fid = h5open("data/bulk.h5", "w")

    # Bulk data
    data = create_group(fid, "bulk")

    write(data, "header", string.(collect(keys(bulk))))


    # Bulktext metadata
    text = create_group(fid, "bulktext")
    bulktext_methods_g = 
    bulktext_units_g = 

    close(fid)

## --- Normalize remaining compositions to 100%
    

    # Check that elements were normalized correctly
    # for i in eachindex(t)
    #     !t[i] && continue

    #     samplesum = 0.0
    #     for j in eachindex(allelements)
    #         samplesum = nanadd(samplesum, bulknew[string(allelements[j])][i])
    #     end
    #     @assert isapprox(samplesum, 100)
    # end

    # Save to file
    matwrite("data/bulknew_norm.mat", bulknew)


## --- End of File