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


## --- Get rock type matches for all samples
    # # DIY and save to a file
    # bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name])
    # bulkrocktype = lowercase.(bulktext.elements.Type[bulktext.index.Type])
    # bulkmaterial = lowercase.(bulktext.elements.Material[bulktext.index.Material])

    # bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial; unmultimatch=false, 
    #     inclusive=false, source=:earthchem
    # )

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

    # Or save yourself the time and load from a file
    fid = h5open("output/bulk_unrestricted_types.h5", "r")
    header = read(fid["bulk_cats_head"])
    data = read(fid["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    close(fid)


## --- Define and organize elements to correct
    majors, minors = get_elements()

    # Collect major and minor elements together. Make sure to keep rock type!
    # CaCO3 values will be preserved in the final bulk.h5 file, but CaCO3 is converted
    # to CaO and CO2 for analyzed wt.% calculations
    allelements = [majors; minors]
    extendelements = [allelements; [:CaCO3, :H2O, :CO2]]
    allkeys = [extendelements; [:Type, :Latitude, :Longitude, :Age]]
    ndata = length(allelements)

    # Get all units present in bulktext
    presentunits = join(collect(keys(bulktext_units.index)), " ")

    # Define densities for relevant major element oxides (g/cm³)
    density = (Al2O3=3.95, K2O=2.35, MgO=3.58, Na2O=2.27, P2O5=2.39, SiO2=2.65, TiO2=4.23, 
        CaO=3.34, FeOT=5.74
    )


## --- Convert all units to wt.%
    # This method is resistant to changes in element order
    for i in extendelements
        # Check unit is present, and convert indices to a list of units
        unit_i = string(i) * "_Unit"
        @assert contains(presentunits, unit_i) "$unit_i not present"
        bulkunits = bulktext_units.units[bulktext_units.index[Symbol(unit_i)]]

        # Convert units of all sample data
        standardize_units!(bulk[i], bulkunits, density, elem=i)
    end


## --- Carbonates: Convert CaCO₃ → CaO + CO₂, and get CO₂ from CaO
    @warn "Did you remember to re-load bulk before running this block?" 

    # Pre-compute some conversion factors
    CaCO3_to_CaO = (40.078+15.999)/(40.078+15.999*3+12.01)
    CaCO3_to_CO2 = (15.999*2+12.01)/(40.078+15.999*3+12.01)

    CaO_to_CO2 = (15.999*2+12.01)/(40.078+15.999)
    CaO_to_CaCO3 = (40.078+15.999*3+12.01)/(40.078+15.999)

    # Convert CaCO₃ to CaO and CO₂
    # Replace the existing CaO value with the computed CaO value
    # Add the existing CO2 value and the computed CO2 values 
    @turbo for i in eachindex(bulk.CaCO3)
        bulk.CaO[i] = ifelse(bulk.CaCO3[i] > 0, CaCO3_to_CaO * bulk.CaCO3[i], bulk.CaO[i])
        bulk.CO2[i] = ifelse(bulk.CaCO3[i] > 0, CaCO3_to_CO2 * bulk.CaCO3[i] + bulk.CO2[i],
            bulk.CO2[i]
        )
    end

    # FOR CARBONATES ONLY: Convert CaO to CaCO₃ and CO₂
    #
    # There are a lot of limestones with only CaO... this seems perhaps not correct, but
    # there are no CO2 values, so they all have analyzed wt.%'s that are too low :(
    #
    # It doesn't matter if existing CaCO₃ values are replaced with CaCO₃ values calculated
    # from CaO, because the CaO values were calculated from the existing CaCO₃ values so 
    # everything works (paste this code into the REPL to check):
    #
    # a = @. !isnan(bulk.CaCO3) & bulk_cats.carb;
    # b = @. !isapprox(bulk.CaCO3[a], CaCO3[a]);
    # [bulk.CaCO3[a][b] CaCO3[a][b] bulk.CaO[a][b]]
    @turbo for i in eachindex(bulk.CaO)
        bulk.CO2[i] = ifelse(bulk_cats.carb[i], bulk.CaO[i]*CaO_to_CO2, NaN)
        bulk.CaCO3[i] = ifelse(bulk_cats.carb[i], bulk.CaO[i]*CaO_to_CaCO3, NaN)
    end


## --- Gypsum: get H₂O from CaO and assume SO₄
    @warn "Did you remember to re-load bulk before running this block?" 

    # Pre-compute conversion factors
    CaO_to_H2O = 2*(2*1.00784+15.999)/(40.078+15.999)
    CaO_to_SO4 = (32.065+4*15.999)/(40.078+15.999)

    # Preallocate
    SO4 = Array{Float64}(undef, length(bulk.SiO2), 1)

    # Compute
    bulkrockname = bulktext.elements.Rock_Name[bulktext.index.Rock_Name]
    isgypsum = bulkrockname.=="GYPSUM"
    @turbo for i in eachindex(bulk.H2O)
        bulk.H2O[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_H2O, 0)
        SO4[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_SO4, 0)
    end


## --- Restrict bulk to whole rock analysis of continental crust
    # Samples must be on land; same methodology used to generate continental lat, lon 
    # points for Macrostrat API query
    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, bulk.Latitude, bulk.Longitude)
    abovesea = elev .> 0

    # Compute the larger of H₂O + CO₂ or LOI to add to the total wt.%
    bulkweight = Array{Float64}(undef, length(bulk.SiO2), 1)
    @turbo for i in eachindex(bulkweight)
        volatiles = nanadd(bulk.CO2[i], bulk.H2O[i])
        bulkweight[i] = ifelse(bulk.Loi[i] > volatiles, bulk.Loi[i], volatiles)
    end
    zeronan!(bulkweight)

    # Add calculated SO₄ to the total wt.%
    bulkweight .+= SO4

    # Samples must be 84-104 total wt.% analyzed
    t = falses(length(bulkweight))
    bounds = (84, 104)

    p = Progress(length(bulkweight), desc="Calculating wt.% ...")
    @time @inbounds for i in eachindex(bulkweight)
        if abovesea[i]
            bulkweight[i] += nansum([bulk[j][i] for j in allelements])
            t[i] = (bounds[1] <= bulkweight[i] <= bounds[2])
        else
            t[i] = false
        end
        # bulkweight[i] += nansum([bulk[j][i] for j in allelements])
        next!(p)
    end


## --- Print to terminal and normalize compositions
    nsamples = round(count(t)/length(t)*100, digits=2)
    @info "Saving $(count(t)) samples ($nsamples%) with $(bounds[1])-$(bounds[2]) analyzed wt.%"

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


## --- Restrict bulktext to data we're keeping in the final output file
    # Methods
    restrictmethods = Symbol.(string.(extendelements) .* "_Meth")
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


## --- Initialize metadata for rock name / type matching
    # Name + Type + Material from index / element format of bulktext
    bulkrockname = lowercase.(bulktext.elements.Rock_Name[bulktext.index.Rock_Name])
    bulkrocktype = lowercase.(bulktext.elements.Type[bulktext.index.Type])
    bulkmaterial = lowercase.(bulktext.elements.Material[bulktext.index.Material])

    # Rock names
    rocknames = get_rock_class(true)
    rocknames = unique((rocknames.sed..., rocknames.met..., rocknames.ign...))

    # Rock subtypes, major types include minors (for class_up function call)
    typelist = get_rock_class(false, true)
    classnames = string.(collect(keys(typelist)))


## --- Find rock type / name matches
    # Match rock types
    bulk_cats = match_rocktype(bulkrockname, bulkrocktype, bulkmaterial; unmultimatch=false, 
        inclusive=false, source=:earthchem
    )

    # Preallocate for rock names
    bulk_lookup = NamedTuple{keys(name_cats)}([falses(length(bulkrockname)) 
        for _ in eachindex(name_cats)]
    )

    p = Progress(length(rocknames), desc="Finding EarthChem samples for each rock name")
    for i in eachindex(rocknames)
        bulk_lookup[i] .= find_earthchem(rocknames[i], bulkrockname, bulkrocktype, 
            bulkmaterial
        )

        # If no matches, jump up a class. Find everything within that class
        if count(bulk_lookup[i]) == 0
            # If the name is the name of a class, DON'T jump up a class
            if rocknames[i] in classnames
                searchlist = typelist[Symbol(rocknames[i])]
            else
                searchlist = typelist[class_up(typelist, rocknames[i])]
            end

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