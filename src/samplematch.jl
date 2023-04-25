## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem
    using ProgressMeter

    # Local utilities
    include("Utilities.jl")         # Depending on REPL or script, only one will fail


## --- Load Earthchem bulk geochemical data
    bulk_raw = matopen("data/bulk.mat")
    bulk_dict = read(bulk_raw, "bulk")
    close(bulk_raw)
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

    # Filter ages younger than 0 or greater than the age of the earth
    invalid_age = vcat(findall(>(4000), bulk.Age), findall(<(0), bulk.Age))
    bulk.Age[invalid_age] .= NaN

    # Get rock types
    bulk_cats = match_earthchem(bulk.Type, major=false)

## --- Load Earthchem metadata
    bulktext_raw = matopen("data/bulktext.mat")
    bulktext_dict = read(bulktext_raw, "bulktext")
    close(bulktext_raw)
    bulktext = NamedTuple{Tuple(Symbol.(keys(bulktext_dict)))}(values(bulktext_dict))
    
    # Preallocate for parsing---make a NamedTuple with just the things we want
    npoints = length(bulktext.index["Composition"])

    # Ignoring sparse arrays for now because they make me sad
    # Correct indices for 0-index offset
    composition_idx = Int.(bulktext.index["Composition"] .+ 1)
    reference_idx = Int.(bulktext.index["Reference"] .+ 1)
    rockname_idx = Int.(bulktext.index["Rock_Name"] .+ 1)
    source_idx = Int.(bulktext.index["Source"] .+ 1)
    type_idx = Int.(bulktext.index["Type"] .+ 1)
    material_idx = Int.(bulktext.index["Material"] .+ 1)

    # Parse numeric codes in index into arrays
    bulk_composition = lowercase.(string.(bulktext.Composition[composition_idx]))
    bulk_reference = lowercase.(string.(bulktext.Reference[reference_idx]))
    bulk_rockname = lowercase.(string.(bulktext.Rock_Name[rockname_idx]))
    bulk_source = lowercase.(string.(bulktext.Source[source_idx]))
    bulk_type = lowercase.(string.(bulktext.Type[type_idx]))
    bulk_material = lowercase.(string.(bulktext.Material[material_idx]))

    # All relevant rock type identifiers
    bulk_lith = bulk_type .* " " .* bulk_material .* " " .* bulk_rockname

## --- To do:
    #=
    Quality control on the rock names--some very mafic granites in this dataset
    Or honestly maybe just skip the rock names and just go on silica content? Idk
    =#


## --- Calculate mean and standard deviation of major element oxides for each rock type
    # In the future, we will want to separate igneous rocks by silica content--this is fine for now
    
    # Define major elements from Faye and Ødegård 1975

    # Temporary data storage
    geochem = Array{NamedTuple}(undef, length(bulk_cats), 1)
    for i = 1:length(bulk_cats)
        type = bulk_cats[i]
        elem = (
            SiO2 = (m = nanmean(bulk.SiO2[type]), e = nanstd(bulk.SiO2[type])),
            Al2O3 = (m = nanmean(bulk.Al2O3[type]), e = nanstd(bulk.Al2O3[type])),
            Fe2O3T = (m = nanmean(bulk.Fe2O3T[type]), e = nanstd(bulk.Fe2O3T[type])),
            TiO2 = (m = nanmean(bulk.TiO2[type]), e = nanstd(bulk.TiO2[type])),
            MgO = (m = nanmean(bulk.MgO[type]), e = nanstd(bulk.MgO[type])),
            CaO = (m = nanmean(bulk.CaO[type]), e = nanstd(bulk.CaO[type])),
            Na2O = (m = nanmean(bulk.Na2O[type]), e = nanstd(bulk.Na2O[type])),
            K2O = (m = nanmean(bulk.K2O[type]), e = nanstd(bulk.K2O[type])),
            MnO = (m = nanmean(bulk.MnO[type]), e = nanstd(bulk.MnO[type]))
        )

        geochem[i] = elem
    end

    # This is objectively a horrible way to get this data into the tuple
    # This is also definately one option I have to store data!
    # To do: load from a file?
    avg_geochem = (
        alluvium = geochem[1],
        siliciclast = geochem[2],
        shale = geochem[3],
        carb = geochem[4],
        chert = geochem[5],
        evaporite = geochem[6],
        phosphorite = geochem[7],
        coal = geochem[8],
        volcaniclast = geochem[9],
        sed = geochem[10],

        volc = geochem[11],
        plut = geochem[12],
        ign = geochem[13],

        metased = geochem[14],
        metaign = geochem[15],
        met = geochem[16],
    )


## --- Load Macrostrat data
    # Reduced size Macrostrat data file
    retrive_file = load("data/toy_responses.jld")
    responses = retrive_file["responses"]
    elevations = float.(retrive_file["elevations"])
    rocklat = float.(retrive_file["latitude"])
    rocklon = float.(retrive_file["longitude"])    
    npoints = retrive_file["npoints"]

    # Preallocate for parsed Macrostrat data
    rocktype = Array{String}(undef, npoints, 1)
    rockdescrip = Array{String}(undef, npoints, 1)
    rockname = Array{String}(undef, npoints, 1)
    rockstratname = Array{String}(undef, npoints, 1)
    rockcomments = Array{String}(undef, npoints, 1)
    agemax = Array{Float64}(undef, npoints, 1)
    agemin = Array{Float64}(undef, npoints, 1)
    age = Array{Float64}(undef, npoints, 1)

    # Parse responses into preallocated arrays
    for i in eachindex(rocktype)
        rocktype[i] = get_macrostrat_lith(responses[i])
        rockdescrip[i] = get_macrostrat_descrip(responses[i])
        rockname[i] = get_macrostrat_name(responses[i])
        rockstratname[i] = get_macrostrat_strat_name(responses[i])
        rockcomments[i] = get_macrostrat_comments(responses[i])
        agemax[i] = get_macrostrat_max_age(responses[i])
        agemin[i] = get_macrostrat_min_age(responses[i])
    end

    # Age of each sample
    for i in 1:npoints
        age[i] = nanmean([agemax[i], agemin[i]])
    end

    # Filter ages younger than 0 or greater than the age of the earth
    invalid_age = vcat(findall(>(4000), age), findall(<(0), age))
    age[invalid_age] .= NaN

    # Make sure age bounds are in the right order
    for i in 1:npoints
        if agemin[i] > agemax[i]
            tempmax = agemin[i]
            tempmin = agemax[i]
            agemin[i] = tempmax
            agemax[i] = tempmin
        end
    end

    # Convert strings to lowercase so they can be matched to known names of rock types
    rocktype = lowercase.(rocktype)
    rockdescrip = lowercase.(rockdescrip)
    rockname = lowercase.(rockname)
    rockstratname = lowercase.(rockstratname)
    rockcomments = lowercase.(rockcomments)

    # Replace tabs with spaces so they will not be confused with the delimitator if exported
    rocktype = replace.(rocktype, "    " => " ")
    rockdescrip = replace.(rockdescrip, "    " => " ")
    rockname = replace.(rockname, "    " => " ")
    rockstratname = replace.(rockstratname, "    " => " ")
    rockcomments = replace.(rockcomments, "    " => " ")

    # Match data to rock types
    macro_cats = match_rocktype(rocktype, rockname, rockdescrip, major=false)


## --- For each Macrostrat sample, estimate the log likelihood that each EarthChem sample
#      accurately represents that sample in space, time, and geochemistry
# TO DO: likelihood vs log likelihood? Is this already log likelihood?

    # Only compare rocks of the same rock type
    @showprogress 1 "Estimating likelihoods for each sample" for type in eachindex(macro_cats)
        # Cover is not present in earthchem data; skip it
        if type==:cover
            continue
        end
        
        # Get samples and data of this rock type
        # TO DO: is it more memory efficient to not assign these to be variables and just index directly?
        chem_samples = bulk_cats[type]          # EarthChem BitVector
        geochem_data = avg_geochem[type]        # Major element averages
        lat = rocklat[macro_cats[type]]         # Macrostrat latitude
        lon = rocklon[macro_cats[type]]         # Macrostrat longitude
        sample_age = age[macro_cats[type]]      # Macrostrat age

        # For each Macrostrat sample, calculate
        for i in eachindex(lat)
            # Distance (σ = 1.8 arc degrees)
            dist = arcdistance(lat[i], lon[i], bulk.Latitude[chem_samples], bulk.Longitude[chem_samples])
            lh_dist = -(dist.^2)./(1.8^2)

            # Age (σ = 38 Ma)
            # TO DO: AgeEst vs Age? Maybe recalculate AgeEst?
            age_diff = abs.(bulk.AgeEst[chem_samples] .- sample_age[i])
            lh_age = -(age_diff.^2)./(38^2)

            # Geochemistry
            # TO DO: only look at ones that are close in age and space
            # TO DO: replace NaNs with 0s
            lh_geochem = zeros(Float64, count(chem_samples))
            for elem in eachindex(geochem_data)
                # Get all EarthChem samples for that element and this rock type
                # TO DO: memory efficiency of assigning this to a variable vs. straight indexing in?
                bulk_geochem = bulk[elem][chem_samples]

                # Assume Macrostrat sample is the average geochemistry for that rock type
                # TO DO: More granular estimate than just by rock type
                geochem_diff = abs.(bulk_geochem .- geochem_data[elem].m)
                
                # Calculate likelihood for that element and add to the other likelihoods
                lh_elem = -(geochem_diff.^2)./(geochem_data[elem].e^2)
                lh_geochem .+= lh_elem
            end

            # Calculate total likelihood for each EarthChem sample
            # NOTE: this won't work for most samples until I fix the NaN / 0 issue
            lh_total = lh_dist .+ lh_age .+ lh_geochem

            # Get the index of the most likely EarthChem sample
            # TO DO: the most likely sample has the smallest likelihood? because if the 
                # difference is larger than the total likelihood value is larger...
            # TO DO: take the average of some arbitrary percentile
            (val, idx) = findmin(abs.(lh_total[findall(!isnan, lh_total)]))
        end
    end

### --- End of file