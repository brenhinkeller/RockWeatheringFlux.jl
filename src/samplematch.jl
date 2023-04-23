## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem

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
    chem_composition = lowercase.(string.(bulktext.Composition[composition_idx]))
    chem_reference = lowercase.(string.(bulktext.Reference[reference_idx]))
    chem_rockname = lowercase.(string.(bulktext.Rock_Name[rockname_idx]))
    chem_source = lowercase.(string.(bulktext.Source[source_idx]))
    chem_type = lowercase.(string.(bulktext.Type[type_idx]))
    chem_material = lowercase.(string.(bulktext.Material[material_idx]))


    
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


## --- For each sample, estimate the log likelihood that... what?
    # Only compare rocks of the same rock type
    for elem in eachindex(macro_cats)
        # Cover is not present in earthchem data; skip it
        if elem==:cover
            continue
        end
        
        # Get EarthChem samples of this rock type
        chem_samples = bulk_cats[elem]     # EarthChem

        # Mean and standard devation of major element oxides for this rock type
        major_oxides = (
            SiO2 = (m = nanmean(bulk.SiO2[chem_samples]), e = nanstd(bulk.SiO2[chem_samples])),
            Al2O3 = (m = nanmean(bulk.Al2O3[chem_samples]), e = nanstd(bulk.Al2O3[chem_samples])),
            Fe2O3T = (m = nanmean(bulk.Fe2O3T[chem_samples]), e = nanstd(bulk.Fe2O3T[chem_samples])),
            TiO2 = (m = nanmean(bulk.TiO2[chem_samples]), e = nanstd(bulk.TiO2[chem_samples])),
            MgO = (m = nanmean(bulk.MgO[chem_samples]), e = nanstd(bulk.MgO[chem_samples])),
            CaO = (m = nanmean(bulk.CaO[chem_samples]), e = nanstd(bulk.CaO[chem_samples])),
            Na2O = (m = nanmean(bulk.Na2O[chem_samples]), e = nanstd(bulk.Na2O[chem_samples])),
            K2O = (m = nanmean(bulk.K2O[chem_samples]), e = nanstd(bulk.K2O[chem_samples])),
            MnO = (m = nanmean(bulk.MnO[chem_samples]), e = nanstd(bulk.MnO[chem_samples]))
        )

        # Get Macrostrat samples of that rock type
        lat = rocklat[macro_cats[elem]]     # Latitude
        lon = rocklon[macro_cats[elem]]     # Longitude
        smpl_age = age[macro_cats[elem]]    # Age

        for i in eachindex(lat)
            # Distance from point of interest (σ = 1.8 arc degrees)
            dist = arcdistance(lat[i], lon[i], bulk.Latitude, bulk.Longitude)
            lh_dist = -(dist.^2)./(1.8^2)

            # Age from point of interest (σ = 38 Ma)
            # Age vs AgeEst? Latter gives more data...
            age_diff = abs.(bulk.AgeEst .- smpl_age[i])
            lh_age = -(age_diff.^2)./(38^2)

            # Geochemical difference from point of interest
            # Assign a likely 
        end
        
        
    end

### --- End of file