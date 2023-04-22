## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem

    # Local utilities
    include("Utilities.jl")         # Depending on REPL or script, only one will fail
    include("src/Utilities.jl")     # Failure doesn't do anything so just run both every time


## --- Load Earthchem bulk geochemical data
    earthchem_raw = matopen("data/bulk.mat")
    earthchem_dict = read(earthchem_raw, "bulk")
    close(earthchem_raw)
    earthchem = NamedTuple{Tuple(Symbol.(keys(earthchem_dict)))}(values(earthchem_dict))

    # Filter ages younger than 0 or greater than the age of the earth
    invalid_age = vcat(findall(>(4000), earthchem.Age), findall(<(0), earthchem.Age))
    earthchem.Age[invalid_age] .= NaN

    # Get rock types
    echem_cats = match_earthchem(earthchem.Type, major=false)

## --- Load Earthchem metadata
    earthchem_raw = matopen("data/bulktext.mat")
    earthchem_dict = read(earthchem_raw, "bulktext")
    close(earthchem_raw)
    chemtext_unparsed = NamedTuple{Tuple(Symbol.(keys(earthchem_dict)))}(values(earthchem_dict))
    
    # Preallocate for parsing---make a NamedTuple with just the things we want
    npoints = length(chemtext_unparsed.index["Composition"])

    # Ignoring sparse arrays for now because they make me sad
    # Correct indices for 0-index offset
    composition_idx = Int.(chemtext_unparsed.index["Composition"] .+ 1)
    reference_idx = Int.(chemtext_unparsed.index["Reference"] .+ 1)
    rockname_idx = Int.(chemtext_unparsed.index["Rock_Name"] .+ 1)
    source_idx = Int.(chemtext_unparsed.index["Source"] .+ 1)
    type_idx = Int.(chemtext_unparsed.index["Type"] .+ 1)
    material_idx = Int.(chemtext_unparsed.index["Material"] .+ 1)

    # Parse numeric codes in index into arrays
    chem_composition = lowercase.(string.(chemtext_unparsed.Composition[composition_idx]))
    chem_reference = lowercase.(string.(chemtext_unparsed.Reference[reference_idx]))
    chem_rockname = lowercase.(string.(chemtext_unparsed.Rock_Name[rockname_idx]))
    chem_source = lowercase.(string.(chemtext_unparsed.Source[source_idx]))
    chem_type = lowercase.(string.(chemtext_unparsed.Type[type_idx]))
    chem_material = lowercase.(string.(chemtext_unparsed.Material[material_idx]))

    
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
    # All rocks must be the same rock type
    for elem in eachindex(macro_cats)
        # Cover is not present in earthchem data; skip it
        if elem==:cover
            continue
        end
        
        # Get EarthChem samples of that rock type
        chem_samples = echem_cats[elem]      # EarthChem

        # Get Macrostrat samples of that rock type
        lat = rocklat[macro_cats[elem]]     # Latitude
        lon = rocklon[macro_cats[elem]]     # Longitude
        smpl_age = age[macro_cats[elem]]    # Age

        for i in eachindex(lat)
            # Distance from point of interest (σ = 1.8 arc degrees)
            dist = arcdistance(lat[i], lon[i], earthchem.Latitude, earthchem.Longitude)
            lh_dist = -(dist.^2)./(1.8^2)

            # Age from point of interest (σ = 38 Ma)
            # Age vs AgeEst? Latter gives more data...
            age_diff = abs.(earthchem.AgeEst .- smpl_age[i])
            lh_age = -(age_diff.^2)./(38^2)

            # Geochemical difference from point of interest
        end
        
        
    end

### --- End of file