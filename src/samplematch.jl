## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem

    # Local utilities
    include("Utilities.jl")


## --- Load Earthchem bulk geochemical data
    earthchem_raw = matopen("data/bulk.mat")
    earthchem_dict = read(earthchem_raw, "bulk")
    close(earthchem_raw)
    earthchem = NamedTuple{Tuple(Symbol.(keys(earthchem_dict)))}(values(earthchem_dict))

    cats_echem = match_earthchem(earthchem.Type, major=false)   # Get rock types

## --- Load Earthchem metadata
    earthchem_raw = matopen("data/bulktext.mat")
    earthchem_dict = read(earthchem_raw, "bulktext")
    close(earthchem_raw)
    text_echem_unparsed = NamedTuple{Tuple(Symbol.(keys(earthchem_dict)))}(values(earthchem_dict))
    
    # Preallocate for parsing---make a NamedTuple with just the things we want
    npoints = length(text_echem_unparsed.index["Composition"])

    # Ignoring sparse arrays for now because they make me sad
    # Correct indices for 0-index offset
    composition_idx = Int.(text_echem_unparsed.index["Composition"] .+ 1)
    reference_idx = Int.(text_echem_unparsed.index["Reference"] .+ 1)
    rockname_idx = Int.(text_echem_unparsed.index["Rock_Name"] .+ 1)
    source_idx = Int.(text_echem_unparsed.index["Source"] .+ 1)
    type_idx = Int.(text_echem_unparsed.index["Type"] .+ 1)
    material_idx = Int.(text_echem_unparsed.index["Material"] .+ 1)

    # Parse numeric codes in index into arrays
    echem_composition = lowercase.(string.(text_echem_unparsed.Composition[composition_idx]))
    echem_reference = lowercase.(string.(text_echem_unparsed.Reference[reference_idx]))
    echem_rockname = lowercase.(string.(text_echem_unparsed.Rock_Name[rockname_idx]))
    echem_source = lowercase.(string.(text_echem_unparsed.Source[source_idx]))
    echem_type = lowercase.(string.(text_echem_unparsed.Type[type_idx]))
    echem_material = lowercase.(string.(text_echem_unparsed.Material[material_idx]))

    
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

### --- End of file