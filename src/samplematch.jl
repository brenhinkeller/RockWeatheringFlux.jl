## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem
    using ProgressMeter

    # Local utilities
    include("src/Utilities.jl")


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

    # Calculate mean and standard deviation of major element oxides for each rock type
    # TO DO: quality control on bulk rock names: why are there granites with 50% silica?
    # TO DO: separate igneous rocks by silica content
    geochem = major_elements(bulk, bulk_cats)

## --- Load Earthchem metadata
    bulktext_raw = matopen("data/bulktext.mat")
    bulktext_dict = read(bulktext_raw, "bulktext")
    close(bulktext_raw)
    bulktext = NamedTuple{Tuple(Symbol.(keys(bulktext_dict)))}(values(bulktext_dict))

    # Ignoring sparse arrays for now because they make me sad
    # Correct indices for 0-index offset
    composition_idx = Int.(bulktext.index["Composition"] .+ 1)
    reference_idx = Int.(bulktext.index["Reference"] .+ 1)
    rockname_idx = Int.(bulktext.index["Rock_Name"] .+ 1)
    source_idx = Int.(bulktext.index["Source"] .+ 1)
    type_idx = Int.(bulktext.index["Type"] .+ 1)
    material_idx = Int.(bulktext.index["Material"] .+ 1)

    # Parse numeric codes in index into arrays
    # TO DO: NamedTuple?
    bulk_composition = lowercase.(string.(bulktext.Composition[composition_idx]))
    bulk_reference = lowercase.(string.(bulktext.Reference[reference_idx]))
    bulk_rockname = lowercase.(string.(bulktext.Rock_Name[rockname_idx]))
    bulk_source = lowercase.(string.(bulktext.Source[source_idx]))
    bulk_type = lowercase.(string.(bulktext.Type[type_idx]))
    bulk_material = lowercase.(string.(bulktext.Material[material_idx]))

    # All relevant rock type identifiers
    bulk_lith = bulk_type .* " " .* bulk_material .* " " .* bulk_rockname


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
# TO DO: NaNs with 0 for age and location too?

    # All EarthChem indices so we can get back to the full list when we filter data
    bulk_idxs = collect(1:length(bulk.Type))

    # Preallocate to store indices
    matches = Dict{Symbol, Matrix{Int64}}()

    # Only compare rocks of the same rock type
    # TO DO: better progress bar?
    @showprogress 0.5 "Matching lithographic / geochemical samples..." for type in eachindex(macro_cats)
        
        # Cover is not in EarthChem data; skip it
        if type==:cover continue end
        
        # Get samples and data of this rock type
        # TO DO: is it more memory efficient to not assign these to be variables and just index directly?
        chem_samples = bulk_cats[type]              # EarthChem BitVector
        sample_idxs = bulk_idxs[bulk_cats[type]]    # Indices of EarthChem samples
        geochem_data = geochem[type]            # Major element averages
        lat = rocklat[macro_cats[type]]             # Macrostrat latitude
        lon = rocklon[macro_cats[type]]             # Macrostrat longitude
        sample_age = age[macro_cats[type]]          # Macrostrat age
        
        # Preallocate array of index of most likely EarthChem sample for each Macrostrat sample
        matched_sample = Array{Int64}(undef, length(lat), 1)

        # For each Macrostrat sample, calculate
        for j in eachindex(lat)
            # Distance (σ = 1.8 arc degrees)
            dist = arcdistance(lat[j], lon[j], bulk.Latitude[chem_samples], bulk.Longitude[chem_samples])
            lh_dist = -(dist.^2)./(1.8^2)

            # Age (σ = 38 Ma)
            # TO DO: AgeEst vs Age? Maybe recalculate AgeEst?
            age_diff = abs.(bulk.AgeEst[chem_samples] .- sample_age[j])
            lh_age = -(age_diff.^2)./(38^2)

            # Geochemistry
            # TO DO: only look at ones that are close in age and space
            lh_geochem = zeros(Float64, count(chem_samples))
            for elem in eachindex(geochem_data)
                # Get all EarthChem samples for that element and this rock type
                # Replace NaNs with 0
                # TO DO: memory efficiency of assigning this to a variable vs. straight indexing in?
                bulk_geochem = zeronan!(bulk[elem][chem_samples])

                # Assume Macrostrat sample is the average geochemistry for that rock type
                # TO DO: More granular estimate than just by rock type
                geochem_diff = abs.(bulk_geochem .- geochem_data[elem].m)
                
                # Calculate likelihood for that element and add to the other likelihoods
                lh_elem = -(geochem_diff.^2)./(geochem_data[elem].e^2)
                lh_geochem .+= lh_elem
            end

            # Calculate total likelihood for each EarthChem sample
            # This has to be added in steps because nanadd can only do one array at a time
                # TO DO: PR that?
            lh_total = nanadd(lh_dist, lh_age)
            lh_total = nanadd(lh_total, lh_geochem)

            # Get the index of the most likely EarthChem sample
            # TO DO: the most likely sample has the largest likelihood? because if the 
                # difference is larger than the total likelihood value is more negative...
            # TO DO: take the average of some arbitrary percentile
            (val, idx) = findmax(abs.(lh_total[findall(!isnan, lh_total)]))
            matched_sample[j] = sample_idxs[findall(!isnan, lh_total)][idx]
        end

        # Add the list of matches 
        # NOTE: no chert matches because Macrostrat data does not have chert
        setindex!(matches, matched_sample, type)
    end

### --- End of file


#=
Test a subset of potential matches

s = matches[:ign][1]
i = findfirst(macro_cats.ign)

println("
    $(rocktype[i]) 
    $(rockdescrip[i]) 
    $(rockname[i]) 
    $(rockstratname[i]) 
    $(rockcomments[i]) 
    $(age[i])"
)
bulk_lith[s]
bulk.AgeEst[s]
=#