## --- Match Macrostrat samples to the most likely EarthChem sample
    # Packages
    using MAT
    using JLD
    using StatGeochem
    using ProgressMeter

    # Local utilities
    include("Utilities.jl")


## --- Load Earthchem bulk geochemical data
    @info "Loading EarthChem bulk data"
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
    @time geochem = major_elements(bulk, bulk_cats)

## --- Load Earthchem metadata
    @info "Loading EarthChem sample metadata"
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
    @info "Loading Macrostrat lithologic data"
    # Reduced size Macrostrat data file
    macrostrat = importdataset("data/toy_responses.tsv", '\t', importas=:Tuple)
    rocklat = macrostrat.rocklat
    rocklon = macrostrat.rocklon
    age = macrostrat.age

    # Match data to rock types
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip, major=false)


## --- For each Macrostrat sample, estimate the log likelihood that each EarthChem sample
#      accurately represents that sample in space, time, and geochemistry
# TO DO: likelihood vs log likelihood? Is this already log likelihood?
# TO DO: NaNs with 0 for age and location too?
    @info "Matching samples"

    # All EarthChem indices so we can get back to the full list when we filter data
#     bulk_idxs = collect(1:length(bulk.Type))

#     # Preallocate to store indices
#     matches = Dict{Symbol, Matrix{Int64}}()

#     # pls_less_memory(bulk_cats, bulk_idxs, geochem, rocklat, rocklon, age, macro_cats, matches)

#     # Only compare rocks of the same rock type
#     # TO DO: better progress bar?
#     @time @showprogress 0.5 "Matching lithographic / geochemical samples..." for type in eachindex(macro_cats)
        
#     # Cover is not in EarthChem data; skip it
#     if type==:cover continue end
    
#     # Get samples and data of this rock type
#     # Roughly 100k fewer allocations to have these intermediate variables
#     # chem_samples = bulk_cats[type]              # EarthChem BitVector
#     # sample_idxs = bulk_idxs[chem_samples]       # Indices of EarthChem samples for this rock type
#     # geochem_data = geochem[type]                # Major element compositions for this rock type
#     # lat = rocklat[macro_cats[type]]             # Macrostrat latitude
#     # lon = rocklon[macro_cats[type]]             # Macrostrat longitude
#     # sample_age = age[macro_cats[type]]          # Macrostrat age

#     matched_sample = likelihood(bulk_cats, bulk_idxs, macro_cats, type)
    
#     # Preallocate array of index of most likely EarthChem sample for each Macrostrat sample
#     # matched_sample = Array{Int64}(undef, length(lat), 1)

#     # # For each Macrostrat sample, calculate
#     # for j in eachindex(lat)
#     #     # Distance (σ = 1.8 arc degrees)
#     #     dist = arcdistance(lat[j], lon[j], bulk.Latitude[chem_samples], bulk.Longitude[chem_samples])
#     #     lh_dist = -(dist.^2)./(1.8^2)

#     #     # Age (σ = 38 Ma)
#     #     # TO DO: AgeEst vs Age? Maybe recalculate AgeEst?
#     #     age_diff = abs.(bulk.AgeEst[chem_samples] .- sample_age[j])
#     #     lh_age = -(age_diff.^2)./(38^2)

#     #     # Geochemistry
#     #     # TO DO: only look at ones that are close in age and space
#     #     lh_geochem = zeros(Float64, count(chem_samples))
#     #     for elem in eachindex(geochem_data)
#     #         # Get all EarthChem samples for that element and this rock type
#     #         # Replace NaNs with 0
#     #         bulk_geochem = zeronan!(bulk[elem][chem_samples])

#     #         # Assume Macrostrat sample is the average geochemistry for that rock type
#     #         # TO DO: More granular estimate than just by rock type
#     #         geochem_diff = abs.(bulk_geochem .- geochem_data[elem].m)
            
#     #         # Calculate likelihood for that element and add to the other likelihoods
#     #         lh_elem = -(geochem_diff.^2)./(geochem_data[elem].e^2)
#     #         lh_geochem .+= lh_elem
#     #     end

#     #     # lh_geochem = geochem_likelihood(bulk_cats, type, geochem)

#     #     # Calculate total likelihood for each EarthChem sample
#     #     # This has to be added in steps because nanadd can only do one array at a time
#     #     lh_total = nanadd(lh_dist, lh_age)
#     #     lh_total = nanadd(lh_total, lh_geochem)

#     #     # Get the index of the most likely EarthChem sample
#     #     # TO DO: the most likely sample has the largest likelihood? because if the 
#     #         # difference is larger than the total likelihood value is more negative...
#     #     # TO DO: take the average of some arbitrary percentile
#     #     (val, idx) = findmax(abs.(lh_total[findall(!isnan, lh_total)]))
#     #     matched_sample[j] = sample_idxs[findall(!isnan, lh_total)][idx]  # Replaced sample_idxs here, check for bugs
#     # end

#     # Add the list of matches 
#     # NOTE: no chert matches because Macrostrat data does not have chert
#     setindex!(matches, matched_sample, type)
# end

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

## -- reduced loop
bulk_idxs = collect(1:length(bulk.Type))
matches = Dict{Symbol, Matrix{Int64}}()
@time @showprogress "Matching lithographic / geochemical samples..." for type in eachindex(macro_cats)
    # Cover is not in EarthChem data; skip it
    if type==:cover continue end

    # Get intermediate variables
    chem_samples = bulk_cats[type]              # EarthChem BitVector
    sample_idxs = bulk_idxs[chem_samples]       # Indices of EarthChem samples for this rock type
    geochem_data = geochem[type]                # Major element compositions for this rock type
    lat = rocklat[macro_cats[type]]             # Macrostrat latitude
    lon = rocklon[macro_cats[type]]             # Macrostrat longitude
    bulklat = bulk.Latitude[chem_samples]       # EarthChem latitudes
    bulklon = bulk.Longitude[chem_samples]      # EarthChem longitudes
    sample_age = age[macro_cats[type]]          # Macrostrat age

    # Preallocate
    # matched_sample = Array{Int64}(undef, length(lat), 1)

    # Find most likely sample
    matched_sample = likelihood(chem_samples, sample_idxs, geochem_data, lat, lon, sample_age, bulklat, bulklon)
    
    setindex!(matches, matched_sample, type)
end


