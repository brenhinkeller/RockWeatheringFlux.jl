## --- Setup
    # External packages
    using ProgressMeter: @showprogress
    using StatGeochem
    using DelimitedFiles
    using Plots
    using Dates
    # using LoopVectorization

    # File parsing packages
    using JLD
    using HDF5
    using HTTP
    using JSON
    using MAT

    # Local utilities
    include("Utilities.jl")


## --- Get some random points on the continental crust
    # # Generate random points
    # npoints = 100000
    # etopo = get_etopo("elevation")
    # rocklat, rocklon, elevations = gen_continental_points(npoints, etopo)

    # # Initialize
    # zoom = 11
    # savefilename = "responses2"
    # responses = Array{Any}(undef, npoints, 1)

    # # Request data from API
    # @showprogress 5 for i = 1:npoints
    # try
    #     responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
    # catch
    #     @warn "No response from Macrostrat server for coordinate $i/$npoints. Trying again in 5 seconds. \n"
    #     try
    #     # Wait and try again
    #     sleep(5)
    #     responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
    #     catch
    #     # If still nothing, add warning
    #     responses[i] = "No response"
    #     @warn "No data from Macrostrat server for coordinate $i/$npoints\n"
    #     end
    # end
    # sleep(0.05)

    # # Checkpoint save every 10,000 points
    # if mod(i,10000)==0
    #     save("data/$savefilename.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)
    # end
    # end

    # # Save the file
    # save("data/$savefilename.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)
    

## --- Alternatively, load pre-generated data
    @info "Loading pregenerated data. This may take up to 30 minutes on slow machines. Started $(Dates.format(now(), "HH:MM"))"
    retrive_file = load("data/pregenerated_responses.jld")

    responses = retrive_file["responses"]
    elevations = float.(retrive_file["elevations"])

    rocklat = float.(retrive_file["latitude"])
    rocklon = float.(retrive_file["longitude"])    
    npoints = retrive_file["npoints"]


## --- Parse Macrostrat responses from the loaded data
    # Preallocate
    rocktype = Array{String}(undef, npoints, 1)
    rockdescrip = Array{String}(undef, npoints, 1)
    rockname = Array{String}(undef, npoints, 1)
    rockstratname = Array{String}(undef, npoints, 1)
    rockcomments = Array{String}(undef, npoints, 1)
    agemax = Array{Float64}(undef, npoints, 1)
    agemin = Array{Float64}(undef, npoints, 1)

    # Parse responses into preallocated arrays
    for i = eachindex(rocktype)
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


## --- Parse Macrostat data references
    # Preallocate
    authors = Array{String}(undef, npoints, 1)
    years = Array{String}(undef, npoints, 1)
    titles = Array{String}(undef, npoints, 1)
    dois = Array{String}(undef, npoints, 1)

    # Parse into preallocated arrays
    for i = eachindex(authors)
    authors[i] = get_macrostrat_ref_authors(responses[i])
    years[i] =get_macrostrat_ref_authors(responses[i])
    titles[i] = get_macrostrat_ref_authors(responses[i])
    dois[i] = get_macrostrat_ref_authors(responses[i])
    end

    # Write the references to a file
    refstrings = @. authors * " | " * years * " | " * titles * " | " * dois
    writedlm("output/burwellrefs.tsv",unique(refstrings),"\t")


## --- Match the Burwell rocktype with our rock types
    sed, ign, met, cover = match_rocktype(rocktype, rockname, rockdescrip, major=true)

    # Exclude suspected cover from the other three categories, just to be sure
    sed .&= .!cover .& .!ign
    ign .&= .!cover
    met .&= .!cover .& .!sed .& .!ign   # Metaseds / metaigns grouped with sed / ign
    cryst = ign .| met                  # Define crystalline rocks as igneous or non metased/ign metamorphic

    # Figure out how many data points weren't matched
    known = sed .| ign .| met
    matched = known .| cover
    not_matched = .!matched
    multi_matched = (sed .& ign) .| (sed .& met) .| (ign .& met)
    total_known = count(known)

    # Relative abundance
    sed_area_frac = count(sed) / total_known
    ign_area_frac = count(ign) / total_known
    met_area_frac = count(met) / total_known
    cryst_area_frac = count(cryst) / total_known

    # Print to terminal
    @info "not matched = $(count(not_matched)), conflicting matches = $(count(multi_matched))\n"

    @info "Rock type totals and relative abundance (multi matched are counted twice):
    sed = $(sum(sed)) ($(round(sed_area_frac*100, digits=2))%)
    ign = $(sum(ign)) ($(round(ign_area_frac*100, digits=2))%)
    met = $(sum(met)) ($(round(met_area_frac*100, digits=2))%)
    cryst (ign + (met & ~sed)) = $(sum(cryst)) ($(round(cryst_area_frac*100, digits=2))%)\n"


## -- Write the data to .tsv files so we can access or export it
    writedlm("output/notmatched.tsv", hcat(rocklat[not_matched], rocklon[not_matched], rocktype[not_matched], rockname[not_matched], 
        rockdescrip[not_matched], rockstratname[not_matched], rockcomments[not_matched]))
    writedlm("output/multimatched.tsv", hcat(rocklat[multi_matched], rocklon[multi_matched], rocktype[multi_matched], 
        rockname[multi_matched], rockdescrip[multi_matched], rockstratname[multi_matched], rockcomments[multi_matched]))
    writedlm("output/sed.tsv", hcat(rocklat[sed], rocklon[sed], rocktype[sed], rockname[sed], rockdescrip[sed], rockstratname[sed], 
        rockcomments[sed]))
    writedlm("output/ign.tsv", hcat(rocklat[ign], rocklon[ign], rocktype[ign], rockname[ign], rockdescrip[ign], rockstratname[ign], 
        rockcomments[ign]))
    writedlm("output/met.tsv", hcat(rocklat[met], rocklon[met], rocktype[met], rockname[met], rockdescrip[met], rockstratname[met], 
        rockcomments[met]))

    # writedlm("output/plutonic.tsv", hcat(rocktype[plut], rockname[plut], rockdescrip[plut], 
    #   rockstratname[plut], rockcomments[plut]))
    # writedlm("output/intrusive.tsv", hcat(rocktype[hypabyssal .& .~plut], rockname[hypabyssal .& .~plut], rockdescrip[hypabyssal .& .~plut], 
    #   rockstratname[hypabyssal .& .~plut], rockcomments[hypabyssal .& .~plut]))

    t = elevations .> 4000
    writedlm("output/highelev.tsv", hcat(rocklat[t], rocklon[t], rocktype[t], rockname[t], rockdescrip[t], rockstratname[t], 
        rockcomments[t]))


## --- Calculate erosion rate at each coordinate point of interest
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("data/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("data/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point
    rockslope = avg_over_area(srtm15_slope, rocklat, rocklon, srtm15_sf, halfwidth=7)

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)


## --- Save data to a file so we can access it easily without rerunning slow scripts
    filecols = ["sed" "met" "ign" "cover" "lat" "lon" "slope"]
    writedlm("output/responses_parsed.tsv", vcat(filecols, hcat(sed, met, ign, cover, rocklat, rocklon, rockslope)))


## --- Do statistics to get sum (s), mean (m), and standard error (e)
    ersn_sed_s, ersn_sed_m, ersn_sed_e = get_stats(rock_ersn[sed])                # Sedimentary
    ersn_ign_s, ersn_ign_m, ersn_ign_e = get_stats(rock_ersn[ign])                # Igneous
    ersn_met_s, ersn_met_m, ersn_met_e = get_stats(rock_ersn[met])                # Metamorphic
    ersn_cryst_s, ersn_cryst_m, ersn_cryst_e = get_stats(rock_ersn[cryst])        # Crystalline
    ersn_global_s, ersn_global_m, ersn_global_e = get_stats(rock_ersn[matched])   # Global (matched points only)

    # Calculate relative contributions of each rock type to total global erosion
    ersn_rel_sed = ersn_sed_s / ersn_global_s * 100
    ersn_rel_ign = ersn_ign_s / ersn_global_s * 100
    ersn_rel_met = ersn_met_s / ersn_global_s * 100
    ersn_rel_cryst = ersn_cryst_s / ersn_global_s * 100

    # Print to terminal
    @info "Erosion rates by rock type
    global = $(round(ersn_global_m, digits=2)) mm/kyr
    sed = $(round(ersn_sed_m, digits=2)) mm/kyr
    ign = $(round(ersn_ign_m, digits=2)) mm/kyr
    met = $(round(ersn_met_m, digits=2)) mm/kyr
    cryst (ign + (met & ~sed)) = $(round(ersn_cryst_m, digits=2)) mm/kyr"

    @info "Total global erosion by rock type
    global total = $(round(ersn_global_s, sigdigits=3)) mm/kyr
    sed = $(round(ersn_sed_s, sigdigits=3)) mm/kyr ($(round(ersn_rel_sed, digits=2))%)
    ign = $(round(ersn_ign_s, sigdigits=3)) mm/kyr ($(round(ersn_rel_ign, digits=2))%)
    met = $(round(ersn_met_s, sigdigits=3)) mm/kyr ($(round(ersn_rel_met, digits=2))%)
    cryst (ign + (met & ~sed)) = $(round(ersn_cryst_s, sigdigits=3)) mm/kyr ($(round(ersn_rel_cryst, digits=2))%)\n"


## --- Import and parse the bulk EarthChem data into a NamedTuple
    bulk_file = matopen("data/bulk.mat")
    bulk_dict = read(bulk_file, "bulk")
    close(bulk_file)

    bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

    # Match EarthChem Types to rock type names
    bulksed, bulkign, bulkmet = match_earthchem(bulk.Type, major=true)

    bulkcryst = bulkign .| bulkmet                            # Crystalline (ign and met where met excludes metaseds and metaigns)
    bulk_not_matched = .!bulksed .& .!bulkign .& .!bulkmet    # Unmatched samples
    bulk_matched = .!bulk_not_matched                         # All matched samples

    # Print to terminal
    @info "not matched = $(count(bulk_not_matched)) of $(length(bulk.Type)) total ($(round((count(bulk_not_matched))/(length(bulk.Type))*100, digits=2))%)\n"


## --- Find the average combined phosphorus content of each rock type
    # # Average P content (corrected from ppm to wt.%) per type
    # p_sed = nanmean(bulk.P[bulksed]) * 1e-6
    # p_ign = nanmean(bulk.P[bulkign]) * 1e-6
    # p_met = nanmean(bulk.P[bulkmet]) * 1e-6
    # p_cryst = nanmean(bulk.P[bulkcryst]) * 1e-6

    # Calculate wt.% as only P2O5
    pwt_sed = nanmean(bulk.P2O5[bulksed])
    pwt_ign = nanmean(bulk.P2O5[bulkign])
    pwt_met = nanmean(bulk.P2O5[bulkmet])
    pwt_cryst = nanmean(bulk.P2O5[bulkcryst])

    prel_sed = pwt_sed / (pwt_sed+pwt_ign+pwt_met) * 100
    prel_ign = pwt_ign / (pwt_sed+pwt_ign+pwt_met) * 100
    prel_met = pwt_met / (pwt_sed+pwt_ign+pwt_met) * 100
    prel_cryst = pwt_cryst / (pwt_sed+pwt_ign+pwt_met) * 100

    # Print to terminal
    @info "Average phosphorus content by rock type:
    sed = $(round(pwt_sed,digits=2)) wt.%
    ign = $(round(pwt_ign,digits=2)) wt.%
    met = $(round(pwt_met,digits=2)) wt.%
    cryst (ign + (met & ~sed)) = $(round(pwt_cryst,digits=2)) wt.%"

    @info "Relative abundance of phosphorus by rock type:
    sed = $(round(prel_sed,digits=2))%
    ign = $(round(prel_ign,digits=2))%
    met = $(round(prel_met,digits=2))%
    cryst (ign + (met & ~sed)) = $(round(prel_cryst,digits=2))%\n"


## -- Calculate what % of total eroded P comes from each rock type
    # Constants
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    const crustal_density = 2750              # Average crustal density (kg/m³)

    # Crustal areas (m^2)
    sed_area = contl_area * sed_area_frac
    ign_area = contl_area * ign_area_frac
    met_area = contl_area * met_area_frac
    cryst_area = contl_area * cryst_area_frac

    # Calculate P flux by source contributions in kg/yr
    sed_contrib = ersn_sed_m * pwt_sed/100 * sed_area * crustal_density * 1e-6
    ign_contrib = ersn_ign_m * pwt_ign/100 * ign_area * crustal_density * 1e-6
    met_contrib = ersn_met_m * pwt_met/100 * met_area * crustal_density * 1e-6
    cryst_contrib = ersn_cryst_m * pwt_cryst/100 * cryst_area * crustal_density * 1e-6
    global_contrib = sed_contrib + ign_contrib + met_contrib

    # Print to terminal
    @info "Global P flux by source:
    sed = $(round(sed_contrib, sigdigits=3)) kg/yr
    ign = $(round(ign_contrib, sigdigits=3)) kg/yr
    met = $(round(met_contrib, sigdigits=3)) kg/yr
    cryst = $(round(cryst_contrib, sigdigits=3)) kg/yr
    global = $(round(global_contrib, sigdigits=3))\n"


## -- Calculate P content for a more granular rock type catagorization
    bulksed, bulkign, bulkmet, bulkvolc, bulkplut, bulkmetaign, bulkmetased = match_earthchem(bulk.Type, major=false)

    pwt_sed = nanmean(bulk.P2O5[bulksed])
    pwt_ign = nanmean(bulk.P2O5[bulkign])
    pwt_met = nanmean(bulk.P2O5[bulkmet])

    pwt_volc = nanmean(bulk.P2O5[bulkvolc])
    pwt_plut = nanmean(bulk.P2O5[bulkplut])
    pwt_metaign = nanmean(bulk.P2O5[bulkmetaign])
    pwt_metased = nanmean(bulk.P2O5[bulkmetased])
  
## -- End of file