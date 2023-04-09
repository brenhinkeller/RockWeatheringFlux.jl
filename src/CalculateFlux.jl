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


## --- Parse Macrostat data references
    # Preallocate
    authors = Array{String}(undef, npoints, 1)
    years = Array{String}(undef, npoints, 1)
    titles = Array{String}(undef, npoints, 1)
    dois = Array{String}(undef, npoints, 1)

    # Parse into preallocated arrays
    for i = eachindex(authors)
        authors[i] = get_macrostrat_ref_authors(responses[i])
        years[i] =get_macrostrat_ref_year(responses[i])
        titles[i] = get_macrostrat_ref_title(responses[i])
        dois[i] = get_macrostrat_ref_doi(responses[i])
    end

    # Write the references to a file
    refstrings = @. authors * " | " * years * " | " * titles * " | " * dois
    writedlm("output/burwellrefs.tsv",unique(refstrings),"\t")


## --- Match the Burwell rocktype with our rock types
    cats = match_rocktype(rocktype, rockname, rockdescrip, major=true)

    # Exclude suspected cover from the other three categories, just to be sure
    cats.sed .&= .! cats.cover .& .! cats.ign
    cats.ign .&= .! cats.cover

    # Metaseds / metaigns grouped with sed / ign
    cats.met .&= .! cats.cover .& .! cats.sed .& .! cats.ign
    
    # Define crystalline rocks as igneous or non metased/ign metamorphic
    cryst = cats.ign .| cats.met                               

    # Figure out how many data points weren't matched
    known = cats.sed .| cats.ign .| cats.met
    matched = known .| cats.cover
    not_matched = .!matched
    multi_matched = (cats.sed .& cats.ign) .| (cats.sed .& cats.met) .| (cats.ign .& cats.met)
    total_known = count(known)

    # Relative abundance
    sed_area_frac = count(cats.sed) / total_known
    ign_area_frac = count(cats.ign) / total_known
    met_area_frac = count(cats.met) / total_known
    cryst_area_frac = count(cryst) / total_known

    # Print to terminal
    @info "not matched = $(count(not_matched)), conflicting matches = $(count(multi_matched))\n"

    @info "Rock type totals and relative abundance (multi matched are counted twice):
    sed = $(count(cats.sed)) ($(round(sed_area_frac*100, digits=2))%)
    ign = $(count(cats.ign)) ($(round(ign_area_frac*100, digits=2))%)
    met = $(count(cats.met)) ($(round(met_area_frac*100, digits=2))%)
    cryst (ign + (met & ~sed)) = $(count(cryst)) ($(round(cryst_area_frac*100, digits=2))%)\n"


## -- Write the data to .tsv files so we can access or export it
    writedlm("output/notmatched.tsv", hcat(rocklat[not_matched], rocklon[not_matched], 
        rocktype[not_matched], rockname[not_matched], rockdescrip[not_matched], rockstratname[not_matched], 
        rockcomments[not_matched]))
    writedlm("output/multimatched.tsv", hcat(rocklat[multi_matched], rocklon[multi_matched], 
        rocktype[multi_matched], rockname[multi_matched], rockdescrip[multi_matched], 
        rockstratname[multi_matched], rockcomments[multi_matched]))
    writedlm("output/sed.tsv", hcat(rocklat[cats.sed], rocklon[cats.sed], rocktype[cats.sed], 
        rockname[cats.sed], rockdescrip[cats.sed], rockstratname[cats.sed], rockcomments[cats.sed]))
    writedlm("output/ign.tsv", hcat(rocklat[cats.ign], rocklon[cats.ign], rocktype[cats.ign], 
        rockname[cats.ign], rockdescrip[cats.ign], rockstratname[cats.ign], rockcomments[cats.ign]))
    writedlm("output/met.tsv", hcat(rocklat[cats.met], rocklon[cats.met], rocktype[cats.met], 
        rockname[cats.met], rockdescrip[cats.met], rockstratname[cats.met], rockcomments[cats.met]))

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
    writedlm("output/responses_parsed.tsv", vcat(filecols, 
        hcat(cats.sed, cats.met, cats.ign, cats.cover, rocklat, rocklon, rockslope)))


## --- Do statistics to get sum (s), mean (m), and standard error (e)
    ersn_sed_s, ersn_sed_m, ersn_sed_e = get_stats(rock_ersn[cats.sed])           # Sedimentary
    ersn_ign_s, ersn_ign_m, ersn_ign_e = get_stats(rock_ersn[cats.ign])           # Igneous
    ersn_met_s, ersn_met_m, ersn_met_e = get_stats(rock_ersn[cats.met])           # Metamorphic
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
    earthchem_raw = matopen("data/bulk.mat")
    earthchem_dict = read(earthchem_raw, "bulk")
    close(earthchem_raw)

    earthchem = NamedTuple{Tuple(Symbol.(keys(earthchem_dict)))}(values(earthchem_dict))

    # Match EarthChem Types to rock type names
    echem_cats = match_earthchem(earthchem.Type, major=true)

    # To match Macrostrat classification, crystalline rocks are ign and met 
    # Note that met excludes metaseds and metaigns
    echem_cryst = echem_cats.ign .| echem_cats.met

    echem_not_matched = .!echem_cats.sed .& .!echem_cats.ign .& .!echem_cats.met     # Unmatched samples
    echem_matched = .!echem_not_matched                                               # All matched samples

    # Print to terminal
    @info "not matched = $(count(echem_not_matched)) of $(length(earthchem.Type)) total ($(round((count(echem_not_matched))/(length(earthchem.Type))*100, digits=2))%)\n"


## --- Find the average combined phosphorus content of each rock type
    # # Average P content (corrected from ppm to wt.%) per type
    # p_sed = nanmean(bulk.P[bulksed]) * 1e-6
    # p_ign = nanmean(bulk.P[bulkign]) * 1e-6
    # p_met = nanmean(bulk.P[bulkmet]) * 1e-6
    # p_cryst = nanmean(bulk.P[bulkcryst]) * 1e-6

    # Calculate wt.% as only P2O5
    pwt_sed = nanmean(earthchem.P2O5[echem_cats.sed])
    pwt_ign = nanmean(earthchem.P2O5[echem_cats.ign])
    pwt_met = nanmean(earthchem.P2O5[echem_cats.met])
    pwt_cryst = nanmean(earthchem.P2O5[echem_cryst])
    pwt_total = nanmean(pwt_sed + pwt_ign + pwt_met)

    # Print to terminal
    @info "Average phosphorus content by rock type:
    sed = $(round(pwt_sed,digits=2)) wt.%
    ign = $(round(pwt_ign,digits=2)) wt.%
    met = $(round(pwt_met,digits=2)) wt.%
    cryst (ign + (met & ~sed)) = $(round(pwt_cryst,digits=2)) wt.%"

    # Calculate what fraction of phosphorus is in each rock type
    # Why do I do this? I don't think this is telling me anything
    prel_sed = pwt_sed / (pwt_total) * 100
    prel_ign = pwt_ign / (pwt_total) * 100
    prel_met = pwt_met / (pwt_total) * 100
    prel_cryst = pwt_cryst / (pwt_total) * 100
    
    # Print to terminal
    @info "Relative abundance of phosphorus by rock type:
    sed = $(round(prel_sed,digits=2))%
    ign = $(round(prel_ign,digits=2))%
    met = $(round(prel_met,digits=2))%
    cryst (ign + (met & ~sed)) = $(round(prel_cryst,digits=2))%\n"


## -- Calculate what % of total eroded P comes from each rock type
    # Constants
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    const crustal_density = 2750              # Average crustal density (kg/m³)

    # Crustal areas (m²)
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

## -- Get macrostrat rock types
    macro_cats = match_rocktype(rocktype, rockname, rockdescrip, major=false)

## -- Get erosion (m/Myr) for each rock type in Macrostrat
    macro_ersn = (
        siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], evaporite = [0.0], 
            coal = [0.0], sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
    )

    for i in eachindex(pflux_source)
        macro_ersn[i][1] = nanmean(rock_ersn[macro_cats[i]])
    end
    

## -- Get crustal area for each rock type in Macrostrat
    # Exclude suspected cover from the other three categories, just to be sure
    macro_cats.sed .&= .! macro_cats.cover .& .! macro_cats.ign
    macro_cats.ign .&= .! macro_cats.cover

    # Metaseds / metaigns grouped with sed / ign
    macro_cats.met .&= .! macro_cats.cover .& .! macro_cats.sed .& .! macro_cats.ign
    
    # Define crystalline rocks as igneous or non metased metamorphic
    macro_cryst = macro_cats.ign .| macro_cats.met

    # Figure out how many data points weren't matched
    known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    total_known = count(known_rocks)

    matched = known .| macro_cats.cover
    total_matched = count(matched)
    not_matched = .!matched
    multi_matched = ((macro_cats.sed .& macro_cats.ign) .| (macro_cats.sed .& macro_cats.met) 
        .| (macro_cats.ign .& macro_cats.met)
    )
    
    # Area (in m²) of each rock type. 
    # Using all rocks including cover because continental area includes areas of cover?
    crustal_area = (
        siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], evaporite = [0.0], 
            coal = [0.0], sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
        cover = [0.0]
    )
    for i in eachindex(crustal_area)
        crustal_area[i][1] = count(macro_cats[i]) / total_matched * contl_area
    end


## -- Get average sed contributions from EarthChem data
    echem_cats = match_earthchem(earthchem.Type, major=false)

    p_wt = (
        alluvium = [0.0], siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], 
            evaporite = [0.0], phosphorite = [0.0], coal = [0.0], volcaniclast = [0.0], 
            sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
    )
    for i in eachindex(p_wt)
        p_wt[i][1] = nanmean(earthchem.P2O5[echem_cats[i]]) 
    end

    # Calculate P flux by source contributions in kg/yr
    pflux_source = (
        siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], evaporite = [0.0], 
            coal = [0.0], sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
    )

    for i in eachindex(pflux_source)
        pflux_source[i][1] = macro_ersn[i] * p_wt[i] * crustal_area[i] * crustal_density * 1e-6
    end
    pflux_global = pflux_source.sed + pflux_source.ign + pflux_source.met
  
## -- End of file