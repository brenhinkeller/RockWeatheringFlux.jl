## --- Setup
  # External packages
  using ProgressMeter: @showprogress
  using StatGeochem
  using DelimitedFiles
  using Plots
  using Dates
  using LoopVectorization

	# File parsing packages
  using JLD
  using HDF5
  using HTTP
  using JSON
  using MAT

  # Local utilities
  include("Utilities.jl")

#=
## --- Generate some random points and get their lithology from Macrostrat / Burwell API
  # Generate random points on the continental crust
  npoints = 5
	etopo = get_etopo("elevation")
	rocklat, rocklon, elevations = gen_continental_points(npoints, etopo)

  # Initialize
  zoom = 11
	savefilename = "responses3"
  responses = Array{Any}(undef, npoints, 1)

  # Request data from API
  @showprogress 5 for i = 1:npoints
    try
        responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
    catch
      @warn "No response from Macrostrat server for coordinate $i/$npoints. Trying again in 5 seconds. \n"
      try
        # Wait and try again
        sleep(5)
        responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
      catch
        # If still nothing, add warning
        responses[i] = "No response"
        @warn "No data from Macrostrat server for coordinate $i/$npoints\n"
      end
    end
    sleep(0.05)

		# Checkpoint save every 10,000 points
    if mod(i,10000)==0
        save("data/$savefilename.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)
    end
  end

	# Save the file
  save("data/$savefilename.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)
=#

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
  sed .&= .!cover
  ign .&= .!cover
  met .&= .!cover .& .!sed &. .!ign   # Metaseds / metaigns grouped with sed / ign
  cryst = ign .| met                  # Define crystalline rocks as igneous or non metased/ign metamorphic

  # Figure out how many data points weren't matched
  not_matched = .~(sed .| ign .| met .| cover);
  multi_matched = (sed .& ign) .| (sed .& met) .| (ign .& met)

  number_not_matched = sum(not_matched)
  number_multi_matched = sum(multi_matched)
  total_matched = npoints - number_not_matched

  # Relative abundance, significant figures are also completely arbitrary
  sed_abnce = round(sum(sed) / total_matched*100, digits=2)
  ign_abnce = round(sum(ign) / total_matched*100, digits=2)
  met_abnce = round(sum(met) / total_matched*100, digits=2)
  cryst_abnce = round(sum(cryst) / total_matched*100, digits=2)

  # Output data to terminal
  @info "not matched = $number_not_matched, conflicting matches = $number_multi_matched\n"

  @info "Rock type totals and relative abundance (multi matched may be counted twice):
  sed = $(sum(sed)) ($(sed_abnce))
  ign = $(sum(ign)) ($(ign_abnce))
  met = $(sum(met)) ($(met_abnce))
  cryst (ign + (met & ~sed)) = $(sum(cryst)) ($(cryst_abnce))\n"


## -- Write the data to .tsv files so we can access or export it
  writedlm("output/notmatched.tsv", hcat(rocktype[not_matched], rockname[not_matched], rockdescrip[not_matched], 
    rockstratname[not_matched], rockcomments[not_matched]))
  writedlm("output/multimatched.tsv", hcat(rocktype[multi_matched], rockname[multi_matched], rockdescrip[multi_matched], 
    rockstratname[multi_matched], rockcomments[multi_matched]))
  writedlm("output/ignsed.tsv", hcat(rocktype[ign .& sed], rockname[ign .& sed], rockdescrip[ign .& sed], 
    rockstratname[ign .& sed], rockcomments[ign .& sed]))
  writedlm("output/plutonic.tsv", hcat(rocktype[plut], rockname[plut], rockdescrip[plut], 
    rockstratname[plut], rockcomments[plut]))
	writedlm("output/intrusive.tsv", hcat(rocktype[hypabyssal .& .~plut], rockname[hypabyssal .& .~plut], rockdescrip[hypabyssal .& .~plut], 
    rockstratname[hypabyssal .& .~plut], rockcomments[hypabyssal .& .~plut]))

	t = elevations .> 4000
	writedlm("output/highelev.tsv", hcat(rocktype[t], rockname[t], rockdescrip[t], rockstratname[t], rockcomments[t]))


## --- Calculate erosion rate at each coordinate point of interest
  # Load the slope variable from the SRTM15+ maxslope file
  srtm15_slope = h5read("data/srtm15plus_maxslope.h5", "vars/slope")
  srtm15_sf = h5read("data/srtm15plus_maxslope.h5", "vars/scalefactor")

  # Get slope at each point (rocklat, rocklon)
  rockslope = avg_over_area(srtm15_slope, rocklat, rocklon, srtm15_sf, halfwidth=7)

  # All erosion rates (mm/kyr)
  rock_ersn = emmkyr.(rockslope)

  # Do statistics
  ersn_sed_s, ersn_sed_m, ersn_sed_e = get_stats(rock_ersn[sed])            # Sedimentary
  ersn_sed_s, ersn_ign_m, ersn_ign_e = get_stats(rock_ersn[ign])            # Igneous
  ersn_met_s, ersn_met_m, ersn_met_e = get_stats(rock_ersn[met])            # Metamorphic
  ersn_cryst_s, ersn_cryst_m, ersn_cryst_e = get_stats(rock_ersn[cryst])    # All crystalline (igneous and metamorphic)
  ersn_global_s, ersn_global_m, ersn_global_e = get_stats(rock_ersn)        # Global

  # Calculate relative contributions of each rock type to total global erosion
  ersn_rel_sed = ersn_sed_s / ersn_global_s * 100
  ersn_rel_ign = ersn_ign_s / ersn_global_s * 100
  ersn_rel_met = ersn_met_s / ersn_global_s * 100
  ersn_rel_cryst = ersn_cryst_s / ersn_global_s * 100

  # Print to terminal (more arbitrary significant figure use)
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

  # Interpret the rock type codes using the list in RockNameInference.m
  n_bulk = length(bulk.Type)
  bulksed = vec(fill(false, n_bulk))
  bulkign = vec(fill(false, n_bulk))
  bulkmet = vec(fill(false, n_bulk))
  unmatched = vec(fill(false, n_bulk))

  @inbounds for i in eachindex(bulk.Type)
    code = modf(bulk.Type[i])[2]

    if code==1
      bulksed[i] = true
    elseif code==2
      bulkign[i] = true
    elseif code==3
      bulkmet[i] = true
    else
      unmatched[i] = true
    end
  end

  n_unmatched = count(unmatched)
  @info "$n_unmatched of $n_bulk total samples ($(round(n_unmatched/n_bulk*100, digits=2))%) were not matched to a rock type.\n"


## --- Find the average combined P and P2O5 content of each rock type
  p2o5_sed = nanmean(bulk.P2O5[bulksed])
  p2o5_ign = nanmean(bulk.P2O5[bulkign])
  p2o5_met = nanmean(bulk.P2O5[bulkmet])
  
  p2o5_cryst = nanmean(bulk.P2O5[bulkmet]) #problem
  p_cryst = p_ign + p_met #problem

  p2o5_total = p2o5_sed + p2o5_ign + p2o5_met
  

  p_sed = nanmean(bulk.P[bulksed]) * 1e-6
  p_ign = nanmean(bulk.P[bulkign]) * 1e-6
  p_met = nanmean(bulk.P[bulkmet]) * 1e-6
  p_total = p_sed + p_ign + p_met

  # Calculate combined wt.% and relative abundance of P + P2O5 by rock type
  pwt_sed = p2o5_sed+p_sed
  pwt_ign = p2o5_ign+p_ign
  pwt_met = p2o5_met+p_met
  pwt_cryst = p2o5_cryst+p_cryst

  prel_sed = pwt_sed / (p2o5_total+p_total) * 100
  prel_ign = pwt_ign / (p2o5_total+p_total) * 100
  prel_met = pwt_met / (p2o5_total+p_total) * 100
  prel_cryst = pwt_cryst / (p2o5_total+p_total) * 100

  # Print to terminal
  @info "Average combined P2O5 and P content by rock type:
  sed = $(round(pwt_sed,digits=2)) wt.%
  ign = $(round(pwt_ign,digits=2)) wt.%
  met = $(round(pwt_met,digits=2)) wt.%
  cryst (ign + (met & ~sed)) = $(round(pwt_cryst,digits=2)) wt.%"

  @info "Relative abundance of combined P2O5 and P by rock type
  sed = $(round(prel_sed,digits=2))%
  ign = $(round(prel_ign,digits=2))%
  met = $(round(prel_met,digits=2))%
  cryst (ign + (met & ~sed)) = $(round(prel_cryst,digits=2))%\n"


## -- Calculate what % of total eroded P comes from each rock type
  # get P into  kg/yr and compare to pre-existing estimates
  # ersn_sed_m * pwt_sed # m/Myr * wt% # If you multiply by area of the continents in m2 then you get m3/myr, multiply by density (say 2750 for average crust) then you have kg/myr, then can convert that to kg P/myr
  # ersn_ign_m * pwt_ign
  # ersn_met_m * pwt_met

## -- End of file