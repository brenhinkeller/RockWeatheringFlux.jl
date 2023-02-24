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
  include("src/Utilities.jl")


## --- Generate some random points on the continental crust
  npoints = 5
	etopo = get_etopo("elevation")
	rocklat, rocklon, elevations = gen_continental_points(npoints, etopo)


## -- Get data from Macrostrat / Burwell API
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


## --- Load data that was just generated
#=
  I'm like 90% sure that we don't actually have to do this... we just finished generating the data. It's alrady
  in these arrays.
=#
#=
  @info "Loading generated data. This make take up to $(0.0001*npoints) minutes. Started $(Dates.format(now(), "HH:MM"))"
  retrive_file = load("data/$savefilename.jld")

  responses = retrive_file["responses"]
  elevations = float.(retrive_file["elevations"])
  
  rocklat = float.(retrive_file["latitude"])
  rocklon = float.(retrive_file["longitude"])    
  npoints = retrive_file["npoints"]
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


## --- Define rock types based on names or partial names of different rock types from GetBurwellBulkAge.m
  # Sedimentary
  sedtypes = ["sediment", "fluv", " clast", "siliciclast", "conglomerat", "gravel", "sand", "psamm", "arenit", "arkos", "silt", 
    "mud", "marl", "clay", "shale", "wacke", "argillite", "argillaceous", "pelit", "pebble", "carbonate", "limestone", "dolo", 
    "caliche", "chalk", "travertine", "tavertine", "teravertine", "tufa", "evaporite", " salt", "salt flat", "gypsum", "boulder", 
    "diamict", "tillite", "stream", "beach", "terrace", "chert", "banded iron", "coal", "anthracite", "marine deposits", "turbidite", 
    "flysch", "paleosol"]

  # Igneous
  volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme", "pipe", "basalt", "andesit", "dacit", "rhyolit", 
    "pillow", "carbonatite", "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", "fonolito", 
    "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", "komatiite",]
  pluttypes = ["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", "monzonit", "syenit", "peridot", "dunit", 
    "harzburg", "anorthosite", "mangerite", "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", 
    "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",]
  hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", "dolerit", "diabase", "porphyry", 
    "microgranite"]
  igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",],volctypes,pluttypes,hypabyssaltypes)

  # Metamorphic
  metasedtypes = ["para", "metased", "meta-sed", "schist", "quartzite", "marble", "skarn", "slate", "phyllite",]
  metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone", "eclogite", "metabasite", "migma",]
  mettypes = vcat(metasedtypes, metaigntypes, ["gneiss", "granulit", "hornfels", "granofels", "mylonit", "meta", "cataclasite", 
    "melange", "gouge", "tecton", "calc silicate"])
  lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", "gossan", "alter", "hydrothermal", 
    "palagonite",]
  highgradetypes = ["crystalline", "basement", "marble", "skarn", "blueschist", "gneiss", "amphibolit", "eclogite", "granulit", 
    "hornfels", "granofels", "sanidinite", "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", 
    "harzburg", "high grade metamorphic"]

  # Other
  covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", "laterite", "surficial deposits", "talus", 
      "scree", "mass-wasting", "slide", "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel", "debris"]
  cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]


## --- Match the Burwell rocktype with our rock types defined above
  # Preallocate arrays for each rock type
  sed = fill(false,npoints)
  ign = fill(false,npoints)
  met = fill(false,npoints)
	volc = fill(false,npoints)
	plut = fill(false,npoints)
	hypabyssal = fill(false,npoints)
	metaign = fill(false,npoints)
	metased = fill(false,npoints)
	lowgrade = fill(false,npoints)
	highgrade = fill(false,npoints)
	cover = fill(false,npoints)

  #=
    I am going to try to optimize this later. However, the first step is to assume that it works and move on. For the sake
    of time, it is not useful to me to figure out how this is matching and exactly what it is doing. As long as it does not
    take too long, this is good enough for now.

    Also this doesn't even run this slowly it just looks sad
  =#

  # Try the "major:{...}" type first
  for i = 1:length(sedtypes)
      sed = sed .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], sedtypes[i]) : false )
  end
  for i = 1:length(igntypes)
      ign = ign .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], igntypes[i]) : false )
  end
  for i = 1:length(mettypes)
      met = met .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], mettypes[i]) : false )
  end
  for i = 1:length(covertypes)
      cover = cover .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], covertypes[i]) : false )
  end
  for i = 1:length(volctypes)
      volc = volc .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], volctypes[i]) : false )
  end
  for i = 1:length(pluttypes)
      plut = plut .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], pluttypes[i]) : false )
  end
  for i = 1:length(metaigntypes)
      metaign = metaign .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], metaigntypes[i]) : false )
  end
  for i = 1:length(metasedtypes)
      metased = metased .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], metasedtypes[i]) : false )
  end
  for i = 1:length(lowgradetypes)
      lowgrade = lowgrade .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], lowgradetypes[i]) : false )
  end
  for i = 1:length(highgradetypes)
      highgrade = highgrade .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], highgradetypes[i]) : false )
  end
  for i = 1:length(hypabyssaltypes)
      hypabyssal = hypabyssal .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? containsi(x[1], hypabyssaltypes[i]) : false )
  end

  # Then check the rest of rocktype
  not_matched = .~(sed .| ign .| met .| cover)
  for i = 1:length(sedtypes)
    sed[not_matched] = sed[not_matched] .| containsi.(rocktype[not_matched],sedtypes[i])
  end
  for i = 1:length(igntypes)
    ign[not_matched] = ign[not_matched] .| containsi.(rocktype[not_matched],igntypes[i])
  end
  for i = 1:length(mettypes)
    met[not_matched] = met[not_matched] .| containsi.(rocktype[not_matched],mettypes[i])
  end
  for i = 1:length(covertypes)
    cover[not_matched] = cover[not_matched] .| containsi.(rocktype[not_matched],covertypes[i])
  end

	not_matched = .~(sed .| cover .| volc .| plut .| metaign .| metased .| lowgrade .| highgrade)
	for i = 1:length(volctypes)
      volc[not_matched] = volc[not_matched] .| containsi.(rocktype[not_matched],volctypes[i])
    end
	for i = 1:length(pluttypes)
      plut[not_matched] = plut[not_matched] .| containsi.(rocktype[not_matched],pluttypes[i])
    end
	for i = 1:length(metaigntypes)
	  metaign[not_matched] = metaign[not_matched] .| containsi.(rocktype[not_matched],metaigntypes[i])
	end
	for i = 1:length(metasedtypes)
	  metased[not_matched] = metased[not_matched] .| containsi.(rocktype[not_matched],metasedtypes[i])
	end
	for i = 1:length(lowgradetypes)
	  lowgrade[not_matched] = lowgrade[not_matched] .| containsi.(rocktype[not_matched],lowgradetypes[i])
	end
	for i = 1:length(highgradetypes)
	  highgrade[not_matched] = highgrade[not_matched] .| containsi.(rocktype[not_matched],highgradetypes[i])
	end
	for i = 1:length(hypabyssaltypes)
	  hypabyssal[not_matched] = hypabyssal[not_matched] .| containsi.(rocktype[not_matched],hypabyssaltypes[i])
	end


  # Then rockname
  not_matched = .~(sed .| ign .| met .| cover)
  for i = 1:length(sedtypes)
    sed[not_matched] = sed[not_matched] .| containsi.(rockname[not_matched],sedtypes[i])
  end
  for i = 1:length(igntypes)
    ign[not_matched] = ign[not_matched] .| containsi.(rockname[not_matched],igntypes[i])
  end
  for i = 1:length(mettypes)
    met[not_matched] = met[not_matched] .| containsi.(rockname[not_matched],mettypes[i])
  end
	for i = 1:length(covertypes)
      cover[not_matched] = cover[not_matched] .| containsi.(rockname[not_matched],covertypes[i])
  end

	not_matched = .~(sed .| cover .| volc .| plut .| metaign .| metased .| lowgrade .| highgrade)
	for i = 1:length(volctypes)
      volc[not_matched] = volc[not_matched] .| containsi.(rockname[not_matched],volctypes[i])
    end
	for i = 1:length(pluttypes)
      plut[not_matched] = plut[not_matched] .| containsi.(rockname[not_matched],pluttypes[i])
    end
	for i = 1:length(metaigntypes)
	  metaign[not_matched] = metaign[not_matched] .| containsi.(rockname[not_matched],metaigntypes[i])
	end
	for i = 1:length(metasedtypes)
	  metased[not_matched] = metased[not_matched] .| containsi.(rockname[not_matched],metasedtypes[i])
	end
	for i = 1:length(lowgradetypes)
	  lowgrade[not_matched] = lowgrade[not_matched] .| containsi.(rockname[not_matched],lowgradetypes[i])
	end
	for i = 1:length(highgradetypes)
	  highgrade[not_matched] = highgrade[not_matched] .| containsi.(rockname[not_matched],highgradetypes[i])
	end
	for i = 1:length(hypabyssaltypes)
	  hypabyssal[not_matched] = hypabyssal[not_matched] .| containsi.(rockname[not_matched],hypabyssaltypes[i])
	end


  # Then rockdescrip
  not_matched = .~(sed .| ign .| met .| cover)
  for i = 1:length(sedtypes)
    sed[not_matched] = sed[not_matched] .| containsi.(rockdescrip[not_matched],sedtypes[i])
  end
  for i = 1:length(igntypes)
    ign[not_matched] = ign[not_matched] .| containsi.(rockdescrip[not_matched],igntypes[i])
  end
  for i = 1:length(mettypes)
    met[not_matched] = met[not_matched] .| containsi.(rockdescrip[not_matched],mettypes[i])
  end
	for i = 1:length(covertypes)
	  cover[not_matched] = cover[not_matched] .| containsi.(rockdescrip[not_matched],covertypes[i])
	end

	not_matched = .~(sed .| cover .| volc .| plut .| metaign .| metased .| lowgrade .| highgrade)
	for i = 1:length(volctypes)
	  volc[not_matched] = volc[not_matched] .| containsi.(rockdescrip[not_matched],volctypes[i])
	end
	for i = 1:length(pluttypes)
	  plut[not_matched] = plut[not_matched] .| containsi.(rockdescrip[not_matched],pluttypes[i])
	end
	for i = 1:length(metaigntypes)
	  metaign[not_matched] = metaign[not_matched] .| containsi.(rockdescrip[not_matched],metaigntypes[i])
	end
	for i = 1:length(metasedtypes)
	  metased[not_matched] = metased[not_matched] .| containsi.(rockdescrip[not_matched],metasedtypes[i])
	end
	for i = 1:length(lowgradetypes)
	  lowgrade[not_matched] = lowgrade[not_matched] .| containsi.(rockdescrip[not_matched],lowgradetypes[i])
	end
	for i = 1:length(highgradetypes)
	  highgrade[not_matched] = highgrade[not_matched] .| containsi.(rockdescrip[not_matched],highgradetypes[i])
	end
	for i = 1:length(hypabyssaltypes)
	  hypabyssal[not_matched] = hypabyssal[not_matched] .| containsi.(rockdescrip[not_matched],hypabyssaltypes[i])
	end

	# Exclude suspected cover from the other three categories, just to be sure
	ign = ign .& .~cover
	met = met .& .~cover
	sed = sed .& .~cover
	cryst = ign .| (met .& .~sed)

  not_matched = .~(sed .| ign .| met .| cover);
  multi_matched = (sed .& ign) .| (sed .& met) .| (ign .& met)

  number_not_matched = sum(not_matched)
  number_multi_matched = sum(multi_matched)

  @info "not matched = $number_not_matched, conflicting matches = $number_multi_matched\n"


## --- Output results of matching to the terminal and .tsv files
  # Print information to terminal
  @info "sed = $(sum(sed)), ign + (met & !sed) = $(sum(ign .| (met .& .~sed)))"

  # Write the data to files so we can access or export it
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


## --- Calculate erosion rate at each point
  # Load the slope variable from the SRTM15+ maxslope file
  srtm15_slope = h5read("data/srtm15plus_maxslope.h5", "vars/slope")
  srtm15_sf = h5read("data/srtm15plus_maxslope.h5", "vars/scalefactor")

  # Get slope at each point (rocklat, rocklon)
  rockslope = avg_over_area(srtm15_slope, rocklat, rocklon, srtm15_sf, halfwidth=7)

  # Erosion rates (mm/kyr)
  rock_ersn = emmkyr.(rockslope)

  # Do statistics
  # After running the proof of concent, just make all the BitMatrices Bitvectors when they're assigned
  sed_ersn_sum = nansum(rock_ersn[vec(sed)])      # Sedimentary
  ign_ersn_sum = nansum(rock_ersn[vec(ign)])      # Igneous
  met_ersn_sum = nansum(rock_ersn[vec(met)])      # Metamorphic
  cryst_ersn_sum = nansum(rock_ersn[vec(cryst)])  # All crystalline (igneous and metamorphic)

  sed_ersn_mean = nanmean(rock_ersn[vec(sed)])      # Sedimentary
  ign_ersn_sum = nanmean(rock_ersn[vec(ign)])      # Igneous
  met_ersn_sum = nanmean(rock_ersn[vec(met)])      # Metamorphic
  cryst_ersn_sum = nanmean(rock_ersn[vec(cryst)])  # All crystalline (igneous and metamorphic)


## --- Import and parse the bulk EarthChem data
  #This currently doesn't work because I can't add MAT.jl to epidote (lol!)

  # Open and read file
  bulk_file = matopen("data/bulk.mat")
  bulk_dict = read(bulk_file, "bulk")
  close(bulk_file)

  # Parse into a NamedTuple
  bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

## --- 
