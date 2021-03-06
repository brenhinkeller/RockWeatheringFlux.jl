	1+1
## --- Setup

    # External packages
    using ProgressMeter: @showprogress
    using Plots
    using JLD, HDF5
    try
        using StatGeochem
    catch
        Pkg.clone("https://github.com/brenhinkeller/StatGeochem.jl")
        using StatGeochem
    end

    # Local utilities
    using HTTP, JSON
    include("RockWeatheringFluxUtilities.jl")

## --- Generate some random points on a sphere

    npoints = 50000;

    rocklat = Array{Float64}(0)
    rocklon = Array{Float64}(0)
    etopo = get_etopo("elevation")
    while length(rocklat) < npoints

        # Generate some random latitudes and longitudes with uniform
        #  spatial density on the globe
        (randlat, randlon) = random_lat_lon(npoints)

        # Find which points are above sea level
        elevations = find_etopoelev(etopo,randlat,randlon)
        abovesea = elevations .> 0

        # Concatenate together all the points that represent exposed crust
        rocklat = vcat(rocklat,randlat[abovesea])
        rocklon = vcat(rocklon,randlon[abovesea])
    end

    rocklat = rocklat[1:npoints]
    rocklon = rocklon[1:npoints]
    elevations = find_etopoelev(etopo,rocklat,rocklon)


## -- Get data from Macrostrat / Burwell API

    zoom = 11
	savefilename = "responses2"
    responses = Array{Any}(length(elevations))
    @showprogress 5 for i = 1:length(elevations)
        try
            responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
        catch
            print("Warning: no data from Macrostrat server \n")
            try
                # Wait and try again
                sleep(10)
                responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
            catch
                # if still nothing, add warnin
                responses[i] = "No response"
                print("Warning: no response from Macrostrat server\n")
            end
        end
        sleep(0.05)
        if mod(i,10000)==0
            save("data/$savefilename.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)
        end
    end

    save("data/$savefilename.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)

## --- Load from saved version (if applicable)

    retrive_file = load("data/responses.jld")
    elevations = retrive_file["elevations"]
    rocklat = retrive_file["latitude"]
    rocklon = retrive_file["longitude"]
    responses = retrive_file["responses"]
    npoints = retrive_file["npoints"]

## ---
	retrive_file = load("data/responses2.jld")
	elevations = vcat(elevations, retrive_file["elevations"])
	rocklat = vcat(rocklat, retrive_file["latitude"])
	rocklon = vcat(rocklon, retrive_file["longitude"])
	responses = vcat(responses, retrive_file["responses"])
	npoints = npoints + retrive_file["npoints"]

## --- Parse the macrostrat responses

    # Variables for parsed responses
    rocktype = Array{String}(length(elevations))
    rockdescrip = Array{String}(length(elevations))
    rockname = Array{String}(length(elevations))
    rockstratname = Array{String}(length(elevations))
    rockcomments = Array{String}(length(elevations))
	agemax = Array{Float64}(length(elevations))
	agemin = Array{Float64}(length(elevations))

    # Parse saved responses
    for i = 1:length(elevations)
        rocktype[i] = get_macrostrat_lith(responses[i])
        rockdescrip[i] = get_macrostrat_descrip(responses[i])
        rockname[i] = get_macrostrat_name(responses[i])
        rockstratname[i] = get_macrostrat_strat_name(responses[i])
        rockcomments[i] = get_macrostrat_comments(responses[i])
		agemax[i] = get_macrostrat_max_age(responses[i])
		agemin[i] = get_macrostrat_min_age(responses[i])
    end

    # Convert to lowercase to match names of rock types
    rocktype = lowercase.(rocktype)
    rockdescrip = lowercase.(rockdescrip)
    rockname = lowercase.(rockname)
    rockstratname = lowercase.(rockstratname)
    rockcomments = lowercase.(rockcomments)

    # Replacing tabs with spaces so as not to be confused with the delimitator if exported
    rocktype = replace.(rocktype, "    ", " ")
    rockdescrip = replace.(rockdescrip, "    ", " ")
    rockname = replace.(rockname, "    ", " ")
    rockstratname = replace.(rockstratname, "    ", " ")
    rockcomments = replace.(rockcomments, "    ", " ")

    #= Use this link to check the information on certain points in the macrostrat
    url = "https://macrostrat.org/api/mobile/map_query?lat=$(rocklat[end])&lng=$(rocklon[end])&z=11" =#

## ---

	# Names or partial names of different rock types (from GetBurwellBulkAge.m)

	covertypes = ["cover", "unconsolidated", "quaternary", "lluv", "soil", "regolith", "laterite", "surficial deposits", "talus", "scree", "mass-wasting", "slide", "peat", "swamp", "marsh", "water", "ice", "glaci", "till", "loess", "gravel"] # "debris",
    sedtypes = ["sediment", "fluv", " clast", "siliciclast", "conglomerat", "gravel", "sand", "psamm", "arenit", "arkos", "silt", "mud", "marl", "clay", "shale", "wacke", "argillite", "argillaceous", "pelit", "pebble", "carbonate", "limestone", "dolo", "caliche", "chalk", "travertine", "tavertine", "teravertine", "tufa", "evaporite", " salt", "salt flat", "gypsum", "boulder", "diamict", "tillite", "stream", "beach", "terrace", "chert", "banded iron", "coal", "anthracite", "marine deposits", "turbidite", "flysch", "paleosol"]

    volctypes = ["volcan", "lava", "lahar", "ignimbrite", "ashfall", "tuff", "diatreme", "pipe", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", "latite", "basanite", "phonolite", "fonolito", "trachyte", "palagonite", "mugearite", "kimberlite", "ultramafitite", "komatiite",]
	pluttypes = ["pluton", "batholith", "granit", "tonalit", "gabbro", "norite", "diorit", "monzonit", "syenit", "peridot", "dunit", "harzburg", "anorthosite", "mangerite", "charnockite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", "porphyry", "megacryst", "rapakivi", "bronzitite", "alaskite", "troctolite",]
	hypabyssaltypes = ["intrus", "hypabyssal", "sill", "dike", "stock", "laccolith", "lopolith", "dolerit", "diabase", "porphyry", "microgranite"]

	igntypes = vcat(["igneous", "silicic ", "mafic", "felsic", "basite",],volctypes,pluttypes,hypabyssaltypes)

	metasedtypes = ["para", "metased", "meta-sed", "schist", "quartzite", "marble", "skarn", "slate", "phyllite",]
	metaigntypes = ["ortho", "metaign", "meta-ign", "serpentin", "amphibolit", "greenstone", "eclogite", "metabasite", "migma",]
    mettypes = vcat(metasedtypes, metaigntypes, ["gneiss", "granulit", "hornfels", "granofels", "mylonit", "meta", "cataclasite", "melange", "gouge", "tecton", "calc silicate"])

	lowgradetypes = ["slate", "phyllite", "serpentin", "greenstone", "greenschist", "zeolite", "gossan", "alter", "hydrothermal", "palagonite",]
	highgradetypes = ["crystalline", "basement", "marble", "skarn", "blueschist", "gneiss", "amphibolit", "eclogite", "granulit", "hornfels", "granofels", "sanidinite", "migma", "enderbite", "anorthosite", "charnockite", "pyroxenite", "peridot", "dunit", "harzburg", "high grade metamorphic"]
	cataclastictypes = ["mylonit", "cataclasite", "melange", "gouge", "tecton",]

    # Allocate arrays for each sample for each rock type
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


    # Check which burwell "lith" rocktypes match one of the rock types
    # Try the "major:: {...}" type first, if present
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
    multi_matched = (sed .& ign) .| (sed .& met) .| (ign .& met);

    number_not_matched = sum(not_matched)
    number_multi_matched = sum(multi_matched)

    print("not matched = $number_not_matched, conflicting matches = $number_multi_matched\n")

## ---

    # Print proportions of sed vs ign, met, and nonsed rocks in responses
    print("sed = ", sum(sed), ", ign + (met & ~sed)= ", sum(ign .| (met .& .~sed)))

    # Create a file to check matching errors
    writedlm("notmatched.tsv", hcat(rocktype[not_matched], rockname[not_matched], rockdescrip[not_matched], rockstratname[not_matched], rockcomments[not_matched]))

    writedlm("multimatched.tsv", hcat(rocktype[multi_matched], rockname[multi_matched], rockdescrip[multi_matched], rockstratname[multi_matched], rockcomments[multi_matched]))

    writedlm("ignsed.tsv", hcat(rocktype[ign .& sed], rockname[ign .& sed], rockdescrip[ign .& sed], rockstratname[ign .& sed], rockcomments[ign .& sed]))

	writedlm("plutonic.tsv", hcat(rocktype[plut], rockname[plut], rockdescrip[plut], rockstratname[plut], rockcomments[plut]))

	writedlm("intrusive.tsv", hcat(rocktype[hypabyssal .& .~plut], rockname[hypabyssal .& .~plut], rockdescrip[hypabyssal .& .~plut], rockstratname[hypabyssal .& .~plut], rockcomments[hypabyssal .& .~plut]))

	t = elevations .> 4000
	writedlm("highelev.tsv", hcat(rocktype[t], rockname[t], rockdescrip[t], rockstratname[t], rockcomments[t]))

## --- Kernel function for resampling
	function unif_norm(rng, mu, sigma)
		return mu + (2rand(rng)-1)*sigma[1] + randn(rng)*sigma[2]
	end

## --- Plot abundance of sed, highgrade+plutonic and lowgrade+volcanic

	nresampled = 10^6
	edges = 0:40:4000
	centers = cntr(edges)
	nbins = length(centers)

	t = sed .& .~cover
	age = (agemax[t] + agemin[t])/2
	# sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	sigma = (agemax[t] - agemin[t])/2
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkblue), linewidth=0.1, label="")
	plot!(h,ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true, xlabel="Age (Ma)")
	savefig(h,"sed.pdf")

	t = ign .& .~(cover .| sed)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkred), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"ign.pdf")

	t = (volc .| lowgrade) .& .~(cover .| sed)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:red), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"volc.pdf")

	t = hypabyssal .& .~(cover .| sed)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkred), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"shallow intrusive.pdf")

	t = plut .& .~(cover .| sed)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkred), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"plut.pdf")

	t = highgrade .& .~(cover .| sed)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkred), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"highgrade.pdf")

	t = (plut .| highgrade) .& .~(cover .| sed)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkred), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"basement (plut+highgrade).pdf")

	t = .~(cover .| sed)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkred), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"nonsed.pdf")

	t = .~(cover .| sed .| metased .| volc .| hypabyssal .| not_matched)
	age = (agemax[t] + agemin[t])/2
	sigma = (agemax[t] - agemin[t])/2/sqrt(2*log(2))
	# ageresampled = bsresample(age,sigma,nresampled,kernel=uniform)
	ageresampled = bsresample(age,
		collect(zip(sigma,age*0.025)),
		nresampled,
		kernel=unif_norm
	)
	ageresampled = ageresampled[(ageresampled .> 0) .& (ageresampled .< 4500)]
	h = histogram(ageresampled, bins=edges, fill=(0,:darkred), linewidth=0.1, label="")
	plot!(h, ylabel="N", ylims=(0,nresampled/nbins*6), xlims=(0,4000), xflip=true,  xlabel="Age (Ma)")
	savefig(h,"nonsednonvolc.pdf")


## --- get references

	authors = get_macrostrat_ref_authors.(responses)
	years = get_macrostrat_ref_year.(responses)
	titles = get_macrostrat_ref_title.(responses)
	dois = get_macrostrat_ref_doi.(responses)

	refstrings = authors .* "|" .* years .* "|" .* titles .* "|" .* dois
	writedlm("burwellrefs.tsv",unique(refstrings),"\t")

## --- Find slope and erosion rate

    function Emmkyr(slp)
        return 10^(slp*0.00567517 + 0.971075)
    end

    slope = get_srtm15plus_aveslope("slope")
    rockslope = find_srtm15plus_aveslope_around(slope,rocklat,rocklon, halfwidth=7, max_allowed_slope=0xffff)

	# Plot slope-area histogram
	h = histogram(rockslope, label="",xlims=[0,500])
	plot!(h,xlabel="Slope (m/km)",ylabel="N")
	savefig(h,"Slope_Histogram.pdf")
	display(h)

    rockEmmkyr = Emmkyr.(rockslope)

    sedEsum = nansum(rockEmmkyr[sed])
    crystEsum = nansum(rockEmmkyr[cryst])

    sedEmean = nanmean(rockEmmkyr[sed])
    crystEmean = nanmean(rockEmmkyr[cryst])
	Emean = nanmean(rockEmmkyr)

    propotionSedArea = sum(sed)/(sum(sed)+sum(ign .| (met .& .~sed)));
    print("$(round(100*propotionSedArea,3))% sed area, $(round(100*(1-propotionSedArea),3))% crystalline area\n")
    print("sed E sum: $sedEsum cryst E sum: $crystEsum\n")
    print("$(round(100*sedEsum/(sedEsum+crystEsum),3))% sed flux, $(round(100*crystEsum/(sedEsum+crystEsum),3))% crystalline flux\n")
    print("sed mean E: $(round(sedEmean,3)), crystalline mean E: $(round(crystEmean,3)), overall mean: $(round(Emean,3))\n")

## --- Proportion vs slope

nbins = 50
xmin = 0
xmax = 500

prop_cryst = Array{Float64}(nbins)
prop_sed = Array{Float64}(nbins)
binedges = range(xmin,xmax,length=nbins+1)
bincenters = cntr(binedges)
for i=1:nbins
	t = (rockslope.>binedges[i]) .& (rockslope.<binedges[i+1])
	prop_cryst[i] = sum(cryst[t])/sum(cryst[t] .| sed[t])
	prop_sed[i] = sum(sed[t])/sum(cryst[t] .| sed[t])
end

h = plot(bincenters,prop_sed,label = "", xlabel="Slope (m/km)", ylabel="Proportion of Sedimentary Bedrock")
# plot!(h,bincenters,prop_cryst+prop_sed,label="Sedimentary Basement")
plot!(legend=:topleft)
savefig(h,"Proportion_sedimentary_bedrock_slope.pdf")
display(h)

## --- Proportion vs elevation

nbins = 50
xmin = 0
xmax = 5000

prop_cryst = Array{Float64}(nbins)
prop_sed = Array{Float64}(nbins)
binedges = range(xmin,xmax,length=nbins+1)
bincenters = cntr(binedges)
for i=1:nbins
	t = (elevations.>binedges[i]) .& (elevations.<binedges[i+1])
	prop_cryst[i] = sum(cryst[t])/sum(cryst[t] .| sed[t])
	prop_sed[i] = sum(sed[t])/sum(cryst[t] .| sed[t])
end

h = plot(bincenters,prop_sed,label = "", xlabel="Elevation (m)", ylabel="Proportion of Sedimentary Bedrock")
# plot!(h,bincenters,prop_cryst+prop_sed,label="Sedimentary Basement")
savefig(h,"Proportion_sedimentary_bedrock_elevation.pdf")
display(h)
## --- For comparison

    # Phanerozoic
    sedimentSubduction_km3_yr = 1.64 # Clift, 2009
    continentSedimentation_km3_yr = 0.8616 # Macrostrat
    continentalErosion_km3_yr = sedimentSubduction_km3_yr + sedimentSubduction_km3_yr

    continentalErosionQuaternary_km3_yr = 6.5 #\cite{Clift:2009fm}  / Milliman (1990)
    continentalArea_km2 =  1.4894E8 # km^2

    continentalErosion_km_yr = continentalErosion_km3_yr / continentalArea_km2
    continentalErosion_mm_kyr = continentalErosion_km_yr * 1000 * 1000 * 1000

### ---
    using StatGeochem
    elev = get_srtm15plus("elevation")

    function movmean33(x::Array{<:Any,2})
        # Simple moving average of 9-element box for 2-d array
        m = Array{Float64}(size(x))
        for i = 2:size(x,1)-1
            for j = 2:size(x,2)-1
                adj = [x[i,j],x[i+1,j+1],x[i,j+1],x[i-1,j+1],x[i-1,j],x[i-1,j-1],x[i,j-1],x[i+1,j-1],x[i+1,j]]
                m[i,j] = mean(adj)
            end
        end
        m[1,:] = m[2,:]
        m[end,:] = m[end-1,:]
        m[:,1] = m[:,2]
        m[:,end] = m[:,end-1]
        return m;
    end

    elev_smoothed = movmean33(elev)

    # Save results
    fid = h5open("srtm15plus_elev_smoothed.h5","w");
    g = g_create(fid, "vars");
    print("Saving to HDF5:\n")
    g["y_lat_cntr"] = y_lat_cntr
    g["x_lon_cntr"] = x_lon_cntr
    g["cellsize"] = cellsize
    g["scalefactor"] = scalefactor
    @time g["elevation","compress",3] = elev_smoothed;
    close(fid)

### ---
	elev = elev_smoothed
	inbasin = fill(false,size(elev))

	inbasin[elev .< -100] = true

	function fill_uphill(elev,inbasin,i,j,depth=0,tolerance=15)

        inbasin[i,j] = true
        depth += 1

        # If we're about to run in the recursion limit or the edge of the map, stop
        if depth > 10000 || i<2 || i > size(elev,1)-1 || j<2 || j > size(elev,2)-1
            return inbasin
        end

		# Otherwise, check each adjoining square recursively
		e = elev[i,j] - tolerance

		if elev[i+1,j+1] >= e && ~inbasin[i+1,j+1]
			inbasin = fill_uphill(elev,inbasin,i+1,j+1,depth)
		end
		if elev[i,j+1] >= e && ~inbasin[i,j+1]
			inbasin = fill_uphill(elev,inbasin,i,j+1,depth)
		end
		if elev[i-1,j+1] >= e && ~inbasin[i-1,j+1]
			inbasin = fill_uphill(elev,inbasin,i-1,j+1,depth)
		end
		if elev[i-1,j] >= e && ~inbasin[i-1,j]
			inbasin = fill_uphill(elev,inbasin,i-1,j,depth)
		end
		if elev[i-1,j-1] >= e && ~inbasin[i-1,j-1]
			inbasin = fill_uphill(elev,inbasin,i-1,j-1,depth)
		end
		if elev[i,j-1] >= e && ~inbasin[i,j-1]
			inbasin = fill_uphill(elev,inbasin,i,j-1,depth)
		end
		if elev[i+1,j-1] >= e && ~inbasin[i+1,j-1]
			inbasin = fill_uphill(elev,inbasin,i+1,j-1,depth)
		end
        if elev[i+1,j] >= e && ~inbasin[i+1,j]
            inbasin = fill_uphill(elev,inbasin,i+1,j,depth)
        end
		return inbasin
	end

	@showprogress for i=2:size(inbasin,1)-1
		for j=2:size(inbasin,2)-1
			if inbasin[i,j]
                # inbasin = fill_uphill(elev,inbasin,i,j)
				if ~inbasin[i+1,j] || ~inbasin[i+1,j+1] || ~inbasin[i,j+1] || ~inbasin[i-1,j+1] || ~inbasin[i-1,j] || ~inbasin[i-1,j-1]  || ~inbasin[i,j-1]  || ~inbasin[i+1,j-1]
					inbasin = fill_uphill(elev,inbasin,i,j)
				end
			end
		end
	end

    heatmap(downsample(inbasin,10))

    1+1
