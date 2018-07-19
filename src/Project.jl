## --- Setup

    # External packages
    using ProgressMeter: @showprogress
    using Plots
    using JLD
    try
        using StatGeochem
    catch
        Pkg.clone("https://github.com/brenhinkeller/StatGeochem.jl")
        using StatGeochem
    end

    # Local utilities
    include("src/Utilities.jl")

## --- Generate some random points on a sphere

    npoints = 50000;
    (randlat, randlon) = random_lat_lon(npoints)

## --- Let's find the geology at one of these points

    etopo = get_etopo("elevation")
    elevations = find_etopoelev(etopo,randlat,randlon)

## --- Check which points are above sea level

    function abovesea(elevations,randlat,randlon)
        lat = Array{Float64}(size(randlat))
        lon = Array{Float64}(size(randlon))
        j = 0
        for i = 1:length(elevations)
            if elevations[i] > 0
                j += 1;
                lat[j] = randlat[i]
                lon[j] = randlon[i]
            end
        end
        return (lat[1:j], lon[1:j])
    end

    randlat = [];
    randlon = [];
    i = 0
    while i < npoints

        (temp_randlat, temp_randlon) = random_lat_lon(npoints)

        elevations = findEtopoElevation(temp_randlat,temp_randlon)

        (lat,lon) = abovesea(elevations,temp_randlat,temp_randlon)

        randlat = vcat(randlat,lat)
        randlon = vcat(randlon,lon)

        i = length(randlat)
    end

    randlat = randlat[1:npoints];
    randlon = randlon[1:npoints];

## -- Try it

    zoom = 11

    responses = Array{Any}(length(elevations))
    @showprogress 5 for i = 1:length(elevations)
        try
            responses[i] = query_macrostrat(randlat[i], randlon[i], zoom)
        catch
            print("Warning: no data from Macrostrat server \n")
            try
                # Wait and try again
                sleep(10)
                responses[i] = query_macrostrat(randlat[i], randlon[i], zoom)
            catch
                # if still nothing, add warnin
                responses[i] = "No response"
                print("Warning: no response from Macrostrat server\n")
            end
        end
        sleep(0.1)
        if mod(i,10000)==0
            save("data/Responses.jld", "responses", responses, "elevations", elevations, "latitude", randlat, "longitude", randlon, "npoints", npoints)
        end
    end

    save("data/Responses.jld", "responses", responses, "elevations", elevations, "latitude", randlat, "longitude", randlon, "npoints", npoints)

## ---
    # Load from saved version
    retrive_file = load("data/Responses.jld")
    elevations = retrive_file["elevations"]
    randlat = retrive_file["latitude"]
    randlon = retrive_file["longitude"]
    responses = retrive_file["responses"]
    npoints = retrive_file["npoints"]

    # Variables for parsed responses
    rocktype = Array{String}(length(elevations))
    rockdescrip = Array{String}(length(elevations))
    rockname = Array{String}(length(elevations))
    rockstratname = Array{String}(length(elevations))
    rockcomments = Array{String}(length(elevations))

    # Parse saved responses
    for i = 1:length(elevations)
        rocktype[i] = get_macrostrat_lith(responses[i])
        rockdescrip[i] = get_macrostrat_descrip(responses[i])
        rockname[i] = get_macrostrat_name(responses[i])
        rockstratname[i] = get_macrostrat_strat_name(responses[i])
        rockcomments[i] = get_macrostrat_comments(responses[i])
    end

    # Convert to lowercase to match names of rock types
    rocktype = lowercase.(rocktype)
    rockdescrip = lowercase.(rockdescrip)
    rockname = lowercase.(rockname)
    rockstratname = lowercase.(rockstratname)
    rockcomments = lowercase.(rockcomments)

    # Replacing tabs with spaces so as not to be confused with the delimitator
    rocktype = replace.(rocktype, "    ", " ")
    rockdescrip = replace.(rockdescrip, "    ", " ")
    rockname = replace.(rockname, "    ", " ")
    rockstratname = replace.(rockstratname, "    ", " ")
    rockcomments = replace.(rockcomments, "    ", " ")

## ---

    #= Use this link to check the information on certain points in the macrostrat
    url = "https://macrostrat.org/api/mobile/map_query?lat=$(randlat[end])&lng=$(randlon[end])&z=11" =#

## ---

    # Names or partial names of different rock types (from GetBurwellBulkAge.m)
    sedtypes = ["fluv", "clast", "conglomerat", "gravel", "sand", "psamm", "arenit", "arkos", "silt", "mud", "marl", "clay", "shale", "wacke", "argillite", "argillaceous","pelit", "pebble", "mass-wasting", "carbonate", "limestone", "dolo", "chalk", "travertine", "tavertine", "tufa", "evaporite", "salt", "gypsum", "boulder", "gravel", "glaci", "till", "loess", "lluv", "regolith", "debris", "fill", "slide", "unconsolidated", "talus", "stream", "beach", "terrace", "chert", "banded iron", "coal", "anthracite", "peat", "sediment", "laterite", "surficial deposits", "marine deposits", "turbidite", "flysch"];
    igntypes = ["volcanic", "extrusive", "tuff", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", "lahar", "lava", "latite", "basanite", "phonolite", "trachyte", "ignimbrite", "palagonite", "mugearite", "pipe", "plutonic", "intrusive", "granit", "tonalit", "gabbro", "diorit", "monzonit", "syenit", "peridot", "dunit", "harzburg", "dolerit", "diabase", "charnockite", "hypabyssal", "norite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", "porphyry", "megacryst", "bronzitite", "alaskite", "troctolite", "igneous", "silicic ", "mafic", "felsic"];
    mettypes = ["para", "metased", "schist", "quartzite", "marble", "slate", "phyllite", "ortho", "metaign", "serpentin", "amphibolit", "greenstone", "eclogite", "basite", "ultramafitite", "meta", "migma", "gneiss", "granulit", "hornfels", "granofels", "mylonit", "cataclasite", "melange", "gouge", "tecton", "calc silicate", "crystalline"];

    # Check which burwell "lith" rocktypes match one of the rock types
    sed = fill(false,npoints)
    for i = 1:length(sedtypes)
      sed = sed .| contains.(rocktype,sedtypes[i])
    end

    ign = fill(false,npoints)
    for i = 1:length(igntypes)
      ign = ign .| contains.(rocktype,igntypes[i])
    end

    met = fill(false,npoints)
    for i = 1:length(mettypes)
      met = met .| contains.(rocktype,mettypes[i])
    end

    # If we don't find matches in rocktype, try rockname
    not_matched = .~(sed .| ign .| met)
    for i = 1:length(sedtypes)
      sed[not_matched] = sed[not_matched] .| contains.(rockname[not_matched],sedtypes[i])
    end
    for i = 1:length(igntypes)
      ign[not_matched] = ign[not_matched] .| contains.(rockname[not_matched],igntypes[i])
    end
    for i = 1:length(mettypes)
      met[not_matched] = met[not_matched] .| contains.(rockname[not_matched],mettypes[i])
    end

    # Then rockdescrip
    not_matched = .~(sed .| ign .| met)
    for i = 1:length(sedtypes)
      sed[not_matched] = sed[not_matched] .| contains.(rockdescrip[not_matched],sedtypes[i])
    end
    for i = 1:length(igntypes)
      ign[not_matched] = ign[not_matched] .| contains.(rockdescrip[not_matched],igntypes[i])
    end
    for i = 1:length(mettypes)
      met[not_matched] = met[not_matched] .| contains.(rockdescrip[not_matched],mettypes[i])
    end

    not_matched = .~(sed .| ign .| met);
    multi_matched = (sed .& ign) .| (sed .& met) .| (ign .& met);

    number_not_matched = sum(not_matched)
    number_multi_matched = sum(multi_matched)

    print("not matched = $number_not_matched, conflicting matches = $number_multi_matched\n")

## ---

    # Print proportions of sed vs ign, met, and nonsed rocks in responses
    print("sed = ", sum(sed), ", ign + met + ~.sed = ", sum((ign .| met) .& .~sed))

    # Create a file to check matching errors
    writedlm("notmatched.tsv", hcat(rocktype[not_matched], rockname[not_matched], rockdescrip[not_matched], rockstratname[not_matched], rockcomments[not_matched]))
    readdlm("notmatched.tsv")

    writedlm("multimatched.tsv", hcat(rocktype[multi_matched], rockname[multi_matched], rockdescrip[multi_matched], rockstratname[multi_matched], rockcomments[multi_matched]))
    readdlm("multimatched.tsv")

## ---
    # heatmap(etopoelev)
