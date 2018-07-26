
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
    # using HTTP, JSON
    include("Utilities.jl")

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
        elevations = find_etopoelev(randlat,randlon)
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
        sleep(0.1)
        if mod(i,10000)==0
            save("data/responses.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)
        end
    end

    save("data/responses.jld", "responses", responses, "elevations", elevations, "latitude", rocklat, "longitude", rocklon, "npoints", npoints)

## --- Load from saved version (if applicable)

    retrive_file = load("data/responses.jld")
    elevations = retrive_file["elevations"]
    rocklat = retrive_file["latitude"]
    rocklon = retrive_file["longitude"]
    responses = retrive_file["responses"]
    npoints = retrive_file["npoints"]

## --- Parse the macrostrat repponses

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
    sedtypes = ["fluv", " clast", "siliciclast", "conglomerat", "gravel", "sand", "psamm", "arenit", "arkos", "silt", "mud", "marl", "clay", "shale", "wacke", "argillite", "argillaceous","pelit", "pebble", "mass-wasting", "carbonate", "limestone", "dolo", "caliche", "chalk", "travertine", "tavertine", "teravertine", "tufa", "evaporite", " salt", "salt flat", "gypsum", "boulder", "gravel", "glaci", "till", "loess", "lluv", "regolith", "soil", "debris", "slide", "unconsolidated", "talus", "stream", "beach", "terrace", "chert", "banded iron", "coal", "anthracite", "peat", "swamp", "marsh", "sediment", "laterite", "surficial deposits", "marine deposits", "turbidite", "flysch"];
    igntypes = ["volcanic", "extrusive", "tuff", "basalt", "andesit", "dacit", "rhyolit", "pillow", "carbonatite", "tephra", "obsidian", "ash", "scoria", "pumice", "cinder", "lahar", "lava", "latite", "basanite", "phonolite", "fonolito", "trachyte", "ignimbrite", "palagonite", "mugearite", "pipe", "plutonic", "intrusive", "granit", "tonalit", "gabbro", "diorit", "monzonit", "syenit", "peridot", "dunit", "harzburg", "dolerit", "diabase", "charnockite", "hypabyssal", "norite", "pegmatite", "aplite", "trond", "essexite", "pyroxenite", "adamellite", "porphyry", "megacryst", "bronzitite", "alaskite", "troctolite", "igneous", "silicic ", "mafic", "felsic"];
    mettypes = ["para", "metased", "schist", "quartzite", "marble", "slate", "phyllite", "ortho", "metaign", "serpentin", "amphibolit", "greenstone", "eclogite", "basite", "ultramafitite", "meta", "migma", "gneiss", "granulit", "hornfels", "granofels", "mylonit", "cataclasite", "melange", "gouge", "tecton", "calc silicate", "crystalline basement"];

    # Allocate arrays for each sample for each rock type
    sed = fill(false,npoints)
    ign = fill(false,npoints)
    met = fill(false,npoints)

    # Check which burwell "lith" rocktypes match one of the rock types
    # Try the "major:: {...}" type first, if present
    for i = 1:length(sedtypes)
      sed = sed .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? contains(x[1], sedtypes[i]) : false )
    end
    for i = 1:length(igntypes)
      ign = ign .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? contains(x[1], igntypes[i]) : false )
    end
    for i = 1:length(mettypes)
      met = met .|  ( match.(r"major.*?{(.*?)}", rocktype) .|> x -> isa(x,RegexMatch) ? contains(x[1], mettypes[i]) : false )
    end

    # Then check the rest of rocktype
    not_matched = .~(sed .| ign .| met)
    for i = 1:length(sedtypes)
      sed[not_matched] = sed[not_matched] .| contains.(rocktype[not_matched],sedtypes[i])
    end
    for i = 1:length(igntypes)
      ign[not_matched] = ign[not_matched] .| contains.(rocktype[not_matched],igntypes[i])
    end
    for i = 1:length(mettypes)
      met[not_matched] = met[not_matched] .| contains.(rocktype[not_matched],mettypes[i])
    end

    # Then rockname
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
    print("sed = ", sum(sed), ", ign + (met & ~sed)= ", sum(ign .| (met .& .~sed)))

    # Create a file to check matching errors
    writedlm("notmatched.tsv", hcat(rocktype[not_matched], rockname[not_matched], rockdescrip[not_matched], rockstratname[not_matched], rockcomments[not_matched]))

    writedlm("multimatched.tsv", hcat(rocktype[multi_matched], rockname[multi_matched], rockdescrip[multi_matched], rockstratname[multi_matched], rockcomments[multi_matched]))

    writedlm("ignsed.tsv", hcat(rocktype[ign .& sed], rockname[ign .& sed], rockdescrip[ign .& sed], rockstratname[ign .& sed], rockcomments[ign .& sed]))

## --- Find slope and erosion rate

    function Emmkyr(slp)
        return 10^(slp*0.00567517 + 0.971075)
    end

    slope = get_srtm15plus_aveslope("slope")
    rockslope = find_srtm15plus_aveslope_around(slope,rocklat,rocklon, halfwidth_lat=0, halfwidth_lon=0)
    # rockslope = find_srtm15plus_aveslope(slope,rocklat,rocklon)
    # rockslope[rockslope.>400] = NaN;
    rockEmmkyr = Emmkyr.(rockslope)

    sedEsum = nansum(rockEmmkyr[sed])
    crystEsum = nansum(rockEmmkyr[ign .| (met .& .~sed)])

    sedEmean = nanmean(rockEmmkyr[sed])
    crystEmean = nanmean(rockEmmkyr[ign .| (met .& .~sed)])

    propotionSedArea = sum(sed)/(sum(sed)+sum(ign .| (met .& .~sed)));
    print("$(round(100*propotionSedArea,3))% sed area, $(round(100*(1-propotionSedArea),3))% crystalline area\n")
    print("sed E sum: $sedEsum cryst E sum: $crystEsum\n")
    print("$(round(100*sedEsum/(sedEsum+crystEsum),3))% sed flux, $(round(100*crystEsum/(sedEsum+crystEsum),3))% crystalline flux\n")
    print("sed mean E: $sedEmean crystalline mean E: $crystEmean\n")

## --- For comparison

    # Phanerozoic
    sedimentSubduction_km3_yr = 1.64 # Clift, 2009
    continentSedimentation_km3_yr = 0.8616 # Macrostrat
    continentalErosion_km3_yr = sedimentSubduction_km3_yr + sedimentSubduction_km3_yr

    continentalErosionQuaternary_km3_yr = 6.5 #\cite{Clift:2009fm}  / Milliman (1990)
    continentalArea_km2 =  1.4894E8 # km^2

    continentalErosion_km_yr = continentalErosion_km3_yr / continentalArea_km2
    continentalErosion_mm_kyr = continentalErosion_km_yr * 1000 * 1000 * 1000
