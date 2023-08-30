# Functions and scripts for creating and modifying files, generally to be used only
# once, but worth saving.

## --- Save Macrostrat and EarthChem matches to files
    # Note that this script has not been tested in entirety; the file writing was done
    # after the matching had already been run in SampleMatch.jl

    # using HDF5
    # using StatGeochem
    # include("utilities/Utilities.jl")

    # # Macrostrat
    # fid = h5open("$macrostrat_io", "r")
    # macrostrat = (
    #     rocktype = read(fid["rocktype"]),
    #     rockname = read(fid["rockname"]),
    #     rockdescrip = read(fid["rockdescrip"]),
    # )
    # close(fid)
    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
    #     macrostrat.rockdescrip, unmultimatch=false, inclusive=false, source=:macrostrat
    # )
    
    # # EarthChem
    # @info "Loading EarthChem data $(Dates.format(now(), "HH:MM"))"
    # fid = h5open("output/bulk.h5", "r")
    # header = read(fid["bulk"]["header"])
    # data = read(fid["bulk"]["data"])
    # bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])


    # path = fid["bulktext"]["sampledata"]
    # header = read(path["header"])
    # index = read(path["index"])
    # target = ["Rock_Name", "Type", "Material"]
    # targetind = [findall(==(i), header)[1] for i in target]
    # bulktext = NamedTuple{Tuple(Symbol.(target))}(
    #     [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] 
    #         for i in eachindex(target)
    #     ]
    # )
    # close(fid)

    # bulk_cats = match_rocktype(bulktext.Rock_Name, bulktext.Type, bulktext.Material; 
    #     unmultimatch=false, inclusive=false, source=:earthchem
    # )

    # Get rock names for each Macrostrat sample
    # name_cats = match_rockname(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    # rocknames = string.(keys(name_cats))

    # # Get EarthChem samples for each rock name
    # typelist = get_rock_class(false, true)      # Subtypes, major types include minors
    # nbulk = length(bulktext.Rock_Name)
    # bulk_lookup = NamedTuple{keys(name_cats)}([falses(nbulk) for _ in eachindex(name_cats)])

    # p = Progress(length(rocknames)+1, desc="Finding EarthChem samples for each rock name")
    # next!(p)
    # for i in eachindex(rocknames)
    #     bulk_lookup[i] .= find_earthchem(rocknames[i], bulktext.Rock_Name, bulktext.Type, 
    #         bulktext.Material
    #     )

    #     # If no matches, jump up a class
    #     if count(bulk_lookup[i]) == 0
    #         newsearch = class_up(typelist, rocknames[i])
    #         newsearch==:carb && (newsearch=="carbonate")    # No carbonatites!
    #         bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
    #             bulktext.Type, bulktext.Material
    #         )

    #         # If still no matches, jump up a class again
    #         if count(bulk_lookup[i]) == 0
    #             newsearch = class_up(typelist, string(newsearch))
    #             bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
    #                 bulktext.Type, bulktext.Material
    #             )
    #         end
    #     end
    #     next!(p)
    # end

    # fid = h5open("output/matches.h5", "w")
    # g = create_group(fid, "vars")

    # # Macrostrat rock types
    # a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    # for i in eachindex(keys(macro_cats))
    #     for j in eachindex(macro_cats[i])
    #         a[j,i] = ifelse(macro_cats[i][j], 1, 0)
    #     end
    # end
    # g["macro_cats"] = a
    # g["macro_cats_head"] = string.(collect(keys(macro_cats)))

    # # EarthChem rock types
    # a = Array{Int64}(undef, length(bulk_cats[1]), length(bulk_cats))
    # for i in eachindex(keys(bulk_cats))
    #     for j in eachindex(bulk_cats[i])
    #         a[j,i] = ifelse(bulk_cats[i][j], 1, 0)
    #     end
    # end
    # g["bulk_cats"] = a
    # g["bulk_cats_head"] = string.(collect(keys(bulk_cats))) 

    # # Macrostrat rock names
    # a = Array{Int64}(undef, length(name_cats[1]), length(name_cats))
    # for i in eachindex(keys(name_cats))
    #     for j in eachindex(name_cats[i])
    #         a[j,i] = ifelse(name_cats[i][j], 1, 0)
    #     end
    # end
    # g["name_cats"] = a
    # g["name_cats_head"] = string.(collect(keys(name_cats))) 

    # # EarthChem rock names
    # a = Array{Int64}(undef, length(bulk_lookup[1]), length(bulk_lookup))
    # for i in eachindex(keys(bulk_lookup))
    #     for j in eachindex(bulk_lookup[i])
    #         a[j,i] = ifelse(bulk_lookup[i][j], 1, 0)
    #     end
    # end
    # g["bulk_lookup"] = a
    # g["bulk_lookup_head"] = string.(collect(keys(bulk_lookup)))

    # close(fid)

    
## --- Add matches to Macrostrat file
    using HDF5
    using StatGeochem
    using Static
    using LoopVectorization
    using Measurements
    include("utilities/Utilities.jl")

    # Macrostrat
    mfid = h5open("$macrostrat_io", "r")
    # macrostrat = (
    #     rocktype = read(mfid["rocktype"]),
    #     rockname = read(mfid["rockname"]),
    #     rockdescrip = read(mfid["rockdescrip"]),
    # )
    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
    #     macrostrat.rockdescrip, unmultimatch=false, inclusive=false, source=:macrostrat
    # )
    # 
    # name_cats = match_rockname(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    # rocknames = string.(keys(name_cats))
    
    # Types
    fid = h5open("output/matches.h5", "r")

    # New Macrostrat
    newfid = h5open("output/250K_responses_new.h5", "w")
    g = create_group(newfid, "vars")
        copy_object(mfid["age"], g, "age")
        copy_object(mfid["agemax"], g, "agemax")
        copy_object(mfid["agemin"], g, "agemin")
        copy_object(mfid["elevation"], g, "elevation")
        copy_object(mfid["npoints"], g, "npoints")
        copy_object(mfid["reference"], g, "reference")
        copy_object(mfid["rockcomments"], g, "rockcomments")
        copy_object(mfid["rockdescrip"], g, "rockdescrip")
        copy_object(mfid["rocklat"], g, "rocklat")
        copy_object(mfid["rocklon"], g, "rocklon")
        copy_object(mfid["rockname"], g, "rockname")
        copy_object(mfid["rockstratname"], g, "rockstratname")
        copy_object(mfid["rocktype"], g, "rocktype")
    g = create_group(newfid, "type")
        copy_object(fid["vars"]["macro_cats"], g, "macro_cats")
        copy_object(fid["vars"]["macro_cats_head"], g, "macro_cats_head")
        copy_object(fid["vars"]["name_cats"], g, "name_cats")
        copy_object(fid["vars"]["name_cats_head"], g, "name_cats_head")
    close(fid)
    close(newfid)
    close(mfid)
    

## --- Add matches to EarthChem file
    using HDF5
    using StatGeochem
    include("utilities/Utilities.jl")

    bfid = h5open("output/bulk.h5", "r")

    # header = read(bfid["bulk"]["header"])
    # data = read(bfid["bulk"]["data"])
    # bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    # 
    # path = bfid["bulktext"]["sampledata"]
    # header = read(path["header"])
    # index = read(path["index"])
    # target = ["Rock_Name", "Type", "Material"]
    # targetind = [findall(==(i), header)[1] for i in target]
    # bulktext = NamedTuple{Tuple(Symbol.(target))}(
    #     [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] 
    #         for i in eachindex(target)
    #     ]
    # )
    # 
    # bulk_cats = match_rocktype(bulktext.Rock_Name, bulktext.Type, bulktext.Material; 
    #     unmultimatch=false, inclusive=false, source=:earthchem
    # )
    # 
    # typelist = get_rock_class(false, true)      # Subtypes, major types include minors
    # nbulk = length(bulktext.Rock_Name)
    # bulk_lookup = NamedTuple{keys(name_cats)}([falses(nbulk) for _ in eachindex(name_cats)])
    # 
    # p = Progress(length(rocknames)+1, desc="Finding EarthChem samples for each rock name")
    # next!(p)
    # for i in eachindex(rocknames)
    #     bulk_lookup[i] .= find_earthchem(rocknames[i], bulktext.Rock_Name, bulktext.Type, 
    #         bulktext.Material
    #     )
    # 
    #     # If no matches, jump up a class
    #     if count(bulk_lookup[i]) == 0
    #         newsearch = class_up(typelist, rocknames[i])
    #         newsearch==:carb && (newsearch=="carbonate")    # No carbonatites!
    #         bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
    #             bulktext.Type, bulktext.Material
    #         )
    # 
    #         # If still no matches, jump up a class again
    #         if count(bulk_lookup[i]) == 0
    #             newsearch = class_up(typelist, string(newsearch))
    #             bulk_lookup[i] .= find_earthchem(string(newsearch), bulktext.Rock_Name, 
    #                 bulktext.Type, bulktext.Material
    #             )
    #         end
    #     end
    #     next!(p)
    # end

    # Types
    fid = h5open("output/matches.h5", "r")

    # New EarthChem
    newfid = h5open("output/bulk_new.h5", "w")
    copy_object(bfid["bulktext"], newfid, "bulktext")
    g = create_group(newfid, "bulk")
        copy_object(bfid["bulk"]["data"], g, "data")
        copy_object(bfid["bulk"]["header"], g, "header")
    g = create_group(newfid, "bulktypes")
        copy_object(fid["vars"]["bulk_cats"], g, "bulk_cats")
        copy_object(fid["vars"]["bulk_cats_head"], g, "bulk_cats_head")
        copy_object(fid["vars"]["bulk_lookup"], g, "bulk_lookup")
        copy_object(fid["vars"]["bulk_lookup_head"], g, "bulk_lookup_head")
    close(fid)
    close(newfid)
    close(bfid)
    

## --- Calculate inverse spatial weight for each rock name
    # # This current method means samples without a latitude or longitude will never get picked.
    # # In theory, I could undo that by assigning a weight to samples without spatial data;
    # # probably the mean of the weights (although this could cause some weird statistical
    # # behavior depending on the distribution...). In any case, it may not end up mattering
    # # if there's enough samples to choose from.
    # using HDF5
    # using StatGeochem

    # # EarthChem
    # fid = h5open("output/bulk.h5", "r")
    #     header = read(fid["bulk"]["header"])
    #     data = read(fid["bulk"]["data"])
    # close(fid)
    # bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    # nbulk = length(bulk[1])

    # # Matched rock names--no sense calculating spatial weights for rocks that aren't matched
    # fid = h5open("output/matches.h5", "r")
    #     rocknames = read(fid["vars"]["bulk_lookup_head"])
    #     data = read(fid["vars"]["bulk_lookup"])
    #     data = @. data > 0
    # close(fid)
    # bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([data[:,i] for i in eachindex(rocknames)])

    # # Preallocate
    # @info "Calculating inverse spatial weights for each rock name"
    # spatial_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([fill(NaN, nbulk) for _ in eachindex(rocknames)])
    # for n in eachindex(keys(spatial_lookup))
    #     println("$n ($n/$(length(rocknames))) \n")
    #     spatial_lookup[n][bulk_lookup[n]] .= invweight_location(bulk.Latitude[bulk_lookup[n]], 
    #         bulk.Latitude[bulk_lookup[n]]
    #     )
    # end

    # # Save to file
    # A = Array{Float64}(undef, nbulk, length(spatial_lookup))
    # for i in eachindex(keys(spatial_lookup))
    #     A[:,i] = spatial_lookup[i]
    # end

    # fid = h5open("output/invspatial.h5", "w")
    #     fid["header"] = rocknames
    #     fid["k"] = A
    # close(fid)


## --- Save basin lat / lon coordinates to a file
    # Get stuff
    run(`gunzip data/octopus/crn_basins_global.kml.gz`);
    (str, isfirstcoord, nbasins, subbasins) = load_octopus_kml("data/octopus/crn_basins_global.kml");
    run(`gzip data/octopus/crn_basins_global.kml`);  

    (basin_polygon_n, basin_polygon_lat, basin_polygon_lon) = parse_octopus_polygon_outlines(str,isfirstcoord)

    # File time!
    fid = h5open("output/basin_coordinates.h5", "w")
    g = create_group(fid, "vars")
    for i in eachindex(basin_polygon_lat, basin_polygon_lon)
        g[lpad(i, 4, "0")] = hcat(basin_polygon_lat[i], basin_polygon_lon[i])
    end
    close(fid)

## --- End of File