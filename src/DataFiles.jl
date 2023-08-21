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

    
## --- Calculate inverse spatial weight for each rock name
    # This current method means samples without a latitude or longitude will never get picked.
    # In theory, I could undo that by assigning a weight to samples without spatial data;
    # probably the mean of the weights (although this could cause some weird statistical
    # behavior depending on the distribution...). In any case, it may not end up mattering
    # if there's enough samples to choose from.
    using HDF5
    using StatGeochem

    # EarthChem
    fid = h5open("output/bulk.h5", "r")
        header = read(fid["bulk"]["header"])
        data = read(fid["bulk"]["data"])
    close(fid)
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    nbulk = length(bulk[1])

    # Matched rock names , really only needed for the names and not the matches
    fid = h5open("output/matches.h5", "r")
        rocknames = read(fid["vars"]["name_cats_head"])
    close(fid)

    # Preallocate
    spatial_lookup = NamedTuple{Tuple(rocknames)}([fill(NaN, nbulk) for _ in eachindex(rocknames)])

    @info "Calculating inverse spatial weights for each rock name"
    for n in eachindex(keys(spatial_lookup))
        println("$n ($n/$(length(rocknames))) \n")
        spatial_lookup[n][bulk_lookup[n]] .= invweight_location(bulk.Latitude[bulk_lookup[n]], 
            bulk.Latitude[bulk_lookup[n]]
        )
    end

    # Save to file
    A = Array{Float64}(undef, nbulk, length(spatial_lookup))
    for i in eachindex(keys(spatial_lookup))
        A[:,i] = spatial_lookup[i]
    end

    fid = h5open("output/invspatial.h5", "w")
        fid["header"] = rocknames
        fid["k"] = A
    close(fid)


## --- End of File