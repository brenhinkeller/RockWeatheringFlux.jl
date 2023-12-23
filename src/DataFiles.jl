# Create a duplicate Macrostrat / Burwell file with new rock type and rock name matches

## --- Set up
    using RockWeatheringFlux
    using HDF5


## --- Match rock names
    @info "Re-matching samples from $macrostrat_io"

    # Macrostrat
    oldfid = h5open(macrostrat_io, "r")

    # Get matches
    macrostrat = (
        rocktype = read(oldfid["vars"]["rocktype"]),
        rockname = read(oldfid["vars"]["rockname"]),
        rockdescrip = read(oldfid["vars"]["rockdescrip"]),
    )
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)


## --- Put data in new file
    # Create a new file and put the old file in it
    fid, fex = split("$macrostrat_io", ".")
    newfid = h5open(fid * "_2.h5" , "w")
    copy_object(oldfid, "vars", newfid, "vars")    # The stuff that isn't changed

    # Rock types and rock names
    bulktypes = create_group(newfid, "type")
    a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    for i in eachindex(keys(macro_cats))
        for j in eachindex(macro_cats[i])
            a[j,i] = ifelse(macro_cats[i][j], 1, 0)
        end
    end
    bulktypes["macro_cats"] = a
    bulktypes["macro_cats_head"] = string.(collect(keys(macro_cats))) 

    close(newfid)
    close(oldfid)


## --- Or if you're working with really old files that don't have the updated structure:
    # newfid = h5open(macrostrat_io, "w")

    # # Copy over existing data
    # g = create_group(newfid, "vars")
    #     copy_object(oldfid["age"], g, "age")
    #     copy_object(oldfid["agemax"], g, "agemax")
    #     copy_object(oldfid["agemin"], g, "agemin")
    #     copy_object(oldfid["elevation"], g, "elevation")
    #     copy_object(oldfid["npoints"], g, "npoints")
    #     copy_object(oldfid["reference"], g, "reference")
    #     copy_object(oldfid["rockcomments"], g, "rockcomments")
    #     copy_object(oldfid["rockdescrip"], g, "rockdescrip")
    #     copy_object(oldfid["rocklat"], g, "rocklat")
    #     copy_object(oldfid["rocklon"], g, "rocklon")
    #     copy_object(oldfid["rockname"], g, "rockname")
    #     copy_object(oldfid["rockstratname"], g, "rockstratname")
    #     copy_object(oldfid["rocktype"], g, "rocktype")

    # # Rock types and rock names
    # bulktypes = create_group(newfid, "type")

    # # Rock types
    # a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    # for i in eachindex(keys(macro_cats))
    #     for j in eachindex(macro_cats[i])
    #         a[j,i] = ifelse(macro_cats[i][j], 1, 0)
    #     end
    # end
    # bulktypes["macro_cats"] = a
    # bulktypes["macro_cats_head"] = string.(collect(keys(macro_cats))) 

    # close(newfid)


## --- End of File