# Functions and scripts for creating and modifying files, generally to be used only
# once, but worth saving.
    
## --- Add matches to Macrostrat file
    using HDF5
    using StatGeochem
    using Static
    using LoopVectorization
    using Measurements
    using ProgressMeter
    
    include("utilities/Utilities.jl")

    # Macrostrat
    # mfid = h5open("$macrostrat_io", "r")
    mfid = h5open("output/250K_responses.h5", "r")

    # Get matches
    macrostrat = (
        rocktype = read(mfid["vars"]["rocktype"]),
        rockname = read(mfid["vars"]["rockname"]),
        rockdescrip = read(mfid["vars"]["rockdescrip"]),
    )
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, 
        macrostrat.rockdescrip, unmultimatch=false, inclusive=false, source=:macrostrat
    )
    
    name_cats = match_rockname(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    
    # New Macrostrat file
    newfid = h5open("output/temp_new_macrostrat.h5", "w")

    # Copy over existing data
    g = create_group(newfid, "vars")
        copy_object(mfid["vars"]["age"], g, "age")
        copy_object(mfid["vars"]["agemax"], g, "agemax")
        copy_object(mfid["vars"]["agemin"], g, "agemin")
        copy_object(mfid["vars"]["elevation"], g, "elevation")
        copy_object(mfid["vars"]["npoints"], g, "npoints")
        copy_object(mfid["vars"]["reference"], g, "reference")
        copy_object(mfid["vars"]["rockcomments"], g, "rockcomments")
        copy_object(mfid["vars"]["rockdescrip"], g, "rockdescrip")
        copy_object(mfid["vars"]["rocklat"], g, "rocklat")
        copy_object(mfid["vars"]["rocklon"], g, "rocklon")
        copy_object(mfid["vars"]["rockname"], g, "rockname")
        copy_object(mfid["vars"]["rockstratname"], g, "rockstratname")
        copy_object(mfid["vars"]["rocktype"], g, "rocktype")

    # Rock types and rock names
    bulktypes = create_group(newfid, "type")

    # Rock types
    a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    for i in eachindex(keys(macro_cats))
        for j in eachindex(macro_cats[i])
            a[j,i] = ifelse(macro_cats[i][j], 1, 0)
        end
    end
    bulktypes["macro_cats"] = a
    bulktypes["macro_cats_head"] = string.(collect(keys(macro_cats))) 

    # Rock names
    a = Array{Int64}(undef, length(name_cats[1]), length(name_cats))
    for i in eachindex(keys(name_cats))
        for j in eachindex(name_cats[i])
            a[j,i] = ifelse(name_cats[i][j], 1, 0)
        end
    end
    bulktypes["name_cats"] = a
    bulktypes["name_cats_head"] = string.(collect(keys(name_cats)))

    close(newfid)
    close(mfid)


## --- End of File