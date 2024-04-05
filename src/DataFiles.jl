# Create a duplicate Macrostrat / Burwell file with new rock type and rock name matches

## --- Set up
    using RockWeatheringFlux
    using HDF5

    # Load file 
    @info "Re-matching samples from $macrostrat_io"
    oldfid = h5open(macrostrat_io, "r")

    # Open a new file
    fid = split(macrostrat_io, ".")[1]
    newfid = h5open(fid * "_2.h5" , "w")


## --- VERSION 1: FILE MATCHES EXISTING FORMAT
    # Match rock class
    macrostrat = (
        rocktype = read(oldfid["vars"]["rocktype"]),
        rockname = read(oldfid["vars"]["rockname"]),
        rockdescrip = read(oldfid["vars"]["rockdescrip"]),
    )
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    metamorph_cats = find_metamorphics(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    # Copy old data over to the new file 
    copy_object(oldfid, "vars", newfid, "vars")    # The stuff that isn't changed


## --- VERSION 2: FILE DOES NOT MATCH EXISTING FORMAT
    # # Match rock class
    # macrostrat = (
    #     rocktype = read(oldfid["rocktype"]),
    #     rockname = read(oldfid["rockname"]),
    #     rockdescrip = read(oldfid["rockdescrip"]),
    # )
    # macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)
    # metamorph_cats = find_metamorphics(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    # # Create vars group and copy over static data
    # g = create_group(newfid, "vars")
    #     copy_object(oldfid["age"], g, "age")
    #     copy_object(oldfid["agemax"], g, "agemax")
    #     copy_object(oldfid["agemin"], g, "agemin")
    #     copy_object(oldfid["elevation"], g, "elevation")
    #     copy_object(oldfid["npoints"], g, "npoints")
    #     copy_object(oldfid["reference"], g, "reference")
    #     copy_object(oldfid["rockcomments"], g, "rockcomments")
    #     copy_object(oldfid["rockdescrip"], g, "rockdescrip")
    #     copy_object(oldfid["rockname"], g, "rockname")
    #     copy_object(oldfid["rockstratname"], g, "rockstratname")
    #     copy_object(oldfid["rocktype"], g, "rocktype")

    # # Trim rocklat and rocklon, if necessary 
    # npoints = length(macrostrat.rocktype)
    # g["rocklat"] = read(oldfid["rocklat"])[1:npoints]
    # g["rocklon"] = read(oldfid["rocklon"])[1:npoints]
    

## --- BOTH VERSIONS: Copy new rock classes over to the new file 
    bulktypes = create_group(newfid, "type")

    a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
    for i in eachindex(keys(macro_cats))
        for j in eachindex(macro_cats[i])
            a[j,i] = ifelse(macro_cats[i][j], 1, 0)
        end
    end

    b = similar(a)
    for i in eachindex(keys(metamorph_cats))
        for j in eachindex(metamorph_cats[i])
            b[j,i] = ifelse(metamorph_cats[i][j], 1, 0)
        end
    end

    bulktypes["macro_cats"] = a
    bulktypes["metamorphic_cats"] = b
    bulktypes["macro_cats_head"] = string.(collect(keys(macro_cats))) 

    close(newfid)
    close(oldfid)


## --- Rename the files so I don't have to 
    fid = split("$macrostrat_io", ".")[1]
    run(`mv $fid.h5 $(fid)_old.h5`)
    run(`mv $(fid)_2.h5 $fid.h5`)

    
## --- End of File