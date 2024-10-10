# Depending on the number of points requested from the Macrostrat / Burwell API, this code
# may take multiple days to run. Recommend running from the command line as a background
# process: nohup julia --project="Project.toml" src/ParseMacrostrat.jl &

## --- Set up
    # Packages
    using RockWeatheringFlux
    using Dates
    using DelimitedFiles, HDF5, JLD, HTTP, JSON


## --- Generate random points on the continental crust
    @info "Script started."
    npoints = RockWeatheringFlux.N
    saveinterval = 100_000
    savepts = round(Int, npoints / saveinterval)
    # FIXME: files
    etopo = h5read(etopo_home, "vars/elevation")
    rocklat, rocklon, elevations = gen_continental_points(npoints, etopo)

    savefilename = "responses"


## --- Get data for each point from the Burwell / Macrostrat API
    start = now()
    @info """
    Starting API calls: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))

    Okay, shouldn't take long. Between an hour and, um, 11 months*. Somewhere in there.

     * A set of 100,000 samples takes 4-5 hours.
    """

    responses = Array{Union{Dict{String, Any}, String}}(undef, npoints, 1)
    skipcount = 0
    for i = 1:npoints
        try
            responses[i] = query_macrostrat(rocklat[i], rocklon[i])
        catch
            try
                # Wait and try again
                sleep(1)
                responses[i] = query_macrostrat(rocklat[i], rocklon[i])
            catch
                # If still nothing, move on
                responses[i] = "No response"
                global skipcount += 1
            end
        end
        sleep(0.05)

        # Checkpoint save and garbage collect at save point intervals
        if i % saveinterval == 0
            GC.gc()
            parsed = parse_macrostrat_responses(responses, i)
            fid = h5open("$macrostrat_parsed$i.h5", "w")
            g = create_group(fid, "vars")
                g["rocklat"] = rocklat
                g["rocklon"] = rocklon
                g["elevation"] = elevations
                g["agemax"] = parsed.agemax
                g["agemin"] = parsed.agemin
                g["age"] = parsed.age
                g["rocktype"] = parsed.rocktype
                g["rockname"] = parsed.rockname
                g["rockdescrip"] = parsed.rockdescrip
                g["rockstratname"] = parsed.rockstratname
                g["rockcomments"] = parsed.rockcomments
                g["reference"] = parsed.refstrings
                g["npoints"] = npoints
            close(fid)

            # Print info to terminal
            point = round(Int, i/saveinterval)
            day = Dates.Date(now())
            time = Dates.format(now(), "HH:MM")
            println("Save point $point/$savepts at $day $time. No response from $skipcount / $i samples.")
        end
    end

    # Save unparsed output, just in case
    save(macrostrat_raw, "responses", responses, "elevations", elevations, 
        "latitude", rocklat, "longitude", rocklon, "npoints", npoints
    )

    # Print success to terminal
    stop = now()
    @info """
    API calls finished at $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).

    Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """


## --- Parse Macrostrat responses and write data to a file
    @info "Parsing data: $(Dates.Date(now())) $(Dates.format(now(), "HH:MM"))"
    parsed = parse_macrostrat_responses(responses, npoints)

    # Write data to file
    @info "Writing to file: $macrostrat_parsed$npoints.h5 at $(Dates.Date(now())) $(Dates.format(now(), "HH:MM"))"
    fid = h5open("$macrostrat_parsed$npoints.h5", "w")
    g = create_group(fid, "vars")
        g["rocklat"] = rocklat
        g["rocklon"] = rocklon
        g["elevation"] = elevations
        g["agemax"] = parsed.agemax
        g["agemin"] = parsed.agemin
        g["age"] = parsed.age
        g["rocktype"] = parsed.rocktype
        g["rockname"] = parsed.rockname
        g["rockdescrip"] = parsed.rockdescrip
        g["rockstratname"] = parsed.rockstratname
        g["rockcomments"] = parsed.rockcomments
        g["reference"] = parsed.refstrings
        g["npoints"] = npoints

    bulktypes = create_group(fid, "type")
        # Rock types
        macro_cats = match_rocktype(parsed.rocktype, parsed.rockname, parsed.rockdescrip, showprogress=false)
        a = Array{Int64}(undef, length(macro_cats[1]), length(macro_cats))
        for i in eachindex(keys(macro_cats))
            for j in eachindex(macro_cats[i])
                a[j,i] = ifelse(macro_cats[i][j], 1, 0)
            end
        end
        bulktypes["macro_cats"] = a
        bulktypes["macro_cats_head"] = string.(collect(keys(macro_cats))) 
        
        # Metamorphic rocks 
        metamorph_cats = find_metamorphics(parsed.rocktype, parsed.rockname, parsed.rockdescrip, showprogress=false)
        a = Array{Int64}(undef, length(metamorph_cats[1]), length(metamorph_cats))
        for i in eachindex(keys(metamorph_cats))
            for j in eachindex(metamorph_cats[i])
                a[j,i] = ifelse(metamorph_cats[i][j], 1, 0)
            end
        end
        bulktypes["metamorphic_cats"] = a
    
    close(fid)

    @info "Data saved successfully! End of process."


## --- Just see if everything makes sense 
    # using Plots 
    # fid = h5open(macrostrat_io, "r")
    # rocklat = read(fid["vars"]["rocklat"])[t]
    # rocklon = read(fid["vars"]["rocklon"])[t]
    # close(fid)

    # mapplot(rocklon, rocklat, 
    #     label="",
    #     markersize=1, 
    #     msc=:auto
    # )


## -- End of file
