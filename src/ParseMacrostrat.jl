# Depending on the number of points requested from the Macrostrat / Burwell API, this code
# may take multiple days to run. Recommend running from the command line as a background
# process: nohup julia src/ParseMacrostrat.jl &

## --- Set up
    # Packages
    using StatGeochem
    using LoopVectorization
    using Static
    using Measurements
    using Dates
    using DelimitedFiles
    using JLD
	using HDF5
	using HTTP
	using JSON
    
    # Local utilities
    include("utilities/Utilities.jl")


## --- Generate random points on the continental crust
    npoints = 1_000_000
    savepts = round(Int, npoints / 10_000)
    etopo = get_etopo("elevation")
    rocklat, rocklon, elevations = gen_continental_points(npoints, etopo)

    
## --- Alternatively, restart process from intermediate file
    # retrive_file = load("output/macrostrat/responses250000.jld")
    # responses = retrive_file["responses"]
    # elevations = float.(retrive_file["elevations"])
    # rocklat = float.(retrive_file["latitude"])
    # rocklon = float.(retrive_file["longitude"])    
    # npoints = retrive_file["npoints"]

    # retrive_file = h5open("output/macrostrat/responses250000.h5", "r")

    # close(retrive_file)


## --- Start timer, estimate runtime, and print to terminal
    start = now()
    @info """
    Starting API calls: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))

    Okay, shouldn't take long. Between an hour* and, um, 11 months*. Somewhere in there.

     * Possible minimum run time: $(canonicalize((ceil(Second(npoints/6), Dates.Hour))))
    ** Possible maximum run time: $(canonicalize((ceil(Second(npoints/3), Dates.Hour))))
    """


## --- Get data for each point from the Burwell / Macrostrat API
    savefilename = "responses"
    zoom = 11
    responses = Array{Any}(undef, npoints, 1)
    for i = 1:npoints
        try
            responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
        catch
            try
                # Wait and try again
                sleep(10)
                responses[i] = query_macrostrat(rocklat[i], rocklon[i], zoom)
            catch
                # If still nothing, move on
                responses[i] = "No response"
            end
        end
        sleep(0.05)

        # Checkpoint save and garbage collect every 10,000 points
        if i % 100_000 == 0
            GC.gc()
            parsed = parse_macrostrat_responses(responses, i)
            fid = h5open("output/$savefilename$i.h5", "w")
                fid["rocklat"] = rocklat
                fid["rocklon"] = rocklon
                fid["elevation"] = elevations
                fid["agemax"] = parsed.agemax
                fid["agemin"] = parsed.agemin
                fid["age"] = parsed.age
                fid["rocktype"] = parsed.rocktype
                fid["rockname"] = parsed.rockname
                fid["rockdescrip"] = parsed.rockdescrip
                fid["rockstratname"] = parsed.rockstratname
                fid["rockcomments"] = parsed.rockcomments
                fid["reference"] = parsed.refstrings
                fid["npoints"] = npoints
            close(fid)

            println("Save point $(round(Int, i/100000))/$savepts at $(Dates.Date(now())) $(
                Dates.format(now(), "HH:MM"))"
            )
        end
    end

    # Final save
    # save("output/$savefilename.jld", "responses", responses, "elevations", elevations, 
    #     "latitude", rocklat, "longitude", rocklon, "npoints", npoints
    # )

    # Print success to terminal
    stop = now()
    @info """
    API calls finished at $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).

    Total runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """


## --- Load data from the .jld file
	# @info "Loading Macrostrat file"
    # retrive_file = load("output/pregenerated_responses.jld")
    # @info "Success!"
    
    # responses = retrive_file["responses"]
    # elevations = float.(retrive_file["elevations"])
    # rocklat = float.(retrive_file["latitude"])
    # rocklon = float.(retrive_file["longitude"])    
    # npoints = retrive_file["npoints"]


## --- Parse Macrostrat responses and write data to .h5 file
    @info "Starting to parse data: $(Dates.Date(now())) $(Dates.format(now(), "HH:MM"))"
    parsed = parse_macrostrat_responses(responses, npoints)

    @info "Writing to file: $(Dates.Date(now())) $(Dates.format(now(), "HH:MM"))"
    fid = h5open("output/$savefilename.h5", "w")
        fid["rocklat"] = rocklat
        fid["rocklon"] = rocklon
        fid["elevation"] = elevations
        fid["agemax"] = parsed.agemax
        fid["agemin"] = parsed.agemin
        fid["age"] = parsed.age
        fid["rocktype"] = parsed.rocktype
        fid["rockname"] = parsed.rockname
        fid["rockdescrip"] = parsed.rockdescrip
        fid["rockstratname"] = parsed.rockstratname
        fid["rockcomments"] = parsed.rockcomments
        fid["reference"] = parsed.refstrings
        fid["npoints"] = npoints
    close(fid)

    @info "Data saved successfully! End of process."


## -- End of file
