## --- Set up
    # Packages
    using StatGeochem
    using DelimitedFiles
    using JLD
	using HDF5
	using HTTP
	using JSON
    
    # Local utilities
    include("Utilities.jl")


## --- Generate random points on the continental crust
    npoints = 1_000_000
    etopo = get_etopo("elevation")
    rocklat, rocklon, elevations = gen_continental_points(npoints, etopo)


## --- Get data for each point from the Burwell / Macrostrat API
    zoom = 11
    responses = Array{Any}(undef, npoints, 1)
    @showprogress 5 for i = 1:npoints
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

        # Checkpoint save every 10,000 points
        (i % 10000 == 0) && save("data/$savefilename.jld", "responses", responses, 
            "elevations", elevations, "latitude", rocklat, "longitude", rocklon, 
            "npoints", npoints
        )
    end

    # Final save
    save("data/$savefilename.jld", "responses", responses, "elevations", elevations, 
        "latitude", rocklat, "longitude", rocklon, "npoints", npoints
    )


## --- Load data from the .jld file
	# @info "Loading Macrostrat file"
    # retrive_file = load("output/pregenerated_responses.jld")
    # @info "Success!"
    
    # responses = retrive_file["responses"]
    # elevations = float.(retrive_file["elevations"])
    # rocklat = float.(retrive_file["latitude"])
    # rocklon = float.(retrive_file["longitude"])    
    # npoints = retrive_file["npoints"]


## --- Parse Macrostrat responses from the loaded data
	@info "Starting to parse data"
    # Preallocate
    rocktype = Array{String}(undef, npoints, 1)
    rockdescrip = Array{String}(undef, npoints, 1)
    rockname = Array{String}(undef, npoints, 1)
    rockstratname = Array{String}(undef, npoints, 1)
    rockcomments = Array{String}(undef, npoints, 1)
    agemax = Array{Float64}(undef, npoints, 1)
    agemin = Array{Float64}(undef, npoints, 1)
    age = Array{Float64}(undef, npoints, 1)

    # Parse responses into preallocated arrays
    for i in eachindex(rocktype)
        rocktype[i] = get_macrostrat_lith(responses[i])
        rockdescrip[i] = get_macrostrat_descrip(responses[i])
        rockname[i] = get_macrostrat_name(responses[i])
        rockstratname[i] = get_macrostrat_strat_name(responses[i])
        rockcomments[i] = get_macrostrat_comments(responses[i])
        agemax[i] = get_macrostrat_max_age(responses[i])
        agemin[i] = get_macrostrat_min_age(responses[i])
    end

    # Age of each sample
    @info "Checking age constraints"
    for i in 1:npoints
        age[i] = nanmean([agemax[i], agemin[i]])
    end

    # Filter ages younger than 0 or greater than the age of the earth
    invalid_age = vcat(findall(>(4000), age), findall(<(0), age))
    age[invalid_age] .= NaN

    # Make sure age bounds are in the right order
    for i in 1:npoints
        if agemin[i] > agemax[i]
            tempmax = agemin[i]
            tempmin = agemax[i]
            agemin[i] = tempmax
            agemax[i] = tempmin
        end
    end

    # Convert strings to lowercase so they can be matched to known names of rock types
    rocktype = lowercase.(rocktype)
    rockdescrip = lowercase.(rockdescrip)
    rockname = lowercase.(rockname)
    rockstratname = lowercase.(rockstratname)
    rockcomments = lowercase.(rockcomments)

    # Replace tabs with spaces so they will not be confused with the delimitator when exported
    rocktype = replace.(rocktype, "    " => " ")
    rockdescrip = replace.(rockdescrip, "    " => " ")
    rockname = replace.(rockname, "    " => " ")
    rockstratname = replace.(rockstratname, "    " => " ")
    rockcomments = replace.(rockcomments, "    " => " ")


## --- Parse Macrostat data references
	@info "Parsing references"
    # Preallocate
    authors = Array{String}(undef, npoints, 1)
    years = Array{String}(undef, npoints, 1)
    titles = Array{String}(undef, npoints, 1)
    dois = Array{String}(undef, npoints, 1)

    # Parse into preallocated arrays
    for i = eachindex(authors)
        authors[i] = get_macrostrat_ref_authors(responses[i])
        years[i] =get_macrostrat_ref_year(responses[i])
        titles[i] = get_macrostrat_ref_title(responses[i])
        dois[i] = get_macrostrat_ref_doi(responses[i])
    end
    refstrings = @. authors * " | " * years * " | " * titles * " | " * dois

## --- Write data to a file
	@info "Writing to file"
    header = ["rocklat" "rocklon" "agemax" "agemin" "age" "rocktype" "rockname" "rockdescrip" "rockstratname" "rockcomments" "reference"]
    writedlm("output/pregenerated_responses.tsv", vcat(header, hcat(rocklat, rocklon, agemax, agemin,
        age, rocktype, rockname, rockdescrip, rockstratname, rockcomments, refstrings))
    )


## -- End of file
