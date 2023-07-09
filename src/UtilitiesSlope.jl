## --- Unzip the kml file
    """
    ```julia
    str, isfirstcoord, nbasins, subbasins = load_octopus_kml(filename)
    ```
    Load OCTOPUS .kml file.

    Return unparsed variables `str`, and the number of basins `nbasins`. Each basin may 
    contain one or more `subbasins`. 
    
    The `isfirstcoord` `BitVector` is used internally to get whole basin data (combining
    subbasins into their parent basins). It cannot be used to index into a `Vector`.
    """
    function load_octopus_kml(filename)
        # Check number and ordering of coordinate lists
        fid = open("data/octopus/crn_basins_global.kml")

        # Count the number of lines in the file
        i = 0
        for line in eachline(fid)
            i += 1
        end
        isfirstcoord = Array{Bool}(undef, i)

        # Parse the coordinates (basin polygon outlines) in each placemark
        i = 0
        ncoords = 0
        lastcoord = 0
        lastplacemark = 0
        seekstart(fid)
        for line in eachline(fid)
            i += 1
            if occursin("</ExtendedData>", line)
                lastplacemark = i
            elseif occursin("<coordinates>", line)
                ncoords += 1
                # If we've had a new placemark more recently than a coordinate, then
                # this is the first coord of the placemark
                if lastplacemark > lastcoord
                    isfirstcoord[ncoords] = true
                else
                    isfirstcoord[ncoords] = false
                end
                lastcoord = i
            end
        end
        isfirstcoord = isfirstcoord[1:ncoords]

        # Read entire string
        seekstart(fid)
        str = read(fid,String)

        # Close the file
        close(fid)

        nbasins = 0
        subbasins = Array{Any}(undef, count(isfirstcoord))
        for i = 1:length(isfirstcoord)
            if isfirstcoord[i]
                nbasins += 1
                subbasins[nbasins] = i
            else
                subbasins[nbasins] = vcat(subbasins[nbasins],i)
            end
        end

        return str, isfirstcoord, nbasins, subbasins
    end


## --- Parse the file
    """
    ```julia
    parse_octopus_kml_variables(str)

    Parse .kml formatted data `str` into a `Dict`.
    ```
    """
    function parse_octopus_kml_variables(str)
        # Find a list of all the string variables
        lm = eachmatch(r"<SimpleField type=\"string\" name=\"(.*?)\"", str);
        stringvars = unique(map(x -> x[1], lm))

        # Find a list of all the numeric variables
        lm = eachmatch(r"<SimpleField type=\"double\" name=\"(.*?)\"", str);
        numvars = unique(map(x -> x[1], lm))

        # Make a dict to store all our parsed data
        data = Dict();

        # Parse all the string variables
        for i = 1:length(stringvars)
            lm = eachmatch(Regex("<SimpleData name=\"$(stringvars[i])\">(.*?)</SimpleData>"),str);
            data[stringvars[i]] = map(x -> x[1], lm);
        end

        # Parse all the numeric variables
        for i = 1:length(numvars)
            lm = eachmatch(Regex("<SimpleData name=\"$(numvars[i])\">(.*?)</SimpleData>"),str);
            data[numvars[i]] = map(x -> parse(Float64, x[1]), lm)
        end

        needsNaNs = ["alcorr", "be10ep", "beprod", "eal_err", "elev_std", "errbe_prod", "al26ep", "alprod", "be10ep_err", 
            "beself", "eal_gcmyr", "erral_ams", "errbe_tot", "al26ep_err", "alself", "be10nc", "besnow", "eal_mmkyr", 
            "erral_muon", "sizemax", "al26nc", "alsnow", "be10nc_err", "betopo", "ebe_err", "erral_prod", "sizemin", 
            "al26nc_err", "altopo", "be10np", "betots", "ebe_gcmyr", "erral_tot", "slp_ave", "al26np", "altots", "be10np_err", 
            "ebe_mmkyr", "errbe_ams", "slp_std", "al26np_err", "area", "becorr", "dbver", "elev_ave", "errbe_muon", ]

        for i = 1:length(needsNaNs)
            data[needsNaNs[i]][data[needsNaNs[i]] .< 0] .= NaN
        end

        return data
    end


## --- Parse the basin polygon outlines
    """
    ```julia
    basin_polygon_n, basin_polygon_lat, basin_polygon_lon = parse_octopus_polygon_outlines(str, isfirstcoord)
    ```

    Return coordinates `basin_polygon_lat`, `basin_polygon_lon` for each basin, and the 
    number of coordinates in each basin `basin_polygon_n`.

    Requires unparsed OCTOPUS data `str`. Does not combine subbasins.
    """
    function parse_octopus_polygon_outlines(str,isfirstcoord)
        # Preallocate
        nbasins = length(isfirstcoord)  # Does not differentiate basins and subbasins
        basin_polygon_lat = Array{Array}(undef, nbasins)
        basin_polygon_lon = Array{Array}(undef, nbasins)
        basin_polygon_n = Array{Int}(undef, nbasins)

        i = 0
        for m in eachmatch(r"<coordinates>\n\t*(.*?)\n",str)
            i += 1 # Number of parsed values
            parsed = delim_string_function(x -> delim_string_parse(x, ',', Float64), m[1], ' ', Array{Float64,1})
            nparsed = length(parsed)

            basin_polygon_n[i] = nparsed;
            basin_polygon_lon[i] = parsed .|> x -> x[1]
            basin_polygon_lat[i] = parsed .|> x -> x[2]
        end

        return basin_polygon_n, basin_polygon_lat, basin_polygon_lon
    end


## --- Calculate average slope for each basin
    """
    ```julia
    get_basin_srtm15plus_aveslope(srtm::Dict,nbasins,subbasins,basin_polygon_lat,basin_polygon_lon)
    ```

    Calculate average slope and standard deviation of each basin using the STRTM15+ DEM.
    """
    function get_basin_srtm15plus_aveslope(srtm::NamedTuple,
        nbasins, subbasins, basin_polygon_lat, basin_polygon_lon)

        # Define variables
        slope = srtm.slope
        x_lon_cntr = srtm.x_lon_cntr
        y_lat_cntr = srtm.y_lat_cntr

        # Initialize progress bar
        p = Progress(nbasins+1, desc = "Calculating basin slope:")
        next!(p)

        # Preallocate
        basinslope_avg = Array{Float64}(undef, nbasins)
        basinslope_err = Array{Float64}(undef, nbasins)
        basin_N = Array{Int64}(undef, nbasins)
        
        # Average slope of each basin
        for i = 1:nbasins
            # Preallocate
            rowsinbasin = Array{Int}(undef, 0)
            columnsinbasin = Array{Int}(undef, 0)

            # Indices of SRTM15+ grid points in each subbasin
            for j = 1:length(subbasins[i])
                k = subbasins[i][j]
                subbasin_lon = basin_polygon_lon[k]
                subbasin_lat = basin_polygon_lat[k]
                (column, row) = find_grid_inpolygon(x_lon_cntr, y_lat_cntr, subbasin_lon, subbasin_lat)
                rowsinbasin = vcat(rowsinbasin,row)
                columnsinbasin = vcat(columnsinbasin,column)
            end

            # Slopes of SRTM15+ grid points in basin i
            pointsinbasin = unique(hcat(rowsinbasin,columnsinbasin),dims=1);
            basin_slopes = Array{UInt16}(undef, length(rowsinbasin))
            for j=1:size(pointsinbasin,1)
                basin_slopes[j] = slope[pointsinbasin[j,1],pointsinbasin[j,2]]
            end

            # Average all the slopes
            basinslope_avg[i] = mean(basin_slopes)
            basinslope_err[i] = std(basin_slopes)
            basin_N[i] = length(basin_slopes)

            # Bump progress meter
            next!(p)
        end

        return basinslope_avg, basinslope_err
    end


## --- Curve fitting functions for slope / erosion rate
    """
    ```julia
    linear(x,p)
    ```
    """
    function linear(x,p)
        y = p[1] .+ x * p[2]
    end


## --- Calculate precipitation
    """
    ```julia
    find_precip(lat::AbstractArray, lon::AbstractArray, precip::AbstractMatrix)
    ```
    
    Find the precipitation at the coordinates `lat` and `lon`.
    """
    function find_precip(lat::AbstractArray, lon::AbstractArray, precip::AbstractMatrix)
        # Preallocate
        out = Array{Float64}(undef,size(lat))

        # Index into the correct position in the precipitation array
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i] > 88.542 || lat[i] < -88.542
                out[i] = NaN
            else
                row = trunc(Int, mod(lon[i]*192/360 + 0.5, 192.0)) + 1
                col = trunc(Int, 46.99999 - lat[i]*94/177.084) + 1
                out[i] = precip[row, col]
            end
        end

        return out
    end


## --- Functions for dealing with SRTM15
    resourcepath = "data"

    # Read srtm15plus file from HDF5 storage, downloading from cloud if necessary
    function get_srtm15plus_aveslope(varname="")
        # Available variable names: "slope", "y_lat_cntr", "x_lon_cntr",
        # "nanval", "cellsize", "scalefactor", and "reference". Units are
        # meters of elevation and decimal degrees of latitude and longitude

        # Construct file path
        filepath = joinpath(resourcepath,"srtm15plus_aveslope.h5")

        # Download HDF5 file from Google Cloud if necessary
        if ~isfile(filepath)
            @info "Downloading srtm15plus.h5 from google cloud storage\n"
            download("https://storage.googleapis.com/statgeochem/srtm15plus_aveslope.h5", filepath)
        end

        # Read and return the file
        return h5read(filepath, "vars/"*varname)
    end

    # Find the elevation of points at position (lat,lon) on the surface of the
    # Earth, using the SRTM15plus 15-arc-second elevation model.
    function find_srtm15plus_aveslope(srtm15plus,lat,lon)

        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        elseif isa(srtm15plus,Dict)
            data = srtm15plus["slope"]
        elseif isa(srtm15plus, Array)
            data = srtm15plus
        else
            error("wrong srtm15plus variable")
        end

        # Scale factor (cells per degree) = 60 * 4 = 240
        # (15 arc seconds goes into 1 arc degree 240 times)
        sf = 240

        # Create and fill output vector
        out=Array{Float64}(size(lat));
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                out[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                # Note that STRTM15 plus has N+1 columns where N = 360*sf
                row = 1 + round(Int,(90+lat[i])*sf)
                col = 1 + round(Int,(180+lon[i])*sf)
                # Find result by indexing
                res = data[row,col]
                if res > 1000
                    out[i] = NaN
                else
                    out[i] = res
                end
            end
        end

        return out
    end

    function find_srtm15plus_aveslope_around(srtm15plus,lat,lon; halfwidth=1::Integer, max_allowed_slope=1000::Number)

        # Interpret user input
        if length(lat) != length(lon)
            error("lat and lon must be equal length\n")
        elseif isa(srtm15plus,Dict)
            data = srtm15plus["slope"]
        elseif isa(srtm15plus, Array)
            data = srtm15plus
        else
            error("wrong srtm15plus variable")
        end

        # Scale factor (cells per degree) = 60 * 4 = 240
        # (15 arc seconds goes into 1 arc degree 240 times)
        sf = 240

        # Create and fill output vector
        out=Array{Float64}(size(lat));
        for i=1:length(lat)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                out[i] = NaN
            else
                # Convert latitude and longitude into indicies of the elevation map array
                # Note that STRTM15 plus has N+1 columns where N = 360*sf
                row = 1 + round(Int,(90+lat[i])*sf)
                col = 1 + round(Int,(180+lon[i])*sf)

                # Find result by indexing
                k = 0;
                out[i] = 0;
                for r = (row-halfwidth):(row+halfwidth)
                    for c = (col-halfwidth):(col+halfwidth)
                        if r>0 && r<43202
                            res = data[r,mod(c-1,86400)+1]
                            if res < max_allowed_slope
                                k +=1
                                out[i] += res
                            end
                        end
                    end
                end
                out[i] /= k
            end
        end

        return out
    end


## --- Bin means with percentile bin edges
    """
    ```julia
    binmeans_percentile(x::AbstractArray, y::AbstractArray; step::Number=5)
    ```

    The means of `y` binned by x into bins equally spaced by percentile. Returns bin centers,
    means, and standard deviations for each bin.
    """
    function binmeans_percentile(x::AbstractArray, y::AbstractArray; step::Number=5)
        # Get bin edges and centers in terms of percentiles
        binedges = collect(0:step:100)
        binedges[end] != 100 && push!(binedges, 100)

        # Sort data
        p = sortperm(x)
        x = x[p]
        y = y[p]

        # Get percentile indices and bin centers
        npoints = length(y)
        indices = round.(Int, binedges / 100 * npoints)     # Percentile indices
        indices[1] = 1                                      # Make indices index-able
        indx = [x[i] for i in indices]
        bincenters = [(indx[i-1]+indx[i])/2 for i in 2:lastindex(indx)]
        
        # Get means for each percentile
        μ = Array{Float64}(undef, length(bincenters))
        σ = Array{Float64}(undef, length(bincenters))
        for i = 2:lastindex(indices)
            μ[i-1] = nanmean(y[indices[i-1]:indices[i]])
            σ[i-1] = nanstd(y[indices[i-1]:indices[i]])
        end

        return bincenters, μ, σ
    end

## -- End of file