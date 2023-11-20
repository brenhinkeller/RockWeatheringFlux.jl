# Slope and flux computations

## --- Slope and erosion rate relationship

    """
    ```julia
    emmkyr(slp)
    ```

    Find erosion rate in mm/kyr given slope `slp`.
    """
    emmkyr(slp) = exp(slp * (0.0091 ± 0.0095) + (3.1 ± 1.9))

    # Previously 10^(slp*0.00567517 + 0.971075)


## --- Calculate wt.% and flux by rock type

    # """
    # ```julia
    # function flux_source(bulk::AbstractArray, bulkidx::Vector{Int64}, erosion::NamedTuple, 
    #     macro_cats::NamedTuple, crustal_area::NamedTuple; 
    #     unitcodes::AbstractMatrix, unitdecoder::AbstractMatrix, crustal_density::Number=2750, 
    #     elem::String="")
    # ```

    # For a specified element in `bulk`, calculate the average wt.% and flux (kg/yr) by rocktype. 
    # Calculate the total global flux of that element (kg/yr). Return the number of samples `n`.

    # Note that `erosion`, `macro_cats`, and `crustal_area` _must_ contain at minimum the keys:
    # ```
    # :siliciclast, :shale, :carb, :chert, :evaporite, :coal, :sed, :volc, :plut, :ign, :metased, 
    # :metaign, :met
    # ```
    # Keys must be type `Symbol`.

    # ### Optional Keyword Arguments
    # - `crustal_density::Number=2750`: Average crustal density (kg/m³).
    # - `elem::String=""`: Element being analyzed, for terminal printout.

    # # Example
    # ```julia-repl
    # julia> wt, flux, global_flux, n = flux_source(bulk.P2O5, bulkidx, erosion, macro_cats, crustal_area, elem="phosphorus")
    # [ Info: 44307 of 50000 phosphorus samples (1%) are not NaN
    # ```
    # """
    # function flux_source(bulk::AbstractArray, bulkidx::Vector{Int64}, erosion::NamedTuple, 
    #     macro_cats::NamedTuple, crustal_area::NamedTuple; 
    #     crustal_density::Number=2750, elem::String="", printinfo=false)

    #     # Preallocate
    #     allkeys = collect(keys(macro_cats))
    #     deleteat!(allkeys, findall(x->x==:cover,allkeys))       # Do not compute cover

    #     allinitvals = fill(NaN ± NaN, length(allkeys))
    #     npoints = length(bulkidx)

    #     wt = Dict(zip(allkeys, allinitvals))
    #     flux = Dict(zip(allkeys, allinitvals))
    #     bulkdata = Array{Float64}(undef, npoints, 1)

    #     # Get EarthChem samples, if present
    #     for i in eachindex(bulkidx)
    #         (bulkidx[i] != 0) ? (bulkdata[i] = bulk[bulkidx[i]]) : (bulkdata[i] = NaN)
    #     end

    #     # Find how many samples have data for the element of interest
    #     n = length(findall(!isnan, bulkdata))
    #     if printinfo
    #         @info "$n of $npoints $elem samples ($(round(n/npoints*100, sigdigits=3))%) are not NaN"
    #     end

    #     # Calculate average wt.% for each rock type
    #     # TO DO: Maybe set no data to 0 instead of NaN? Would require re-writing a bit...
    #     for i in keys(wt)
    #         wt[i] = nanmean(bulkdata[macro_cats[i]]) ± nanstd(bulkdata[macro_cats[i]])
    #     end
    #     wt = NamedTuple{Tuple(keys(wt))}(values(wt))

    #     # Calculate provenance by rock type
    #     for i in keys(flux)
    #         flux[i] = erosion[i] * crustal_area[i] * wt[i] * crustal_density* 1e-8
    #     end
    #     flux = NamedTuple{Tuple(keys(flux))}(values(flux))

    #     # Compute global flux
    #     global_flux = nansum([flux.sed, flux.ign, flux.met])

    #     return wt, flux, global_flux, n
    # end


## --- Find the average value of slope over an area

    """
    ```julia
    function movingwindow(data::AbstractArray, lat::AbstractArray, lon::AbstractArray, 
        sf::Number=240; 
        n::Number=5, maxpossible::Number=0xffff
    )
    ```

    Find the average value of geospatial `data` in a `n` * `n` km window centered at 
    `lat`ᵢ, `lon`ᵢ.

    ### Optional Kwargs:
    * `sf::Number=240`: Scale factor, or cells per degree, for the geospatial `data`. 
        For 15 arc-second resolution, the scale factor is 240, because 15 arc-seconds go 
        into 1 arc-degree 60 (arc minutes / degree) * 4 (cells / arc minute) = 240 times.

    * `n::Number=5`: Defines the size of the window, in kilometers.

    * `maxpossible::Number=0xffff`: The maximum possible value for the variable of interest; 
        variables with values greater than this are ignored.


    """
    function movingwindow(data::AbstractArray, lat::AbstractArray, lon::AbstractArray, 
            sf::Number=240; n::Number=5, maxpossible::Number=0xffff
        )
        # Interpret user input - make sure DEM is the correct size
        (nrows, ncols) = size(data)
        @assert nrows == 180 * sf + 1   # Latitude
        @assert ncols == 360 * sf + 1   # Longitude

        # Preallocate
        out = fill(NaN ± NaN, length(lat))
        # NShalfwidth = Array{Int64}(undef, nrows, 1)   # Constant (for now?)
        EWhalfwidth = Array{Int64}(undef, nrows, 1)   # Depends on latitude

        # Precalculate the number of grid cells per n×n window at each latitude
        km_per_cell = 1852 * 60 / sf / 1000    # Kilometers per cell (lon) at the equator
        target = n/2
        NShalfwidth = Int(round(target / km_per_cell))

        for i = 1:nrows
            latᵢ = 90 - (1/sf*(i-1))
            gridᵢ = cos(deg2rad(latᵢ)) * km_per_cell
            EWhalfwidth[i] = min(Int(round(target / gridᵢ)), ncols)
        end

        # Find result by indexing into the varname matrix
        for i in eachindex(lat, lon)
            if isnan(lat[i]) || isnan(lon[i]) || lat[i]>90 || lat[i]<-90 || lon[i]>180 || lon[i]<-180
                # Result is NaN if either input is NaN or out of bounds
                continue

            else
                # Convert latitude and longitude into indicies 
                row = 1 + round(Int,(90 + lat[i]) * sf)
                col = 1 + round(Int,(180 + lon[i]) * sf)

                row = (row-EWhalfwidth[row]):(row+EWhalfwidth[row])
                col = (col-NShalfwidth):(col+NShalfwidth)

                # Preallocate
                s = fill(NaN, length(row)*length(col))  # Hold the data values we'll use
                k = 0                                   # Counter

                # Index into the array
                for r in row
                    for c in col
                        # Only do the computation if we are in bounds
                        if 1 <= r <= nrows
                            res = data[r, mod(c-1,ncols-1)+1]

                            # Ignore values that are larger than the max possible value
                            if res < maxpossible
                                k +=1
                                s[k] = res
                            end
                        end
                    end
                end

                # Save the average value
                out[i] = (nansum(s) / k) ± (nanstd(s))
            end
        end

        return out
    end


## --- Parse OCTOPUS data

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


## --- DEM and slope computations

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


    # Read srtm15plus file from HDF5 storage, downloading from cloud if necessary
    resourcepath = "data"
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


## --- Bin means with percentile bin edges

    """
    ```julia
    binmeans_percentile(x::AbstractArray, y::AbstractArray; step::Number=5)
    ```

    The means of `y` binned by x into bins equally spaced by percentile. Returns bin centers,
    means, and standard errors of the mean.

    # Example
    ```julia
    (c, m, ex, ey) = binmeans_percentile(x, y, step=5)
    ```
    """
    function binmeans_percentile(x::AbstractArray, y::AbstractArray; step::Number=5)
        @assert length(x) == length(y) "x and y must be equal lengths."

        # Sort data by x value
        p = sortperm(x)
        x = x[p]
        y = y[p]

        # Calculate bin edges and centers. Last bin may be larger than expected
        binedges = collect(0:step:100)
        if binedges[end] != 100 
            binedges[end] = 100
        end
        
        # Calculate the index of each percentile in x, correcting for zero-indexing
        r = ceil.(Int, binedges ./ 100 .* length(x))

        # Calculate bin centers
        nbins = length(binedges) - 1
        bincenters = [(x[r[i]+1] + x[r[i+1]]) / 2 for i = 1:nbins]

        # Calculate means and variances
        μ = zeros(nbins)        # Mean
        σȳ = zeros(nbins)      # Y standard deviation
        σx̄ = zeros(nbins)      # X standard deviation

        # Since we've sorted the data, figuring out what bin each value belongs to is easy!
        for i = 1:nbins
            yᵢ = y[r[i]+1:r[i+1]]
            n = length(yᵢ)
            
            μ[i] = nanmean(yᵢ)
            σȳ[i] = nanstd(yᵢ) / sqrt(n)
            σx̄[i] = nanstd(x[r[i]+1:r[i+1]]) / sqrt(n)
        end

        return bincenters, μ, σx̄, σȳ
    end

    
## --- End of file