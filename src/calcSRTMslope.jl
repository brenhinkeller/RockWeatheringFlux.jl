## --- Load packages we'll be using and load the SRTM15+ dataset
    using StatGeochem
    using HDF5
    using Dates

"""
```julia
```
Calculate the maximum slope for each point in the SRTM15+ data set. Generates a file that lives in
the resource folder.

This doesn't really make sense in utilities.jl because it's just to generate one of the data files.
Something to think about when I do repo organization...

Also I know I just moved this to a function but in retrospect. Not sure I actually like that!
"""
function max_srtm_slope()
    # Get SRTM15+ file
    @info "Loading SRTM\n"
    srtm = get_srtm15plus()

    # Calculate the maximum slope for each point in the data set
    @info "Calculating slope. This may take up to 30 minutes. Started $(Dates.format(now(), "HH:MM"))"
    slope = maxslope(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], srtm["cellsize"], minmatval=-12000)

    # Save results
    @info "Saving slope to HDF5 file"
    filename = "srtm15plus_maxslope"
    fid = h5open("data/$filename.h5","w")
    g = create_group(fid, "vars")

    # Copy over SRTM15+ location data
    g["y_lat_cntr"] = srtm["y_lat_cntr"]
    g["x_lon_cntr"] = srtm["x_lon_cntr"]
    g["cellsize"] = srtm["cellsize"]
    g["scalefactor"] = srtm["scalefactor"]

    # Add a data set for slope and compress data
    @time g["slope", compress=3] = slope
    close(fid)
end



