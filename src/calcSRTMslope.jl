## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using Dates

    # Get SRTM15+ file
    @info "Loading SRTM"
    srtm = get_srtm15plus()

## --- Calculate maximum slope and save data set
    @info "Calculating slope. This may take up to 30 minutes. Started $(Dates.format(now(), "HH:MM"))"
    slope = maxslope(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], srtm["cellsize"], minmatval=-12000)

    # Save results to the data folder
    @info "Saving slope to HDF5 file"
    filename = "srtm15plus_maxslope"
    fid = h5open("data/$filename.h5","w")
    g = create_group(fid, "vars")

    # Copy over SRTM15+ location data
    g["y_lat_cntr"] = srtm["y_lat_cntr"]
    g["x_lon_cntr"] = srtm["x_lon_cntr"]
    g["cellsize"] = srtm["cellsize"]
    g["scalefactor"] = srtm["scalefactor"]

    # Add a data set for slope and compress data (Takes about 2.5 minutes)
    @time g["slope", compress=3] = slope
    close(fid)

## --- End of file



