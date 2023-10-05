## --- Set up
    # Packages
    using StatGeochem
    using HDF5
    using Dates

    # Get SRTM15+ file
    @info "Loading SRTM"
    srtm = h5read("data/srtm15plus.h5", "vars/")
    # try
    #     srtm = h5read("data/srtm15plus.h5", "vars/")
    # catch
    #     # Note that this will download the SRTM15+ file to \resources, not \data
    #     srtm = get_srtm15plus()
    # end


## --- Calculate maximum slope and save data set
    @info "Maximum slope calculation started $(Dates.format(now(), "HH:MM"))"
    @time slope = maxslope(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], srtm["cellsize"], 
        minmatval=-12000)

    # Save results to the data folder
    @info "Saving slope to HDF5 file"
    filename = "srtm15plus_maxslope"
    fid = h5open("output/$filename.h5","w")
    g = create_group(fid, "vars")

    # Copy over SRTM15+ location data
    g["y_lat_cntr"] = srtm["y_lat_cntr"]
    g["x_lon_cntr"] = srtm["x_lon_cntr"]
    g["cellsize"] = srtm["cellsize"]
    g["scalefactor"] = srtm["scalefactor"]

    # Add a data set for slope and compress data (Takes about 2.5 minutes)
    g["slope", compress=3] = slope
    close(fid)


## --- End of file