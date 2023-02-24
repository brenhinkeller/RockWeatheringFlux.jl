## --- Load packages we'll be using and load the SRTM15+ dataset
    using StatGeochem
    using HDF5

    @info "Loading SRTM\n"
    @time srtm = get_srtm15plus()


## --- Calculate maximum slope
    @info "Calculating slope\n"
    @time slope = maxslope(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], srtm["cellsize"], minmatval=-12000)

    # Save results
    @info "Saving to HDF5 file\n"
    fid = h5open("data/srtm15plus_maxslope.h5","w")
    g = create_group(fid, "vars")

    # Copy over SRTM15+ location data
    g["y_lat_cntr"] = srtm["y_lat_cntr"]
    g["x_lon_cntr"] = srtm["x_lon_cntr"]
    g["cellsize"] = srtm["cellsize"]
    g["scalefactor"] = srtm["scalefactor"]

    # Add a dataset for slope and compress data
    @time g["slope", compress=3] = slope        # thoughts on chunking in addition to compressing? could improve efficiency...
    close(fid)


# # Calculate average slope
# print("Calculating slope\n")
# @time slope = aveslope(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], srtm["cellsize"], minmatval=-12000, maxmatval=9000)
#
# # Save results
# fid = h5open("srtm15plus_aveslope.h5","w");
# g = g_create(fid, "vars");
# print("Saving to HDF5:\n")
# g["y_lat_cntr"] = srtm["y_lat_cntr"]
# g["x_lon_cntr"] = srtm["x_lon_cntr"]
# g["cellsize"] = srtm["cellsize"]
# g["scalefactor"] = srtm["scalefactor"]
# @time g["slope","compress",3] = slope;
# close(fid)



