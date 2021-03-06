using StatGeochem
using HDF5

# Load dataset as dict
print("Loading SRTM\n")
@time srtm = get_srtm15plus()

# # Calculate average slope
# print("Calculating slope\n")
# @time slope = ave_slope_earth(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], srtm["cellsize"], minmatval=-12000, maxmatval=9000)
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


# Calculate maximum slope
print("Calculating slope\n")
@time slope = maxslope(srtm["elevation"], srtm["x_lon_cntr"], srtm["y_lat_cntr"], srtm["cellsize"], minmatval=-12000)

# Save results
fid = h5open("srtm15plus_maxslope.h5","w");
g = g_create(fid, "vars");
print("Saving to HDF5:\n")
g["y_lat_cntr"] = srtm["y_lat_cntr"]
g["x_lon_cntr"] = srtm["x_lon_cntr"]
g["cellsize"] = srtm["cellsize"]
g["scalefactor"] = srtm["scalefactor"]
@time g["slope","compress",3] = slope;
close(fid)
