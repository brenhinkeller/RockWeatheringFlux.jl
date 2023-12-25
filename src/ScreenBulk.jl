# This file is a wrapper for ScreenBulkBase.jl, which reads bulk.mat and bulktext.mat,
# converts all units to wt.%, restricts data to whole-rock analyses, and saves the data as
# a new file.

# Using a wrapper allows me to call ScreenBulkBase.jl in sensitivity tests, meaning I can 
# use the same code to run my tests as I do in the actual computations. This should reduce
# unintended errors from editing or not updating the simulation code!

## --- Set up
    # Packages
    using RockWeatheringFlux
    using MAT, HDF5, Dates
    using LoopVectorization: @turbo
    
    # File name IO
    fileout = "output/bulk.h5"

    # Start timer 
    start = now()
    @info """ Start: $(Dates.Date(start)) $(Dates.format(start, "HH:MM"))
    Output: $fileout
    """

    # Calculate a reasonable assumption for wt.% volatiles
    # lim = (12.01+2*16)/(40.08+12.01+16*3)*100                     # Limestone
    # mag = (12.01+2*16)/(24.869+12.01+16*3)*100                    # Magnesite
    dol = (12.01+2*16)/((24.869+40.08)/2+12.01+16*3)*100          # Dolomite
    gyp = (32.07+16*3+2*(18))/(40.08+32.07+16*4+2*(18))*100       # Gypsum
    bas = (32.07+16*3+0.5*(18))/(40.08+32.07+16*4+0.5*(18))*100   # Bassanite (2CaSO₄⋅H₂O)

    # Call file
    include("ScreenBulkBase.jl")

    # End timer
    stop = now()
    @info """
    Stop: $(Dates.Date(stop)) $(Dates.format(stop, "HH:MM")).
    Program runtime: $(canonicalize(round(stop - start, Dates.Minute))).
    """

## --- End of File

# TO DO: 
    # Re-read in old files, and figure out what to do for parsing / correcting element 
        # oxide weights. For example, some samples may have different values for CaO and 
        # CaCO3, but Mg and MgO may represent the same data.
    
# Read in raw data files. Useful code:

# # Filter ages younger than 0 or greater than the age of the earth
# invalid_age = vcat(findall(>(4000), bulk.Age), findall(<(0), bulk.Age))
# bulk.Age[invalid_age] .= NaN

# # Fill in any missing ages from bounds
# for i in eachindex(bulk.Age)
#     if isnan(bulk.Age[i])
#         bulk.Age[i] = nanmean([bulk.Age_Max[i], bulk.Age_Min[i]])
#     end
# end