## Rho's edits to RockWeatheringFluxProject.jl but I don't have version control yet

## --- Import the packages and utilities we'll neeed
# External packages
using ProgressMeter: @showprogress
using Plots
using JLD, HDF5
using StatGeochem
using HTTP, JSON

# Local utilities
include("RockWeatheringFluxUtilities.jl")


## -- Generate random points on Earth
npoints = 500;
rocklat = Array{Float64}(undef, npoints)
rocklon = Array{Float64}(undef, npoints)
elevations = Array{Float64}(undef, npoints)
etopo = get_etopo("elevation")

function bar(rocklat, rocklon, elevations, etopo)
    i = 1   # Generic counter

    while i < (npoints + 1)
        # Generate a random point
        (randlat, randlon) = randlatlon()

        # Find the elevation
        elev = find_etopoelev(etopo, randlat, randlon)

        # If the point is above sea level, add it to the list of points on exposed crust
        if elev[1] > 0
            rocklat[i] = randlat
            rocklon[i] = randlon
            elevations[i] = elev[1]

            i += 1
        end
    end
end

# this is actually ~10x faster, as it turns out
rocklat = Array{Float64}(undef, 0)
rocklon = Array{Float64}(undef, 0)

function foo(rocklat, rocklon, etopo)
    while length(rocklat) < npoints

        # Generate some random latitudes and longitudes with uniform
        #  spatial density on the globe
        (randlat, randlon) = random_lat_lon(npoints)

        # Find which points are above sea level
        elevations = find_etopoelev(etopo,randlat,randlon)
        abovesea = elevations .> 0

        # Concatenate together all the points that represent exposed crust
        rocklat = vcat(rocklat,randlat[abovesea])
        rocklon = vcat(rocklon,randlon[abovesea])
    end

    rocklat = rocklat[1:npoints]
    rocklon = rocklon[1:npoints]
    elevations = find_etopoelev(etopo,rocklat,rocklon)
end