## --- Setup
    # External packages
    using StatGeochem
    using DelimitedFiles
    # using Dates

    # # File parsing packages
    # using JLD
    # using HDF5
    # using HTTP
    # using JSON
    using MAT

    # Local utilities
    include("Utilities.jl")

## --- Load pre-generated Macrostrat data
    macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)
    macro_cats = match_rocktype(rocktype, rockname, rockdescrip, major=false)

    # Reduce likelihood of double matches
    macro_cats.sed .&= .! macro_cats.cover .& .! macro_cats.ign		# Sed excludes cover and igns
    macro_cats.ign .&= .! macro_cats.cover							# Ign excludes cover
    macro_cats.met .&= .! macro_cats.cover							# Met excludes cover
    
    # Define crystalline rocks as igneous or non metased metamorphic
    macro_cryst = macro_cats.ign .| macro_cats.metaign .| macro_cats.met 

    # Figure out how many data points weren't matched
    known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    total_known = count(known_rocks)

    matched = known_rocks .| macro_cats.cover
    not_matched = .!matched
    multi_matched = ((macro_cats.sed .& macro_cats.ign) .| (macro_cats.sed .& macro_cats.met) 
        .| (macro_cats.ign .& macro_cats.met)
    )
    
    # Print to terminal
    @info "Macrostrat parsing complete!
    not matched = $(count(not_matched)), conflicting matches of major rocktypes = $(count(multi_matched))\n"

## --- Load EarthChem data
    # Bulk data
    bulk_raw = matopen("data/bulk.mat")
    bulk_dict = read(bulk_raw, "bulk")
    close(bulk_raw)
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

    # Matched samples from samplematch.jl
    bulkidx = readdlm("output/matched_bulkidx.tsv")
    

## --- Calculate erosion rate at each coordinate point of interest
	@info "Calculating slope and erosion at each point"
	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("data/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("data/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point
    rockslope = avg_over_area(srtm15_slope, rocklat, rocklon, srtm15_sf, halfwidth=7)

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)
    
    @info "Slope and erosion calculated successfully."


## --- Get erosion (m/Myr) for each rock type in Macrostrat
    macro_ersn = (
        siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], evaporite = [0.0], 
            coal = [0.0], sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
    )
    for i in eachindex(macro_ersn)
        macro_ersn[i][1] = nanmean(rock_ersn[macro_cats[i]])
    end


## --- Get crustal area for each rock type in Macrostrat
    # Constants
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    const crustal_density = 2750              # Average crustal density (kg/m³)

    # Area (in m²) of each rock type. 
    # Using all rocks, assuming equal distribution of rock types under cover
    crustal_area = (
        siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], evaporite = [0.0], 
            coal = [0.0], sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
        cover = [0.0]
    )
    for i in eachindex(crustal_area)
        crustal_area[i][1] = count(macro_cats[i]) / total_known * contl_area
    end

    

## --- Get average sed contributions from EarthChem data
    bulk_cats = match_earthchem(bulk.Type, major=false)

    p_wt = (
        alluvium = [0.0], siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], 
            evaporite = [0.0], phosphorite = [0.0], coal = [0.0], volcaniclast = [0.0], 
            sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
    )
    for i in eachindex(p_wt)
        p_wt[i][1] = nanmean(bulk.P2O5[bulk_cats[i]]) 
    end

    # Calculate P flux by source contributions in kg/yr
    pflux_source = (
        siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], evaporite = [0.0], 
            coal = [0.0], sed = [0.0],
        volc = [0.0], plut = [0.0], ign = [0.0],
        metased = [0.0], metaign = [0.0], met = [0.0],
    )
    for i in eachindex(pflux_source)
        pflux_source[i][1] = macro_ersn[i][1] * p_wt[i][1] * crustal_area[i][1] * crustal_density * 1e-6
    end
    pflux_global = pflux_source.sed[1] + pflux_source.ign[1] + pflux_source.met[1]

## --- End of File
