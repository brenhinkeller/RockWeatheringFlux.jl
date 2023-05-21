## --- Setup
    # External packages
    using StatGeochem
    using DelimitedFiles
    using ProgressMeter
    using LoopVectorization
    # using Dates

    # # File parsing packages
    # using JLD
    using HDF5
    # using HTTP
    # using JSON
    using MAT

    # Local utilities
    include("Utilities.jl")


## --- Load pre-generated Macrostrat data
    @info "Loading Macrostrat data"

    # Load and match
    macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)
    @time macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    # Reduce likelihood of double matches
    macro_cats.sed .&= .! macro_cats.cover .& .! macro_cats.ign		# Sed excludes cover and igns
    macro_cats.ign .&= .! macro_cats.cover							# Ign excludes cover
    macro_cats.met .&= .! macro_cats.cover							# Met excludes cover
    macro_cats.cryst .&= .! macro_cats.metased                      # Cryst excludes metased

    # Figure out how many data points weren't matched
    known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    total_known = count(known_rocks)

    matched = known_rocks .| macro_cats.cover
    not_matched = .!matched
    multi_matched = ((macro_cats.sed .& macro_cats.ign) .| (macro_cats.sed .& macro_cats.met) 
        .| (macro_cats.ign .& macro_cats.met)
    )
    
    # Print to terminal
    # Conflicting matches are ok -- 
    @info "Macrostrat parsing complete!
      not matched = $(count(not_matched))
      conflicting = $(count(multi_matched))
      sed and ign = $(count(macro_cats.sed .& macro_cats.ign))
      sed and met = $(count(macro_cats.sed .& macro_cats.met))
      ign and met = $(count(macro_cats.ign .& macro_cats.met))"


## --- Load EarthChem data
    @info "Loading EarthChem data"

    # Bulk data
    bulk_raw = matopen("data/bulk.mat")
    bulk_dict = read(bulk_raw, "bulk")
    close(bulk_raw)
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

    # Match codes to rock types
    bulk_cats = match_earthchem(bulk.Type)

    # Load indices of matched samples from samplematch.jl
    bulkidx = Int.(vec(readdlm("output/matched_bulkidx.tsv")))
    

## --- Calculate erosion rate at each coordinate point of interest
	# @info "Calculating slope and erosion at each point"
	
    # Load the slope variable from the SRTM15+ maxslope file
    # srtm15_slope = h5read("data/srtm15plus_maxslope.h5", "vars/slope")
    # srtm15_sf = h5read("data/srtm15plus_maxslope.h5", "vars/scalefactor")

    # # Get slope at each coordinate point
    # rockslope = avg_over_area(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
    #     srtm15_sf, halfwidth=7
    # )

    # # Calculate all erosion rates (mm/kyr)
    # # TO DO: update this function with a better erosion estimate
    # rock_ersn = emmkyr.(rockslope)


## --- Preallocate
    allkeys = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal, :sed,
        :volc, :plut, :ign,
        :metased, :metaign, :met,
        :cryst
    )
    allinitvals = fill(NaN, length(allkeys))

    # Dictionaries
    erosion = Dict(zip(allkeys, allinitvals))
    crustal_area = Dict(zip(allkeys, allinitvals))

    # Declare constants
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    const crustal_density = 2750              # Average crustal density (kg/m³)

## --- Calculate
    # Erosion (m/Myr)
    for i in eachidex(allkeys)
        erosion[allkeys[i]] = nanmean(rock_ersn[macro_cats[i]])
    end
    erosion = NamedTuple{Tuple(allkeys)}(values(erosion))

    # Crustal area (m²), assume proportional distribution of rock types under cover
    for i in eachidex(allkeys)
        crustal_area[allkeys[i]] = count(macro_cats[i]) / total_known * contl_area
    end
    crustal_area = NamedTuple{Tuple(allkeys)}(values(crustal_area))
    

## --- Get average sed contributions from EarthChem data
    # p_wt = (
    #     alluvium = [0.0], siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], 
    #         evaporite = [0.0], phosphorite = [0.0], coal = [0.0], volcaniclast = [0.0], 
    #         sed = [0.0],
    #     volc = [0.0], plut = [0.0], ign = [0.0],
    #     metased = [0.0], metaign = [0.0], met = [0.0],
    # )
    # for i in eachindex(p_wt)
    #     p_wt[i][1] = nanmean(bulk.P2O5[bulk_cats[i]]) 
    # end

    # # Calculate P flux by source contributions in kg/yr
    # pflux_source = (
    #     siliciclast = [0.0], shale = [0.0], carb = [0.0], chert = [0.0], evaporite = [0.0], 
    #         coal = [0.0], sed = [0.0],
    #     volc = [0.0], plut = [0.0], ign = [0.0],
    #     metased = [0.0], metaign = [0.0], met = [0.0],
    # )
    # for i in eachindex(pflux_source)
    #     pflux_source[i][1] = macro_ersn[i][1] * p_wt[i][1] * crustal_area[i][1] * crustal_density * 1e-6
    # end
    # pflux_global = pflux_source.sed[1] + pflux_source.ign[1] + pflux_source.met[1]

## --- End of File
