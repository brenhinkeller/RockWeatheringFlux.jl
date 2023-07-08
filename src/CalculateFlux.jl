## --- Setup
    # Packages
    using StatGeochem
    using DelimitedFiles
    using ProgressMeter
    using Measurements
    using HDF5
    using MAT

    # Local utilities
    include("Utilities.jl")
    include("NaNMeasurements.jl")


## --- Load pre-generated Macrostrat data
    @info "Loading Macrostrat data"

    # Load and match
    macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)
    @time macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    # Figure out how many data points weren't matched
    known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    total_known = count(known_rocks)

    matched = known_rocks .| macro_cats.cover
    not_matched = .!matched
    multi_matched = ((macro_cats.sed .& macro_cats.ign) .| (macro_cats.sed .& macro_cats.met) 
        .| (macro_cats.ign .& macro_cats.met)
    )
    
    # Print to terminal
    @info """
    Macrostrat parsing complete!
    not matched = $(count(not_matched))
    """

    if count(multi_matched) > 0
        @warn """
        $(count(multi_matched)) conflicting matches present
        sed and ign = $(count(macro_cats.sed .& macro_cats.ign))
        sed and met = $(count(macro_cats.sed .& macro_cats.met))
        ign and met = $(count(macro_cats.ign .& macro_cats.met))
        """
    end
    

## --- Calculate erosion rate at each coordinate point of interest
	@info "Calculating slope and erosion at each point"
	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("data/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("data/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point
    # Modify this function to return an error as well
    rockslope = avg_over_area(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, halfwidth=7
    )

    # Calculate all erosion rates (mm/kyr)
    rock_ersn = emmkyr.(rockslope)


## --- Get erosion and continental area for each rock type
    # Preallocate
    allkeys = keys(macro_cats)
    allinitvals = fill(NaN, length(allkeys))

    erosion = Dict(zip(allkeys, allinitvals .± allinitvals))
    crustal_area = Dict(zip(allkeys, allinitvals))

    # Average erosion by rock type (m/Myr)
    for i in allkeys
        erosion[i] = nanmean(rock_ersn[macro_cats[i]])
    end
    erosion = NamedTuple{Tuple(keys(erosion))}(values(erosion))

    # Crustal area (m²), assume proportional distribution of rock types under cover
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    for i in allkeys
        crustal_area[i] = count(macro_cats[i]) / total_known * contl_area
    end
    crustal_area = NamedTuple{Tuple(keys(crustal_area))}(values(crustal_area))


## --- Load EarthChem data
    @info "Loading EarthChem data"
    bulk = matread("data/bulk_newunits.mat")
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))

    # Get rock types from codes
    bulk_cats = match_earthchem(bulk.Type)

    # Load indices of matched samples from samplematch.jl
    bulkidx = Int.(vec(readdlm("output/matched_bulkidx.tsv")))


## --- Every element of interest from EarthChem
    majors, minors = get_elements()
    allelements = [majors; minors]
    nelements = length(allelements)
    # strallelements = string.(allelements)

    npoints = length(macrostrat.rocktype)
    subcats = collect(allkeys)
    deleteat!(subcats, findall(x->x==:cover,subcats))       # Do not compute cover


# ## --- Create HDF5 file to store results for eroded material
#     fid = h5open("output/erodedmaterial.h5", "w")

#     # Metadata
#     write(fid, "element_names", strallelements)                           # Names of analyzed elements
#     write(fid, "npoints", npoints)                                    # Total macrostrat samples
#     nsample = create_dataset(fid, "nsamples_byelem", Int, (nelements,))   # Non-NaN samples for each element

#     # Bulk rock global flux
#     bulkrockflux = create_group(fid, "bulkrockflux")
#     write(bulkrockflux, "rocktypes", string.(subcats))
#     bulkflux_val = create_dataset(bulkrockflux, "val", Float64, (length(subcats),))
#     bulkflux_std = create_dataset(bulkrockflux, "std", Float64, (length(subcats),))

#     # For each element in allelements, wt.% and flux by rock subtype, and total global flux
#     elementflux = create_group(fid, "elementflux")

#     # Global flux
#     totalflux = create_group(elementflux, "totalelemflux")
#     totalflux_val = create_dataset(totalflux, "val", Float64, (nelements,))
#     totalflux_std = create_dataset(totalflux, "std", Float64, (nelements,))

#     # Separated by rock subtypes
#     byrocktype = create_group(elementflux, "byrocktype")
#     for i in subcats
#         typegroup = create_group(byrocktype, string(i))

#         create_dataset(typegroup, "wtpct_val", Float64, (nelements,))
#         create_dataset(typegroup, "wtpct_std", Float64, (nelements,))

#         create_dataset(typegroup, "flux_val", Float64, (nelements,))
#         create_dataset(typegroup, "flux_std", Float64, (nelements,))
#     end


# ## --- Calculate undifferentiated (bulk) flux
#     const crustal_density = 2750                                # kg/m³

#     # Keys for erosion and crustal_area are in an arbitrary order. Correct indexing critical
#     bulkflux = Dict(zip(subcats, fill(NaN ± NaN, length(subcats))))
#     for i in subcats
#         bulkflux[i] = erosion[i] * crustal_area[i] * crustal_density * 1e-6
#     end
#     bulkflux = NamedTuple{Tuple(keys(bulkflux))}(values(bulkflux))

#     # Save to file. Keys are again in arbitrary order
#     for i in eachindex(subcats)
#         sub_i = subcats[i]
#         bulkflux_val[i] = bulkflux[sub_i].val
#         bulkflux_std[i] = bulkflux[sub_i].err
#     end

    
# ## --- Calculate flux by element / element oxide
#     # TO DO: I should make sure these aren't also scrambled...
#     for i in eachindex(allelements)
#         # Calculate wt.%, flux, and global flux of each element
#         wt, flux, global_flux, n = flux_source(bulk[allelements[i]], bulkidx, erosion, macro_cats, 
#             crustal_area, elem=strallelements[i]
#         )

#         # Write data to file
#         nsample[i] = n
#         totalflux_val[i] = global_flux.val
#         totalflux_std[i] = global_flux.err

#         for j in eachindex(subcats)
#             typegroup = string(subcats[j])

#             byrocktype[typegroup]["flux_val"][i] = flux[subcats[j]].val
#             byrocktype[typegroup]["flux_std"][i] = flux[subcats[j]].err
            
#             byrocktype[typegroup]["wtpct_val"][i] = wt[subcats[j]].val
#             byrocktype[typegroup]["wtpct_std"][i] = wt[subcats[j]].err
#         end
#     end

#     close(fid)


## --- Calculate denundation at each point
    # Declare constants
    const crustal_density = 2750                                # kg/m³
    const unit_sample_area = (148940000 * 1000000) / npoints    # m²
    # const kg_to_gt = 1e12                                       # Conversion factor

    # Create file to save data
    fid = h5open("output/erodedmaterial_new.h5", "w")

    # Denundation at each point
    sampleflux = Array{Measurement{Float64}}(undef, npoints, 1)
    for i in eachindex(sampleflux)
        sampleflux[i] = rock_ersn[i] * unit_sample_area * crustal_density * 1e-6
    end

    # Save to file
    sampleflux_val, sampleflux_err = unmeasurementify(sampleflux)
    bulk_denundation = create_group(fid, "bulk_denundation")
    write(bulk_denundation, "values", sampleflux_val)
    write(bulk_denundation, "errors", sampleflux_err)


## --- Calculate flux of each element at each point
    # Preallocate file space
    element_flux = create_group(fid, "element_flux")
    elem_vals = create_dataset(element_flux, "values", Float64, (npoints, nelements))
    elem_errs = create_dataset(element_flux, "errors", Float64, (npoints, nelements))
    elem_head = create_dataset(element_flux, "header", string.(allelements))

    # Use rock type wt.% averages
    for i in eachindex(allelements)
        # Compute average wt.% for element i by rock type
        avgwt = Dict(zip(subcats, fill(NaN ± NaN, length(subcats))))
        for j in collect(keys(avgwt))
            avgwt[j] = nanmean(bulk[allelements[i]][bulk_cats[j]]) ± nanstd(bulk[allelements[i]][bulk_cats[j]])
        end

        # Filter rocks of each rock type and compute flux of element i
        elementflux = zeros(Measurement{Float64}, npoints)
        for j in collect(keys(avgwt))
            filter = macro_cats[j]
            for k in eachindex(filter)
                filter[k] && (elementflux[k] = sampleflux[k] * avgwt[j] * 1e-2)
            end
        end

        # Save to file
        elementflux_val, elementflux_err = unmeasurementify(elementflux)
        elem_vals[:,i] = elementflux_val
        elem_errs[:,i] = elementflux_err
    end

    close(fid)

## --- End of File
