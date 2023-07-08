## --- Setup
    # Packages
    using StatGeochem
    using DelimitedFiles
    using ProgressMeter
    using Measurements
    using HDF5
    using MAT
    using DelimitedFiles
    
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
    biglist = [:SiO2,:Al2O3,:Fe2O3T,:TiO2,:MgO,:CaO,:Na2O,:K2O,
        :Ag,:As,:Au,:B,:Ba,:Be,:Bi,:C,:CaCO3,:Cd,:Ce,:Cl,:Co,:Cr2O3,:Cs,:Cu,:Dy,:Er,
        :Eu,:F,:Ga,:Gd,:Hf,:Hg,:Ho,:I,:In,:Ir,:La,:Li,:Lu,:MnO,:Mo,:Nb,:Nd,:NiO,:Os,
        :P2O5,:Pb,:Pd,:Pt,:Pr,:Re,:Rb,:Sb,:Sc,:Se,:S,:Sm,:Sn,:Sr,:Ta,:Tb,:Te,:Th,:Tl,:Tm,
        :U,:V,:W,:Y,:Yb,:Zn,:Zr
    ]
    nelements = length(biglist)
    strbiglist = string.(biglist)

    npoints = length(macrostrat.rocktype)
    subcats = collect(allkeys)
    deleteat!(subcats, findall(x->x==:cover,subcats))       # Do not compute cover


## --- Composition of exposed crust; start by getting samples with complete geochemical data
    # Number of samples measuring over 95% of rock components
    complete = Array{Float64}(undef, length(bulkidx), 1)
    @time for i in eachindex(bulkidx)
        if bulkidx[i] != 0
            for j in biglist
                complete[i] = nanadd(complete[i], bulk[j][bulkidx[i]])
            end
        else
            complete[i] = 0
        end
    end

    above95 = vec(complete .>= 95.0)
    totalgeochem = count(above95) / length(bulkidx) * 100
    if totalgeochem > 50
        @info "$(round(totalgeochem, sigdigits=3))% samples have geochemical data for > 95% total wt.%"
    else
        @warn "$(round(totalgeochem, sigdigits=3))% samples have geochemical data for > 95% total wt.%"
    end

    
## --- Compute average major element geochemistry for major rock types
    # Reduce analyzed data to only data with > 95% total wt. analyzed
    bulkidx95 = bulkidx[above95]

    # Get major elements, avoid hardcoding
    majorelem = collect(keys(major_elements(bulk, bulk_cats.sed, trues(count(bulk_cats.sed)))))

    # Preallocate
    majorcomp = Dict(
        :sed => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem)))),
        :met => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem)))),
        :ign => Dict(
            :fel => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem)))),     # Felsic
            :int => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem)))),     # Intermediate
            :maf => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem)))),     # Mafic
            :all => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem))))      # All igneous
        )
    )

    # Define igneous rock compositions by silica (from Keller and Schoene, 2012)
    ignsilica = (
        fel = (62, 74),      # Felsic (low exclusive, high inclusive)
        int = (51, 62),      # Intermediate
        maf = (43, 51),      # Mafic
        all = (0, 100)       # All igneous
    )

    # Get silica data
    silicadata = Array{Float64}(undef, length(bulkidx95[macro_cats.ign[above95]]), 1)
    for i in eachindex(silicadata)
        silicadata[i] = ifelse(bulkidx95[macro_cats.ign[above95]][i] != 0, 
            bulk.SiO2[bulkidx95[macro_cats.ign[above95]][i]], NaN
        )
    end

    # Compute composition of exposed crust for each rock subtype!
    for i in keys(majorcomp)
        # Igneous rocks separated by silica content
        if i==:ign
            for j in keys(majorcomp[i])
                # Get silica thresholds
                if j==:all
                    t = trues(length(silicadata))
                else
                    bound = ignsilica[j]
                    t = @. bound[1] < silicadata <= bound[2]
                end

                # Get data for each element
                for p in keys(majorcomp[i][j])
                    data = Array{Float64}(undef, length(bulkidx95[macro_cats[i][above95]]), 1)
                    for k in eachindex(data)
                        data[k] = ifelse(bulkidx95[macro_cats[i][above95]][k] != 0, 
                            bulk[p][bulkidx95[macro_cats[i][above95]][k]], NaN
                        )
                    end

                    # Reduce data to appropriate silica content, put in dictionary
                    data = data[t]
                    majorcomp[i][j][p] = nanmean(data) ± nanstd(data)
                end
            end

        # Other rocks, not differentiated
        else
            for j in keys(majorcomp[i])
                # Get data for the current element
                data = Array{Float64}(undef, length(bulkidx95[macro_cats[i][above95]]), 1)
                for k in eachindex(data)
                    data[k] = ifelse(bulkidx95[macro_cats[i][above95]][k] != 0, 
                        bulk[j][bulkidx95[macro_cats[i][above95]][k]], NaN
                    )
                end

                # Put data in dictionary
                majorcomp[i][j] = nanmean(data) ± nanstd(data)
            end
        end
    end


## --- Terminal printouts
    # Samples available for analysis for each major rock type, compared to crustal abundance
    sedtotal = count(macro_cats.sed)
    mettotal = count(macro_cats.met)
    igntotal = count(macro_cats.ign)

    sedfrac = sedtotal / total_known * 100
    metfrac = mettotal / total_known * 100
    ignfrac = igntotal / total_known * 100

    sedmeasure = count(macro_cats.sed .& above95) / sedtotal * 100
    metmeasure = count(macro_cats.met .& above95) / mettotal * 100
    ignmeasure = count(macro_cats.ign .& above95) / igntotal * 100

    @info """
    Samples with >95% measured geochemistry vs. crustal abundance:
    sed: $(round(sedmeasure, sigdigits=3))% of samples, $(round(sedfrac, sigdigits=3))% of rocks
    met: $(round(metmeasure, sigdigits=3))% of samples, $(round(metfrac, sigdigits=3))% of rocks
    ign: $(round(ignmeasure, sigdigits=3))% of samples, $(round(ignfrac, sigdigits=3))% of rocks
    """

    # Calculate crustal abundance of felsic / intermediate / mafic igneous rocks
    allsilicadata = Array{Float64}(undef, length(bulkidx[macro_cats.ign]), 1)
    for i in eachindex(allsilicadata)
        allsilicadata[i] = ifelse(bulkidx[macro_cats.ign][i] != 0, 
            bulk.SiO2[bulkidx[macro_cats.ign][i]], NaN
        )
    end

    ufel = count(@. ignsilica.fel[2] < allsilicadata)                      # Ultra-silicic
    fel = count(@. ignsilica.fel[1] < allsilicadata <= ignsilica.fel[2])   # Felsic
    int = count(@. ignsilica.int[1] < allsilicadata <= ignsilica.int[2])   # Intermediate
    maf = count(@. ignsilica.maf[1] < allsilicadata <= ignsilica.maf[2])   # Mafic
    umaf = count(@. ignsilica.maf[1] > allsilicadata)                      # Ultra-mafic

    allign = ufel + fel + int + maf + umaf      # All non-NaN igneous

    ufelfrac = ufel / allign * 100              # Percent abundances
    felfrac = fel / allign * 100
    intfrac = int / allign * 100
    maffrac = maf / allign * 100
    umaffrac = umaf / allign * 100

    # Samples available for analysis for igneous rock types, compared to crustal abundance
    ufelmeasure = count(@. ignsilica.fel[2] < silicadata)                      # Ultra-silicic
    felmeasure = count(@. ignsilica.fel[1] < silicadata <= ignsilica.fel[2])   # Felsic
    intmeasure = count(@. ignsilica.int[1] < silicadata <= ignsilica.int[2])   # Intermediate
    mafmeasure = count(@. ignsilica.maf[1] < silicadata <= ignsilica.maf[2])   # Mafic
    umafmeasure = count(@. ignsilica.maf[1] > silicadata)                      # Ultra-mafic

    allignmeasure = ufelmeasure + felmeasure + intmeasure + mafmeasure + umafmeasure

    ufelfracmeasure = ufelmeasure / allignmeasure * 100       # Percent abundances
    felfracmeasure = felmeasure / allignmeasure * 100
    intfracmeasure = intmeasure / allignmeasure * 100
    maffracmeasure = mafmeasure / allignmeasure * 100
    umaffracmeasure = umafmeasure / allignmeasure * 100

    @info """
    Igneous rocks by silica content with >95% measured geochemistry vs. crustal abundance:
    > 45%: $(round(umaffracmeasure, sigdigits=3))% of samples, $(round(umaffrac, sigdigits=3))% of rocks
    43-51%: $(round(maffracmeasure, sigdigits=3))% of samples, $(round(maffrac, sigdigits=3))% of rocks
    51-62%: $(round(intfracmeasure, sigdigits=3))% of samples, $(round(intfrac, sigdigits=3))% of rocks
    62-74%: $(round(felfracmeasure, sigdigits=3))% of samples, $(round(felfrac, sigdigits=3))% of rocks
    < 74%: $(round(ufelfracmeasure, sigdigits=3))% of samples, $(round(ufelfrac, sigdigits=3))% of rocks
    """


## --- Export crustal composition results
    # Preallocate
    bigmatrix = Array{Float64}(undef, length(keys(majorcomp[:sed])) + 1, 
        length(keys(majorcomp)) + length(keys(majorcomp[:ign])) - 1
    )
    exportdata = Array{Union{Float64, String}}(undef, size(bigmatrix).+1)

    # Define column and row names
    header = collect(keys(majorcomp))
    deleteat!(header, findall(x->x==:ign,header))
    push!(header, collect(keys(majorcomp[:ign]))...)

    rows = vcat(collect(keys(majorcomp[:sed])), :total)

    # Put data into matrix
    for i in eachindex(header)
        if header[i]==:sed || header[i]==:met
            for j in eachindex(rows[1:end-1])
                # By element
                bigmatrix[j, i] = majorcomp[header[i]][rows[j]].val

                # Total
                bigmatrix[end, i] = nansum(bigmatrix[1:end-1, i])
            end
        else
            for j in eachindex(rows[1:end-1])
                # By element
                bigmatrix[j, i] = majorcomp[:ign][header[i]][rows[j]].val

                # Total
                bigmatrix[end, i] = nansum(bigmatrix[1:end-1, i])
            end
        end
    end

    # Add column and row names
    exportdata[1,1] = ""
    exportdata[2:end, 2:end] = bigmatrix
    exportdata[1, 2:end] = string.(reshape(header, 1, length(header)))
    exportdata[2:end, 1] = string.(rows)

    # Write file
    writedlm("output/exposedcrust.tsv", exportdata)


# ## --- Create HDF5 file to store results for eroded material
#     fid = h5open("output/erodedmaterial.h5", "w")

#     # Metadata
#     write(fid, "element_names", strbiglist)                           # Names of analyzed elements
#     write(fid, "npoints", npoints)                                    # Total macrostrat samples
#     nsample = create_dataset(fid, "nsamples_byelem", Int, (nelements,))   # Non-NaN samples for each element

#     # Bulk rock global flux
#     bulkrockflux = create_group(fid, "bulkrockflux")
#     write(bulkrockflux, "rocktypes", string.(subcats))
#     bulkflux_val = create_dataset(bulkrockflux, "val", Float64, (length(subcats),))
#     bulkflux_std = create_dataset(bulkrockflux, "std", Float64, (length(subcats),))

#     # For each element in biglist, wt.% and flux by rock subtype, and total global flux
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
#     for i in eachindex(biglist)
#         # Calculate wt.%, flux, and global flux of each element
#         wt, flux, global_flux, n = flux_source(bulk[biglist[i]], bulkidx, erosion, macro_cats, 
#             crustal_area, elem=strbiglist[i]
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
    elem_head = create_dataset(element_flux, "header", string.(biglist))

    # Use rock type wt.% averages
    for i in eachindex(biglist)
        # Compute average wt.% for element i by rock type
        avgwt = Dict(zip(subcats, fill(NaN ± NaN, length(subcats))))
        for j in collect(keys(avgwt))
            avgwt[j] = nanmean(bulk[biglist[i]][bulk_cats[j]]) ± nanstd(bulk[biglist[i]][bulk_cats[j]])
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
