## --- Setup
    # Packages
    using StatGeochem
    using DelimitedFiles
    using ProgressMeter
    using LoopVectorization
    using Measurements
    using HDF5
    using MAT
    
    # Local utilities
    include("Utilities.jl")


## --- Load pre-generated Macrostrat data
    @info "Loading Macrostrat data"

    # Load and match
    macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)
    @time macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    minorsed = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal)
    minorign = (:volc, :plut)
    minormet = (:metased, :metaign)
    minortypes = (minorsed..., minorign..., minormet...)

    # Exclude cover from all major and minor rock types
    macro_cats.sed .&= .! macro_cats.cover
    macro_cats.ign .&= .! macro_cats.cover
    macro_cats.met .&= .! macro_cats.cover
    for i in minortypes
        macro_cats[i] .&= .! macro_cats.cover
    end

    # Exclude metamorphic rocks from sed and igns. Class as metased and metaign
    macro_cats.metased .|= (macro_cats.sed .& macro_cats.met)
    macro_cats.metaign .|= (macro_cats.ign .& macro_cats.met)

    macro_cats.sed .&= .! macro_cats.met
    for i in minorsed
        macro_cats[i] .&= .! macro_cats.met
    end

    macro_cats.ign .&= .! macro_cats.met
    for i in minorign
        macro_cats[i] .&= .! macro_cats.met
    end

    # Exclude subtypes from other subtypes
    # Relatively arbitrary exclusion of sed types from each other
    # Re-order the list to change how the exclusion works
    for i in 1:(length(minorsed)-1)
        for j in minorsed[i+1:end]
            macro_cats[minorsed[i]] .&= .! macro_cats[j]
        end
    end

    # Volcanic and plutonic is downgraded to undifferentiated igneous
    macro_cats.volc .&= .! (macro_cats.volc .& macro_cats.plut)
    macro_cats.plut .&= .! (macro_cats.volc .& macro_cats.plut)

    # Metaigneous and metasedimentary is downgraded to undifferentiated metamorphic
    macro_cats.metased .&= .! (macro_cats.metased .& macro_cats.metaign)
    macro_cats.metaign .&= .! (macro_cats.metased .& macro_cats.metaign)

    # Sed / ign rocks are classified as ign
    macro_cats.sed .&= .! macro_cats.ign            
    for i in minorsed
        macro_cats[i] .&= .! macro_cats.ign 
    end

    # Make sure crystalline rocks are only igneous and non-metased metamorphic rocks
    macro_cats.cryst .= (macro_cats.ign .| macro_cats.met)
    macro_cats.cryst .&= .! (macro_cats.sed .| macro_cats.cover .| macro_cats.metased)

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
      not matched = $(count(not_matched))"

    if count(multi_matched) > 0
        @warn "$(count(multi_matched)) conflicting matches present
          sed and ign = $(count(macro_cats.sed .& macro_cats.ign))
          sed and met = $(count(macro_cats.sed .& macro_cats.met))
          ign and met = $(count(macro_cats.ign .& macro_cats.met))"
    end


## --- DEBUGGING: Check that all subtypes are included in the major type
    # # Get BitVectors of all minor types
    # minorsed_bv = falses(length(macro_cats.sed))
    # for i in minorsed
    #     minorsed_bv .|= macro_cats[i]
    # end

    # minorign_bv = falses(length(macro_cats.ign))
    # for i in minorign
    #     minorign_bv .|= macro_cats[i]
    # end

    # minormet_bv = falses(length(macro_cats.met))
    # for i in minormet
    #     minormet_bv .|= macro_cats[i]
    # end

    # # Make sure all minor types are included in the major type
    # for i in eachindex(macro_cats.sed)
    #     if minorsed_bv[i]
    #         @assert macro_cats.sed[i] "Sedimentary index $i"
    #     end
    # end
    # for i in eachindex(macro_cats.ign)
    #     if minorign_bv[i]
    #         @assert macro_cats.ign[i] "Igneous index $i"
    #     end
    # end
    # for i in eachindex(macro_cats.met)
    #     if minormet_bv[i]
    #         @assert macro_cats.met[i] "Metamorphic index $i"
    #     end
    # end
    

## --- Calculate erosion rate at each coordinate point of interest
	@info "Calculating slope and erosion at each point"
	
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("data/srtm15plus_maxslope.h5", "vars/slope")
    srtm15_sf = h5read("data/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point
    rockslope = avg_over_area(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, 
        srtm15_sf, halfwidth=7
    )

    # Calculate all erosion rates (mm/kyr)
    # TO DO: update this function with a better erosion estimate
    # TO DO: update this to be a measurement value (i.e. prop. uncertainty!)
    rock_ersn = emmkyr.(rockslope)


## --- Get erosion and continental area for each rock type
    # Preallocate
    allkeys = keys(macro_cats)
    allinitvals = fill(NaN, length(allkeys))

    erosion = Dict(zip(allkeys, allinitvals))
    crustal_area = Dict(zip(allkeys, allinitvals))

    # Erosion (m/Myr)
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
    
## --- DEBUGGING: Total area of minor rock types should be less than area of major types
    # minorsed_area = 0.0
    # for i in minorsed
    #     minorsed_area += crustal_area[i]
    # end
    # @assert minorsed_area < crustal_area.sed

    # minorign_area = 0.0
    # for i in minorign
    #     minorign_area += crustal_area[i]
    # end
    # @assert minorign_area < crustal_area.ign

    # minormet_area = 0.0
    # for i in minormet
    #     minormet_area += crustal_area[i]
    # end
    # @assert minormet_area < crustal_area.met


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
    ndata = length(biglist)
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

    abovethreshold = complete .>= 95.0
    totalgeochem = count(abovethreshold) / length(bulkidx) * 100
    if totalgeochem > 50
        @info "$(round(totalgeochem, sigdigits=3))% samples measure > 95% total wt.%"
    else
        @warn "$(round(totalgeochem, sigdigits=3))% samples measure > 95% total wt.%"
    end

    # Sed / met / ign breakdown of total rock geochemistry measurements
    sedtotal = count(macro_cats.sed)
    mettotal = count(macro_cats.met)
    igntotal = count(macro_cats.ign)

    sedmeasure = count(macro_cats.sed .& abovethreshold) / sedtotal * 100
    metmeasure = count(macro_cats.met .& abovethreshold) / mettotal * 100
    ignmeasure = count(macro_cats.ign .& abovethreshold) / igntotal * 100

    @info "Samples with >95% measured geochemistry:
      sed: $(round(sedmeasure, sigdigits=3))%
      met: $(round(metmeasure, sigdigits=3))%
      ign: $(round(ignmeasure, sigdigits=3))%
    "
    
## --- Compute average major element geochemistry for major rock types (sanity check)
    # Get major elements, avoid hardcoding
    majorelem = collect(keys(major_elements(bulk, bulk_cats.sed, trues(count(bulk_cats.sed)))))

    # Preallocate
    majorcomp = Dict(
        :sed => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem)))),
        :met => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem)))),
        :ign => Dict(zip(majorelem, fill(NaN±NaN, length(majorelem))))
    )

    for i in keys(majorcomp)
        for j in keys(majorcomp[i])
            # Get data for that element
            data = Array{Float64}(undef, length(bulkidx[macro_cats[i]]), 1)
            for k in eachindex(data)
                if bulkidx[macro_cats[i]][k] != 0
                    data[k] = bulk[j][bulkidx[macro_cats[i]][k]]
                else
                    data[k] = NaN
                end

                # data[k] = ifelse(bulkidx[macro_cats[i]][k] != 0, bulk[j][bulkidx[macro_cats[i]][k]], NaN)
            end

            # Put data in dictionary
            majorcomp[i][j] = nanmean(data) ± nanstd(data)
        end
    end


## --- Create HDF5 file to store results for composition of exposed crust
    fid = h5open("output/exposedcrust.h5", "w")

    close(fid)


## --- Create HDF5 file to store results for eroded material
    fid = h5open("output/erodedmaterial.h5", "w")

    # Metadata
    write(fid, "element_names", strbiglist)                           # Names of analyzed elements
    write(fid, "npoints", npoints)                                    # Total macrostrat samples
    nsample = create_dataset(fid, "nsamples_byelem", Int, (ndata,))   # Non-NaN samples for each element

    # Bulk rock global flux
    bulkrockflux = create_group(fid, "bulkrockflux")
    write(bulkrockflux, "rocktypes", string.(subcats))
    bulkflux_val = create_dataset(bulkrockflux, "val", Float64, (length(subcats),))
    bulkflux_std = create_dataset(bulkrockflux, "std", Float64, (length(subcats),))

    # For each element in biglist, wt.% and flux by rock subtype, and total global flux
    elementflux = create_group(fid, "elementflux")

    # Global flux
    totalflux = create_group(elementflux, "totalelemflux")
    totalflux_val = create_dataset(totalflux, "val", Float64, (ndata,))
    totalflux_std = create_dataset(totalflux, "std", Float64, (ndata,))

    # Separated by rock subtypes
    byrocktype = create_group(elementflux, "byrocktype")
    for i in subcats
        typegroup = create_group(byrocktype, string(i))

        create_dataset(typegroup, "wtpct_val", Float64, (ndata,))
        create_dataset(typegroup, "wtpct_std", Float64, (ndata,))

        create_dataset(typegroup, "flux_val", Float64, (ndata,))
        create_dataset(typegroup, "flux_std", Float64, (ndata,))
    end


## --- Calculate undifferentiated (bulk) flux
    const crustal_density = 2750    # (kg/m³)

    # Keys for erosion and crustal_area are in an arbitrary order. Correct indexing critical
    # TO DO: make this into a measurement type
    bulkflux = Dict(zip(subcats, fill(NaN, length(subcats))))
    for i in subcats
        bulkflux[i] = erosion[i] * crustal_area[i] * crustal_density * 1e-6
    end
    bulkflux = NamedTuple{Tuple(keys(bulkflux))}(values(bulkflux))

    # Save to file. Keys are again in arbitrary order
    for i in eachindex(subcats)
        sub_i = subcats[i]
        bulkflux_val[i] = bulkflux[sub_i]
        # bulkflux_val[i] = bulkflux[sub_i].val
        # bulkflux_std[i] = bulkflux[sub_i].err
    end


## --- DEBUGGING: Total flux of minor rock types should be less than flux of major types
    # minorsed_flux = 0.0
    # for i in minorsed
    #     minorsed_flux += bulkflux[i]
    # end
    # @assert minorsed_flux < bulkflux.sed

    # minorign_flux = 0.0
    # for i in minorign
    #     minorign_flux += bulkflux[i]
    # end
    # @assert minorign_flux < bulkflux.ign

    # minormet_flux = 0.0
    # for i in minormet
    #     minormet_flux += bulkflux[i]
    # end
    # @assert minormet_flux < bulkflux.met


## --- Calculate flux by element / element oxide
    # TO DO: I should make sure these aren't also scrambled...
    for i in eachindex(biglist)
        # Calculate wt.%, flux, and global flux of each element
        wt, flux, global_flux, n = flux_source(bulk[biglist[i]], bulkidx, erosion, macro_cats, 
            crustal_area, elem=strbiglist[i]
        )

        # Write data to file
        nsample[i] = n
        totalflux_val[i] = global_flux.val
        totalflux_std[i] = global_flux.err

        for j in eachindex(subcats)
            typegroup = string(subcats[j])

            byrocktype[typegroup]["flux_val"][i] = flux[subcats[j]].val
            byrocktype[typegroup]["flux_std"][i] = flux[subcats[j]].err
            
            byrocktype[typegroup]["wtpct_val"][i] = wt[subcats[j]].val
            byrocktype[typegroup]["wtpct_std"][i] = wt[subcats[j]].err
        end
    end

    close(fid)


## --- End of File
