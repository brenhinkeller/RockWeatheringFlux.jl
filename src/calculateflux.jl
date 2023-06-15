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

    sedtypes = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal)
    igntypes = (:volc, :plut)
    mettypes = (:metased, :metaign)
    alltypes = (sedtypes..., igntypes..., mettypes...)

    # Exclude cover from all major and minor rock types
    macro_cats.sed .&= .! macro_cats.cover
    macro_cats.ign .&= .! macro_cats.cover
    macro_cats.met .&= .! macro_cats.cover
    for i in alltypes
        macro_cats[i] .&= .! macro_cats.cover
    end

    # Exclude metamorphic rocks from sed and igns. Class as metased and metaign
    macro_cats.metased .|= (macro_cats.sed .& macro_cats.met)
    macro_cats.metaign .|= (macro_cats.ign .& macro_cats.met)

    macro_cats.sed .&= .! macro_cats.met
    for i in sedtypes
        macro_cats[i] .&= .! macro_cats.met
    end

    macro_cats.ign .&= .! macro_cats.met
    for i in igntypes
        macro_cats[i] .&= .! macro_cats.met
    end

    # Definitional exclusions
    macro_cats.cryst .&= .! macro_cats.metased      # Cryst rocks cannot be metasedimentary
    macro_cats.sed .&= .! macro_cats.ign            # Sed / ign rocks are classified as ign
    for i in sedtypes
        macro_cats[i] .&= .! macro_cats.ign 
    end

    # Exclude subtypes from other subtypes
    # Relatively arbitrary exclusion of sed types from each other
    # Re-order the list to change how the exclusion works
    for i in 1:(length(sedtypes)-1)
        for j in sedtypes[i+1:end]
            macro_cats[sedtypes[i]] .&= .! macro_cats[j]
        end
    end

    # Volcanic and plutonic is downgraded to undifferentiated igneous
    macro_cats.volc .&= .! (macro_cats.volc .& macro_cats.plut)
    macro_cats.plut .&= .! (macro_cats.volc .& macro_cats.plut)

    # Metaigneous and metasedimentary is downgraded to undifferentiated metamorphic
    macro_cats.metased .&= .! (macro_cats.metased .& macro_cats.metaign)
    macro_cats.metaign .&= .! (macro_cats.metased .& macro_cats.metaign)

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
    for i in eachindex(allkeys)
        erosion[allkeys[i]] = nanmean(rock_ersn[macro_cats[i]])
    end
    erosion = NamedTuple{Tuple(allkeys)}(values(erosion))

    # Crustal area (m²), assume proportional distribution of rock types under cover
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    for i in eachindex(allkeys)
        crustal_area[allkeys[i]] = count(macro_cats[i]) / total_known * contl_area
    end
    crustal_area = NamedTuple{Tuple(allkeys)}(values(crustal_area))
    

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


## --- Create HDF5 file to store results
    fid = h5open("output/rwf_output4.h5", "w")

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

## --- Get data
    # Calculations by element / compound
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

    # Total bulk rock mass flux by rock type
    # TO DO: make this into a measurement type
    const crustal_density = 2750    # (kg/m³)
    for i in eachindex(subcats)
        bulkflux1 = erosion[i] * crustal_area[i] * crustal_density* 1e-6
        bulkflux_val[i] = bulkflux
        # bulkflux_val[i] = bulkflux.val
        # bulkflux_std[i] = bulkflux.err
    end

    close(fid)

## --- End of File
