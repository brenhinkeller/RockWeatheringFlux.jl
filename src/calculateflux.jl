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

    # Exclude cover from all rock types
    macro_cats.sed .&= .! macro_cats.cover
    macro_cats.ign .&= .! macro_cats.cover
    macro_cats.met .&= .! macro_cats.cover

    # Exclude metamorphic rocks from sed and igns. Class as metased and metaign
    macro_cats.metased .|= (macro_cats.sed .& macro_cats.met)
    macro_cats.metaign .|= (macro_cats.ign .& macro_cats.met)

    macro_cats.sed .&= .! macro_cats.met
    macro_cats.ign .&= .! macro_cats.met

    # Definitional exclusions
    macro_cats.cryst .&= .! macro_cats.metased      # Cryst rocks cannot be metasedimentary
    macro_cats.sed .&= .! macro_cats.ign            # Sed / ign rocks are classified as ign

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
    allkeys = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal, :sed,
        :volc, :plut, :ign,
        :metased, :metaign, :met,
        :cryst
    )
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

    # Bulk data
    bulk_raw = matopen("data/bulk.mat")
    bulk_dict = read(bulk_raw, "bulk")
    close(bulk_raw)
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

    # Match codes to rock types
    bulk_cats = match_earthchem(bulk.Type)

    # Load indices of matched samples from samplematch.jl
    bulkidx = Int.(vec(readdlm("output/matched_bulkidx.tsv")))


## --- Set up arrays and files to get data for Everything
    # Every element in EarthChem
    biglist = sort([:Ag,:Al2O3,:As,:Au,:B,:Ba,:Be,:Bi,:C,:CaCO3,:Cd,:Ce,:Cl,:Co,:Cr2O3,:Cs,
        :Cu,:Dy,:Er,:Eu,:F,:Fe2O3T,:Ga,:Gd,:Hf,:Hg,:Ho,:I,:In,:Ir,:K2O,:La,:Li,:Lu,
        :MgO,:MnO,:Mo,:Na2O,:Nb,:Nd,:NiO,:Os,:P2O5,:Pb,:Pd,:Pt,:Pr,:Re,:Rb,:Sb,:Sc,:Se,
        :S,:SiO2,:Sm,:Sn,:Sr,:Ta,:Tb,:Te,:Th,:TiO2,:Tl,:Tm,:U,:V,:W,:Y,:Yb,:Zn,:Zr
    ])
    ndata = length(biglist)

    # Density in g/cm³, only for the elements stored as vol.%
    biglistdensity = (Al2O3=3.95,K2O=2.35,MgO=3.58,Na2O=2.27,P2O5=2.39,SiO2=2.65,TiO2=4.23)


## --- Get units of each measurement
    # Metadata
    bulktext_raw = matopen("data/bulktext.mat")
    bulktext_dict = read(bulktext_raw, "bulktext")
    close(bulktext_raw)
    bulktext = NamedTuple{Tuple(Symbol.(keys(bulktext_dict)))}(values(bulktext_dict))

    #=
    Things to consider...
        not every measurement has a unit
        units may not be the same for each element...

    =#
    # Add a progress bar to this dear god
    # Also maybe only do this for the data I'm planning on analyzing?
    possibleunit = collect(keys(bulktext.unit))
    totalkeys = ""
    for i in possibleunit
        totalkeys = totalkeys * " " * i
    end

    p = Progress(ndata, desc="Finding units for each EarthChem sample...")
    unitsused = falses(length(bulktext.Units))
    unitsbyelem = map(copy, fill(unitsused, ndata))

    for i in eachindex(biglist)
        # Get unit array and correct for 0-index
        unitname = string(biglist[i]) * "_Unit"
        
        if !contains(totalkeys, unitname)
            @warn "No units for $unitname"
            continue
        end

        unitarray = Int.(bulktext.unit[unitname] .+ 1)
        unitsused = zeros(length(bulktext.Units))
        missinginfo = falses(length(bulk[biglist[i]]))

        # Go through each data point for the element i in bulk
        # Get the unit that it's in
            # If there's a unit present, convert to wt.%
            # If not.... 
                # Need to infer units. Probably what we have to do is create a counter system for each unit
                # Each time there's a unit, add 1 to the counter for that unit
                # After going through all the data that has units, go through the data that doesn't have units
                # And assume that the units are whatever unit was used most frequently

        for j in eachindex(bulk[biglist[i]])
            if isnan(bulk[biglist[i]][j])
                continue
            else
                unitindex = unitarray[j]

                if unitindex != 1
                # If units are present (e.g. not 1) convert the units
                    unitsused[unitindex] += 1
                    unit = bulktext.Units[unitindex]

                    if unit=="WT%"
                        # No need to convert
                        continue

                    elseif unit=="DPM/G"
                        # Disintegrations per minute, not a measure of concentration
                        bulk[biglist[i]][j] = NaN

                    elseif unit=="PPM" || unit=="MICROGRAM PER GRAM" || unit=="MILLIGRAM PER KILOGRAM"
                        # Assume ppm by mass. Conversion is the same for all units
                        bulk[biglist[i]][j] /= 10000

                    elseif unit=="VOL%"
                        # Convert volume to mass using density
                        # TO DO: make this break-proof in case I don't have density measurement
                        bulk[biglist[i]][j] *= biglistdensity[i]

                    elseif unit=="RATIO"
                        # L + ratio + the ratio of what to what??

                    else
                        @warn "Unrecognized unit $unit in $(biglist[i])."
                        missinginfo[j] = true
                    end
                else
                # If units are not present, note the index of that unit
                    missinginfo[j] = true
                end
            end
        end

        # Loop through the list of data without unit information
        next!(p)
    end

    # Currently I just want this to give me a list of all the units that get used
    for i in eachindex(unitsbyelem)
        unitlist = bulktext.Units[unitsbyelem[i]]
        allunit = ""
        for j in unitlist
            allunit = allunit * " " * j
        end
        println("$(biglist[i]) uses $allunit")
    end


    # ndata = length(biglist)
    # strbiglist = string.(biglist)

    # # Define which elements are wt.%, ppm, or ppb
    # biglistunits = (Ag=:ppm,Al2O3=:wt,As=:ppm,Au=:ppb,B=:ppm,Ba=:ppm,Be=:ppm,Bi=:ppm,
    #     C=:ppm,CaCO3=:wt,Cd=:ppm,Ce=:ppm,Cl=:ppm,CoO=:wt,Cr2O3=:wt,Cs=:ppm,Cu=:ppm,
    #     Dy=:ppm,Er=:ppm,Eu=:ppm,F=:ppm,Fe2O3T=:wt,Ga=:ppm,Gd=:ppm,Hf=:ppm,Hg=:ppm,
    #     Ho=:ppm,I=:ppm,In=:ppm,Ir=:ppm,K2O=:wt,La=:ppm,Li=:ppm,Lu=:ppm,MgO=:wt,
    #     MnO=:wt,Mo=:ppm,Na2O=:wt,Nb=:ppm,Nd=:ppm,NiO=:wt,Os=:ppb,P2O5=:wt,Pb=:ppm,
    #     Pd=:ppm,Pt=:ppm,Pr=:ppm,Re=:ppb,Rb=:ppm,Sb=:ppm,Sc=:ppm,Se=:ppm,SO2=:wt,
    #     SiO2=:wt,Sm=:ppm,Sn=:ppm,Sr=:ppm,Ta=:ppm,Tb=:ppm,Te=:ppm,Th=:ppm,TiO2=:wt,
    #     Tl=:ppm,Tm=:ppm,U=:ppm,V=:ppm,W=:ppm,Y=:ppm,Yb=:ppm,Zn=:ppm,Zr=:ppm
    # )

    # Alternative way to store this data... probably more human and less machine readable
    # Elements in EarthChem, by unit
    pct = [:Al2O3, :CaCO3, :CoO, :Cr2O3, :Fe2O3T, :K2O, :MgO, :MnO, :Na2O, :NiO, :P2O5, 
        :SO2, :SiO2, :TiO2
    ]
    ppm = [:Ag, :As, :B, :Ba, :Be, :Bi, :C, :Cd, :Ce, :Cl, :Cs, :Cu, :Dy, :Er, :Eu, :F, 
        :Ga, :Gd, :Hf, :Hg, :Ho, :I, :In,  :La, :Li, :Lu, :Mo, :Nb, :Nd, :Pb, :Pr, :Rb, 
        :Sb, :Sc, :Se, :Sm, :Sn, :Sr, :Ta, :Tb, :Te, :Th, :Tl, :Tm, :U, :V, :W, :Y, :Yb, 
        :Zn, :Zr
    ]
    ppb = [:Pd, :Os, :Ir, :Pt, :Au, :Re]
    biglist = vcat(pct, ppm, ppb)
    ndata = length(biglist)
    strbiglist = string.(biglist)

    # Create file to store data
    npoints = length(macrostrat.rocktype)
    subcats = [:siliciclast, :shale, :carb, :chert, :evaporite, :coal, :sed, :volc, :plut, 
        :ign, :metased, :metaign, :met, :cryst
    ]

    # Create HDF5 file to store results
    fid = h5open("output/rwf_output3.h5", "w")

    # Metadata
    write(fid, "element_names", strbiglist)                     # Names of analyzed elements
    write(fid, "npoints", npoints)                              # Total macrostrat samples
    nsample = create_dataset(fid, "element_n", Int, (ndata,))   # Non-NaN samples for each element
    element_names = create_dataset(fid, "element_names", String, ((),))

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
    @info "Calculating results"

    # Calculations by element / compound
    for i in eachindex(biglist)
        # Make sure units are in wt.%
        # unit = biglistunits[biglist[i]]
        # if unit==:ppm
        #     data = bulk[biglist[i]] ./ 10000
        # elseif unit==:ppb
        #     data = bulk[biglist[i]] ./ 10000000
        # elseif unit==:wt
        #     data = bulk[biglist[i]]
        # else
        #     @warn "$(biglist[i]) Unit not recognized"
        # end



        # Calculate wt.%, flux, and global flux of each element
        wt, flux, global_flux, n = flux_source(data, bulkidx, erosion, macro_cats, 
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
        netflux = erosion[i] * crustal_area[i] * crustal_density* 1e-6
        bulkflux_val[i] = netflux
        # bulkflux_val[i] = netflux.val
        # bulkflux_std[i] = netflux.err
    end

    close(fid)


## --- End of File
