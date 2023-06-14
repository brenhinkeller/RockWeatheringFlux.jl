# This file will read bulk.mat and bulktext.mat and standardize all units as wt.%
# Data will be saved as a new file
    # Currently a new .mat file

# TO DO: 
    # save as an .h5 file. Requires re-writing other files so putting off for now.
    # Re-read in old files, and figure out what to do for parsing / correcting element 
        # oxide weights. For example, some samples may have different values for CaO and 
        # CaCO3, but Mg and MgO represent the same data.

## --- Set up
    # Packages
    using MAT
    using StatGeochem
    
    # Load and parse data files
    bulk = matread("data/bulk.mat")["bulk"]
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))

    bulktext = matread("data/bulktext.mat")["bulktext"]

    # Type stabilize sample metadata
    bulktext["elements"] = unique([String.(bulktext["elements"]); ["Units", "Methods"]])
    for n in bulktext["elements"]
        bulktext[n] = String.(bulktext[n])
    end
    bulktext = NamedTuple{Tuple(Symbol.(keys(bulktext)))}(values(bulktext))

    # Densify the sparse unit arrays, correct for zero-indexing and Float64 type
    for k in keys(bulktext.unit)
        bulktext.unit[k] = collect(bulktext.unit[k])
        bulktext.unit[k] .+= 1.
        bulktext.unit[k] = Int.(bulktext.unit[k])
    end

## --- Define and organize elements to correct
    majors = [:SiO2,:Al2O3,:Fe2O3T,:TiO2,:MgO,:CaO,:Na2O,:K2O,]
    minors = [:Ag,:As,:Au,:B,:Ba,:Be,:Bi,:C,:CaCO3,:Cd,:Ce,:Cl,:Co,:Cr2O3,:Cs,:Cu,:Dy,:Er,
        :Eu,:F,:Ga,:Gd,:Hf,:Hg,:Ho,:I,:In,:Ir,:La,:Li,:Lu,:MnO,:Mo,:Nb,:Nd,:NiO,:Os,
        :P2O5,:Pb,:Pd,:Pt,:Pr,:Re,:Rb,:Sb,:Sc,:Se,:S,:Sm,:Sn,:Sr,:Ta,:Tb,:Te,:Th,:Tl,:Tm,
        :U,:V,:W,:Y,:Yb,:Zn,:Zr
    ]

    # Collect major and minor elements together
    allconstit = vcat(majors, minors)
    ndata = length(allconstit)
    
    # Check that there are units for all elements of interest
    unitslist = string.(allconstit) .* "_Unit"
    possiblekeys = join(collect(keys(bulktext.unit)), " ")
    missingkeys = .! contains.(possiblekeys, unitslist)
    @assert count(missingkeys) == 0 "$(unitslist[missingkeys]) not present"

    # Define densities for relevant major element oxides
    density = (Al2O3=3.95,K2O=2.35,MgO=3.58,Na2O=2.27,P2O5=2.39,SiO2=2.65,TiO2=4.23)


## --- Convert all units to wt.%
    # This is slightly less efficient, but is resistant to changes in element order
    for i in allconstit
        # Convert list of unit indices to a list of units
        unit_i = string(i) * "_Unit"
        bulkunits = bulktext.Units[bulktext.unit[unit_i]]

        # Convert all sample data
        standardize_units!(bulk[i], bulkunits, elem=string(i))
    end

## --- Write corrected data to a new file

## --- End of File