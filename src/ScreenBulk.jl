# This file will read bulk.mat and bulktext.mat and standardize all units as wt.%
# Data will be saved as a new file
    # Currently a new .mat file

# TO DO: 
    # save as an .h5 file. Requires re-writing other files so putting off for now.
    # Re-read in old files, and figure out what to do for parsing / correcting element 
        # oxide weights. For example, some samples may have different values for CaO and 
        # CaCO3, but Mg and MgO represent the same data.
    # Consider making multiple utilities files...
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

## --- Set up
    # Packages
    using MAT
    using StatGeochem
    using LoopVectorization

    # Local utilities
    include("Utilities.jl")
    
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
        bulktext.unit[k] .+= 1.0
        bulktext.unit[k] = Int.(bulktext.unit[k])
    end

## --- Define and organize elements to correct
    majors = [:SiO2,:Al2O3,:FeOT,:TiO2,:MgO,:CaO,:Na2O,:K2O,]
    minors = [:Ag,:As,:Au,:B,:Ba,:Be,:Bi,:C,:CaCO3,:Cd,:Ce,:Cl,:Co,:Cr2O3,:Cs,:Cu,:Dy,:Er,
        :Eu,:F,:Ga,:Gd,:Hf,:Hg,:Ho,:I,:In,:Ir,:La,:Li,:Lu,:MnO,:Mo,:Nb,:Nd,:NiO,:Os,
        :P2O5,:Pb,:Pd,:Pt,:Pr,:Re,:Rb,:Sb,:Sc,:Se,:S,:Sm,:Sn,:Sr,:Ta,:Tb,:Te,:Th,:Tl,:Tm,
        :U,:V,:W,:Y,:Yb,:Zn,:Zr
    ]

    # Collect major and minor elements together
    allconstit = [majors; minors]
    ndata = length(allconstit)

    # Get all units present in bulktext
    presentunits = join(collect(keys(bulktext.unit)), " ")

    # Define densities for relevant major element oxides (g/cm³)
    density = (Al2O3=3.95, K2O=2.35, MgO=3.58, Na2O=2.27, P2O5=2.39, SiO2=2.65, TiO2=4.23, 
        CaO=3.34, FeOT=5.74
    )


## --- Convert all units to wt.%
    # This method is resistant to changes in element order
    for i in allconstit
        # Check unit is present, and convert indices to a list of units
        unit_i = string(i) * "_Unit"
        @assert contains(presentunits, unit_i) "$unit_i not present"
        bulkunits = bulktext.Units[bulktext.unit[unit_i]]

        # Convert units of all sample data
        standardize_units!(bulk[i], bulkunits, density, elem=i)
    end


## --- Restrict data to 84-104 wt.% analyzed
    bulkweight = zeros(length(bulk.SiO2));
    @time @turbo for i in eachindex(allconstit)
        for j in eachindex(bulkweight)
            bulkweight[j] = nanadd(bulkweight[j], bulk[allconstit[i]][j])
        end
    end

    # Restrict data to within bounds and write to file
    bounds = (84, 104)
    t = @. bounds[1] < bulkweight < bounds[2]

    nsamples = round(count(t)/length(t)*100, digits=2)
    @info "Saving $nsamples% samples between $(bounds[1])% and $(bounds[2])% analyzed wt."

    bulknew = Dict()
    for i in allconstit
        bulknew[string(i)] = bulk[i][t]
    end
    matwrite("data/bulknew.mat", bulknew)


## --- End of File