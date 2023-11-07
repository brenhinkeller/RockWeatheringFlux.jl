# Pared down version of ScreenBulk.jl

## --- Load and parse data file
    bulk = matread("data/bulk.mat")["bulk"];
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk));


## --- Load and parse metadata
    bulktext = matread("data/bulktext.mat")["bulktext"]

    # Type stabilize sample metadata
    elements = unique([String.(bulktext["elements"]); ["Units", "Methods"]])
    for n in elements
        bulktext[n] = String.(bulktext[n])
    end

    # Densify sparse arrays, correct for zero-indexing and Float64 type
    for k in keys(bulktext["unit"])
        bulktext["unit"][k] = collect(bulktext["unit"][k])
        bulktext["unit"][k] .+= 1.0
        bulktext["unit"][k] = Int.(bulktext["unit"][k])
    end
    for k in keys(bulktext["index"])
        bulktext["index"][k] = collect(bulktext["index"][k])
        bulktext["index"][k] .+= 1.0
        bulktext["index"][k] = Int.(bulktext["index"][k])
    end

    # Information by element
    bulktext_units = (
        units = bulktext["Units"],
        index = NamedTuple{Tuple(Symbol.(collect(keys(bulktext["unit"]))))}(
            [bulktext["unit"][i] for i in keys(bulktext["unit"])]
        )
    )
    bulktext = (
        elements = NamedTuple{Tuple(Symbol.(bulktext["elements"]))}(
            [bulktext[i] for i in bulktext["elements"]]
        ),
        index = NamedTuple{Tuple(Symbol.(collect(keys((bulktext["index"])))))}(
            [bulktext["index"][i] for i in keys(bulktext["index"])]
        )
    )
    

## --- Or save yourself the time and load from a file
    fid = h5open("output/bulk_unrestricted_types.h5", "r")
    header = read(fid["bulk_cats_head"])
    data = read(fid["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    
    # We want major types to include minor types
    # Remember there are no minor metamorphic rocks anymore! Goodbye! I'll miss you :(
    typelist, minorsed, minorvolc, minorplut, minorign = get_rock_class();
    for type in minorsed
        bulk_cats.sed .|= bulk_cats[type]
    end
    for type in minorvolc
        bulk_cats.volc .|= bulk_cats[type]
    end
    for type in minorplut
        bulk_cats.plut .|= bulk_cats[type]
    end
    for type in minorign
        bulk_cats.ign .|= bulk_cats[type]
    end


## --- Define and organize elements to convert units
    # Get major and minor elements, but remove :Volatiles from before it causes problems
    majors, minors = get_elements()
    allelements = [majors; minors]
    allkeys = [allelements; [:Latitude, :Longitude, :Loc_Prec, :Age, :Age_Max, :Age_Min]]
    deleteat!(allelements, findall(x->x==:Volatiles,allelements))

    # LOI, CaCO3, H2O, and CO2 don't get reported in the output, but we'll want them with 
    # the correct units for conversions 
    source_elements = [allelements; [:Loi, :CaCO3, :H2O, :CO2]]

    # Get all units present in bulktext
    presentunits = join(collect(keys(bulktext_units.index)), " ")

    # Define densities for relevant major element oxides (g/cm³)
    density = (Al2O3=3.95, K2O=2.35, MgO=3.58, Na2O=2.27, P2O5=2.39, SiO2=2.65, TiO2=4.23, 
        CaO=3.34, FeOT=5.74
    )


## --- Convert all units to wt.%
    # This method is resistant to changes in element order
    for i in source_elements
        # Check unit is present, and convert indices to a list of units
        unit_i = string(i) * "_Unit"
        @assert contains(presentunits, unit_i) "$unit_i not present"
        bulkunits = bulktext_units.units[bulktext_units.index[Symbol(unit_i)]]

        # Convert units of all sample data
        standardize_units!(bulk[i], bulkunits, density, elem=i)
    end


## --- Fix some missing data and weird conversions
    # Pre-compute some conversion factors
    CaCO3_to_CaO = (40.078+15.999)/(40.078+15.999*3+12.01)
    CaCO3_to_CO2 = (15.999*2+12.01)/(40.078+15.999*3+12.01)

    CaO_to_CO2 = (15.999*2+12.01)/(40.078+15.999)
    CaO_to_H2O = 2*(2*1.00784+15.999)/(40.078+15.999)
    CaO_to_SO3 = (32.065+3*15.999)/(40.078+15.999)

    # CaCO₃ → CaO + CO₂, and get CO₂ from CaO
    # Replace the existing CaO value with the computed CaO value
    # Add the existing CO2 value to the computed CO2 values 
    @turbo for i in eachindex(bulk.CaCO3)
        bulk.CaO[i] = ifelse(bulk.CaCO3[i] > 0, CaCO3_to_CaO * bulk.CaCO3[i], bulk.CaO[i])
        bulk.CO2[i] = ifelse(bulk.CaCO3[i] > 0, CaCO3_to_CO2 * bulk.CaCO3[i] + bulk.CO2[i],
            bulk.CO2[i]
        )
    end

    # Carbonates, siliciclasts, and shales: Calculate CO₂ from CaO
    # There are a lot of limestones with only CaO :(
    target = bulk_cats.carb .| bulk_cats.siliciclast .| bulk_cats.shale;
    @turbo for i in eachindex(bulk.CaO)
        bulk.CO2[i] += ifelse(target[i], bulk.CaO[i]*CaO_to_CO2, NaN)
    end

    # GYPSUM ONLY: CaO → H₂O, SO₄ (technically SO₃ because already have an oxygen in CaO)
    SO3 = Array{Float64}(undef, length(bulk.SiO2), 1)
    isgypsum = bulktext.elements.Rock_Name[bulktext.index.Rock_Name] .== "GYPSUM"
    @turbo for i in eachindex(bulk.H2O)
        bulk.H2O[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_H2O, NaN)
        SO3[i] = ifelse(isgypsum[i], bulk.CaO[i]*CaO_to_SO3, 0)
    end


## --- Create new volatiles measurement
    # Initialize as the larger of (CO₂ + H₂O) OR LOI, and add calculated SO₃
    volatiles = Array{Float64}(undef, length(bulk.SiO2), 1)

    zeronan!(bulk.CO2)
    zeronan!(bulk.H2O)
    @turbo volatiles .= bulk.CO2 .+ bulk.H2O

    @turbo for i in eachindex(volatiles)
        volatiles[i] = ifelse(bulk.Loi[i] > volatiles[i], bulk.Loi[i], volatiles[i])
    end

    @turbo volatiles .+= SO3


## --- Compute total wt.% analyzed for all samples
    # Exclude rocks below the shelf break
    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, bulk.Latitude, bulk.Longitude)
    abovesea = elev .> -140;    # Shelf break at 140 m

    # Compute bulk analyzed weight for rocks above sea level
    bulkweight = Array{Float64}(undef, length(bulk.SiO2), 1)

    @inbounds for i in eachindex(bulkweight)
        if abovesea[i]
            bulkweight[i] = nansum([bulk[j][i] for j in allelements]) + volatiles[i]
        else
            bulkweight[i] = NaN
        end
    end


## --- Correct for assumed unmeasured volatile loss in sedimentary rocks
    # If the sample is a sedimentary rock with a total analyzed wt.% below 84%, assume 
    # the "missing" data is volatiles that were not included in the analysis
    additional = zeros(length(bulkweight))
    for i in eachindex(bulkweight)
        additional[i] = ifelse(bulk_cats.sed[i] && bulkweight[i] < 100, 
            (100 - bulkweight[i]), 0.0)
    end


## --- End of File