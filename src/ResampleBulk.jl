#  Resample all rock types based on spatial weights (after Keller et al., 2015)
# 
# This will give me something to compare SampleMatch results to.

# May run into some problems with minor sedimentary types, and it's unclear if sorting 
# metamorphic rocks by metased / metaign is doing anything useful

## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5, MAT
    

## --- Load EarthChem data 
    # Filtered and normalized to 100%
    fid = h5open("output/bulk.h5", "r")

    # Bulk
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    # Rock type matches
    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])
    
    # Rock name matches
    # rocknames = read(fid["bulktypes"]["bulk_lookup_head"])
    # data = read(fid["bulktypes"]["bulk_lookup"])
    # data = @. data > 0
    # bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}([data[:,i] 
    #     for i in eachindex(rocknames)])
    
    # Metadata
    path = fid["bulktext"]["sampledata"]
    header = read(path["header"])
    index = read(path["index"])

    target = ["Rock_Name", "Type", "Material"]
    targetind = [findall(==(i), header)[1] for i in target]

    bulktext = NamedTuple{Tuple(Symbol.(target))}(
        [lowercase.(read(path["elements"][target[i]]))[index[:,targetind[i]]] 
            for i in eachindex(target)]
    )

    close(fid)


## --- Definitions we'll need for resampling
    # # Major types are inclusive of all minor types
    # minorsed, minorign, minormet = get_minor_types()
    # for type in minorsed
    #     bulk_cats.sed .|= bulk_cats[type]
    # end
    # for type in minorign
    #     bulk_cats.ign .|= bulk_cats[type]
    # end
    # for type in minormet
    #     bulk_cats.met .|= bulk_cats[type]
    # end

    # # Get the data we'll use in resampling. Exclude type and header row
    # bulkmatrix = float.(unelementify(bulk)[2:end , 1:end .!= end-3])
    # header = collect(keys(bulk))[1:end-4]

    # # Uncertanties are 0, except for SiO₂ which is 1.0, after Keller et al., 2015
    # uncert = zeros(size(bulkmatrix))
    # i = findfirst(x -> x==:SiO2, header)
    # uncert[:,i] .= 1.0


## --- Resample each rock type based on spatial weights
    # types = keys(bulk_cats)

    # # Everyone's favorite file format!
    # fid = h5open("output/resampled/resampled.h5", "w")
    # g_main = create_group(fid, "vars")
    # g = create_group(g_main, "data")

    # for t in types
    #     t==:cover && continue
    #     println("Starting $t")

    #     # Calculate spatial weights, keep ∼1/5 of the data in each resampling
    #     k = invweight_location(bulk.Latitude[bulk_cats[t]], bulk.Longitude[bulk_cats[t]])
    #     p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)

    #     # Get all the data we'll use in the resampling
    #     rockdata = bulkmatrix[bulk_cats[t][:],:]
    #     rockuncert = uncert[bulk_cats[t][:],:]

    #     # Final dataset size should be proportional to the starting dataset, just because
    #     # we have 270K igneous samples and 6 evaporite samples. Minimum 500 samples
    #     nrows = max(count(bulk_cats[t]) * 5, 500)

    #     sim = bsresample(rockdata, rockuncert, nrows, p)

    #     # Save the data to the file, or whatever
    #     g₀ = create_group(g, "$t")
    #     g₀["data"] = sim
    #     g₀["k"] = k
    # end

    # # May as well just resample the everything too, as a treat
    # k = invweight_location(bulk.Latitude, bulk.Longitude)
    # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    # sim = bsresample(bulkmatrix, uncert, 1_500_000, p)

    # # And save
    # g₀ = create_group(g, "bulk")
    # g₀["data"] = sim
    # g₀["k"] = k

    # # Save the header too, we'll want that
    # g_main["header"] = string.(header)

    # close(fid)


## --- Resample all rocks and known rock names
    # Get list of defined rock names and assign each name a number
    rocknames = get_rock_class(major=true, inclusive=false)
    rocknames = unique((rocknames.sed..., rocknames.met..., rocknames.ign...))
    kz = collect(1:length(rocknames))

    # Match rock names
    p = Progress(length(rocknames), desc="Matching EarthChem samples and rock names")
    bulk_lookup = NamedTuple{Tuple(Symbol.(rocknames))}(
        [falses(length(bulktext.Rock_Name)) for _ in eachindex(rocknames)]
    )
    for i in eachindex(rocknames)
        bulk_lookup[i] .= find_earthchem(rocknames[i], bulktext.Rock_Name, bulktext.Type, 
            bulktext.Material)
        next!(p)
    end

    # Find the maximum number of names matched to a single sample
    maxmatches = 0
    for i in eachindex(bulk.SiO2)
        allnames = get_type(bulk_lookup, i, all_keys=true)
        if allnames===nothing
            continue
        elseif length(allnames) > maxmatches
            maxmatches = length(allnames)
        end
    end

    # Store rock name data
    p = Progress(length(bulk.SiO2), desc="Storing rock name data")
    kz_data = zeros(length(bulk.SiO2), maxmatches)
    rocknames = Symbol.(rocknames)
    for i in eachindex(bulk.SiO2)
        allnames = get_type(bulk_lookup, i, all_keys=true)
        if allnames===nothing
            continue
        else
            # Get kz numbers for each matched name
            kzn = zeros(maxmatches)
            for j in eachindex(allnames)
                for k in eachindex(rocknames)
                    if allnames[j] == rocknames[k]
                        kzn[j] = k
                        break
                    end
                end
            end
            kz_data[i,:] .= kzn
        end
        next!(p)
    end

    # Get appropriately sized header
    n = string.(collect(1:maxmatches))
    kz_header = ["Kz_" * s for s in n]

    # Get uncertainties for each element from Keller et al., 2015 data
    plutonic = matread("data/volcanicplutonic/plutonic.mat")["plutonic"];
    errs = NamedTuple{Tuple(Symbol.(keys(plutonic["err"])))}(values(plutonic["err"]))

    major, minor = get_elements()
    elements = [major; minor]
    elem_uncert = Array{Float64}(undef, length(elements), 1)
    for i in eachindex(elem_uncert)
        eᵢ = elements[i]
        if haskey(errs, eᵢ)
            elem_uncert[i] = errs[eᵢ]
        else
            elem_uncert[i] = 0.02   # Assumed error for Volatiles and In
        end
    end

    # Calculate age uncertainty. If min / max bounds are unknown, set uncertainty at 5%
    age_uncert = Array{Float64}(undef, length(bulk.SiO2), 1)
    for i in eachindex(age_uncert)
        age_uncert[i] = (bulk.Age_Max[i] - bulk.Age_Min[i]) / 2
    end
    for i in eachindex(age_uncert)
        age_uncert[i] = ifelse(isnan(age_uncert[i]), bulk.Age[i] * 0.05, age_uncert[i])
    end

    # Get whole-dataset information
    header = [string.(elements); ["Lat", "Lon", "Age"]; kz_header]
    ndata = length(header)
    nsamples = length(bulk.SiO2)
    kz_n = length(kz_header)

    # Uncertainties
    uncert = zeros(nsamples, ndata)
    for i = 1:(ndata-kz_n-3)                    # Elements
        uncert[:,i] .= elem_uncert[i]
    end
    uncert[:,end-kz_n-1] .= bulk.Loc_Prec       # Latitude
    uncert[:,end-kz_n] .= bulk.Loc_Prec         # Longitude

    # Data
    data = Array{Float64}(undef, nsamples, ndata-kz_n)
    for i = 1:(ndata-kz_n-3)
        data[:,i] .= bulk[elements[i]]
    end
    data[:,end-kz_n-1] .= bulk.Latitude
    data[:,end-kz_n] .= bulk.Latitude
    data = hcat(data, kz_data)

    # Resample!
    nsims = Int(1e6)
    k = invweight_location(bulk.Latitude, bulk.Longitude)
    p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
    simout = bsresample(data, uncert, nsims, p)

    # Save to a file 
    fid = h5open("output/resampled/resampled_rocknames.h5", "w")
    g = create_group(fid, "vars")
        g["header"] = header
        g["rocknames_kz"] = string.(rocknames)
        g["k"] = k
        g["simout"] = simout
    close(fid)


## --- Test resampled behavior against VolcanicPlutonic
    # # Load data
    # plutonic = matread("data/volcanicplutonic/plutonic.mat")["plutonic"];
    # plutonic = NamedTuple{Tuple(Symbol.(keys(plutonic)))}(values(plutonic));
    # src = plutonic

    # # volcanic = matread("data/volcanicplutonic/volcanic.mat")["volcanic"];
    # # volcanic = NamedTuple{Tuple(Symbol.(keys(volcanic)))}(values(volcanic));
    # # src = volcanic

    # SiO2min, SiO2max = 40, 80
    # simitems = [:SiO2];
    # nsims = Int(1e7)

    # # Get resampling weights
    # # test = trues(length(src.Latitude))
    # test = @. !isnan(src.Latitude) & !isnan(src.Longitude) & (src.Elevation .> -100);
    # # k = src.k[test]
    # k = invweight_location(src.Latitude[test], src.Longitude[test])
    # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0);

    # # Data matrix
    # data=zeros(length(src.SiO2),length(simitems));
    # for i in eachindex(simitems)
    #     data[:,i] = src[simitems[i]];
    # end
    # data = data[test[:],:];

    # # Uncertainty matrix
    # uncertainty=zeros(length(src.SiO2),length(simitems));
    # for i in eachindex(simitems)
    #     uncertainty[:,i] .= src.CalcAbsErr[string(simitems[i])];
    # end
    # uncertainty = uncertainty[test[:],:];

    # # Run Monte Carlo simulation
    # simout = bsresample(data, uncertainty, nsims, p)

    # # Plot results    
    # c, n = bincounts(simout, SiO2min, SiO2max, 160)
    # n = float(n) ./ nansum(float(n) .* step(c))
    # h = plot(c, n, seriestype=:bar, framestyle=:box, color=:darkblue, linecolor=:darkblue,
    #     label="Plutonic (n = $(count(test)))", ylabel="Abundance", xlabel="SiO2 [wt.%]",
    #     ylims=(0, round(maximum(n), digits=2)+0.01), xlims=(SiO2min, SiO2max)
    # )
    # display(h)


## --- End of file