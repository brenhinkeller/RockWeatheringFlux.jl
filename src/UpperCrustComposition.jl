## --- Set up
    # Packages
    using StatGeochem
    using DelimitedFiles
    using Measurements
    using MAT

    # Local utilities
    include("Utilities.jl")
    include("NaNMeasurements.jl")

## --- Load data
    # EarthChem
    bulk = matread("data/bulk_newunits.mat")
    bulk = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk))
    bulk_cats = match_earthchem(bulk.Type)

    # Matched EarthChem samples
    bulkidx = Int.(vec(readdlm("output/matched_bulkidx.tsv")))

    # Macrostrat
    macrostrat = importdataset("data/pregenerated_responses.tsv", '\t', importas=:Tuple)
    macro_cats = match_rocktype(macrostrat.rocktype, macrostrat.rockname, macrostrat.rockdescrip)

    known_rocks = macro_cats.sed .| macro_cats.ign .| macro_cats.met
    total_known = count(known_rocks)

    # Get elements
    majors, minors = get_elements()
    allelements = [majors; minors]


## --- Composition of exposed crust; start by getting samples with complete geochemical data
    # Number of samples measuring over 95% of rock components
    complete = Array{Float64}(undef, length(bulkidx), 1)
    for i in eachindex(bulkidx)
        if bulkidx[i] != 0
            for j in allelements
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
    majorelem, = get_elements()

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

## --- End of file