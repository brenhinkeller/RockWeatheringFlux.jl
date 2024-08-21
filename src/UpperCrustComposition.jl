## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using StatsBase
    using NNLS                      # Unregistered package


## --- Load data
    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)

    # Elements of interest
    majors, minors = get_elements()
    allelements = [majors; minors]

    # How many samples explain 90% of the matches?
    npoints = unique_sample(mbulk.Sample_ID, 90)
    @info """ Matched sample metadata:
    Total: $(length(unique(mbulk.Sample_ID))) samples 
    90% of matches explained by $npoints samples"""


## --- Compute and export composition of exposed crust!
    # Some elements just aren't measured, but a good whole rock geochemistry should
    # measure the major elements. If it's a NaN, it's probably just not there fr fr
    # All other elements should have NaNs and not zeros though
    for i in majors
        zeronan!(mbulk[i])
    end
    for i in minors 
        nanzero!(mbulk[i])
    end

    # Save to a file 
    result = Array{Float64}(undef, (length(allelements), length(class)))
    result_err = similar(result)
    rows = string.(allelements)
    cols = hcat("element", reshape(string.(collect(keys(class))), 1, :))
    
    for i in eachindex(keys(class))
        result[:,i] .= [nanmean(mbulk[j][class[i]]) for j in allelements]
        result_err[:,i] .= [nanstd(mbulk[j][class[i]])./sqrt(npoints).*2 for j in allelements]
    end
    writedlm("$ucc_out", vcat(cols, hcat(rows, result)))
    writedlm("$ucc_out_err", vcat(cols, hcat(rows, result_err)))

    
## --- Terminal printout ± 2 SEM
    majorcomp = round.([result[:,end][i] for i in eachindex(majors)], digits=1)
    majorcomp_err = round.([result_err[:,end][i]*2 for i in eachindex(majors)], sigdigits=1)

    @info """Bulk crustal composition ($geochem_fid | $macrostrat_io):
      $(join(rpad.(majors, 8), " "))
      $(join(rpad.(majorcomp, 8), " "))
    ± $(join(rpad.(majorcomp_err, 8), " "))

    Total (majors): $(round(nansum(result[:,end][1:length(majors)]), sigdigits=4))%
    Total (major + trace): $(round(nansum(result[:,end]), sigdigits=4))%
    """


## --- Terminal printout to copy paste into the LaTeX-formatting Excel sheet 
    comp = NamedTuple{keys(class)}([(
        comp = round.([result[:,k][i] for i in eachindex(majors)], sigdigits=3),
        sem = round.([result_err[:,k][i] for i in eachindex(majors)], sigdigits=1)
    ) for k in eachindex(keys(class))])

    # Carbonate SiO2 and CaO values for the paper :)
    i_silica = findfirst(x->x==:SiO2, majors)
    i_calcium = findfirst(x->x==:CaO, majors)
    @info """Carbonates:
    SiO2: $(comp.carb.comp[i_silica]) \\pm $(comp.carb.sem[i_silica])
    CaO:  $(comp.carb.comp[i_calcium]) \\pm $(comp.carb.sem[i_calcium])
    """

    # Sedimentary, volcanic, and plutonic report only wt.% without error
    target = (:sed, :volc, :plut,)
    out = fill("", length(majors)+1)
    for t in target
        for i in eachindex(majors) 
            out[i] *= "$(comp[t].comp[i]); "
        end
        out[end] *= "$(round(sum(comp[t].comp), sigdigits=4)); "
    end

    # Whole-earth estimate reported with error
    for i in eachindex(majors)
        out[i] *= "\$ $(comp.bulk.comp[i]) \\pm $(comp.bulk.sem[i]) \$; "
    end
    out[end] *= "$(round(sum(comp.bulk.comp), sigdigits=4)); "

    println("composition of sed / volc / plut / bulk")
    for i in eachindex(out)
        println("$(out[i])")
    end


## --- Compute mixing ratios for other crust composition estimates
    # Hypothesis: other estimates can be approximated by mixing our estimates in different 
    # proportions. Separate into sed / plut / volc.

    # In order: SiO2, Al2O3, FeOT, TiO2, MgO, CaO, Na2O, K2O, Volatiles
    this_study = NamedTuple{keys(class)}([result[:,k][i] for i in eachindex(majors)] 
        for k in eachindex(keys(class))
    )
    pease = [59.5, 17.3, 6.53, 0.901, 2.88, 5.39, 4.14, 3.24, 0.129]
    rudnickgao = [66.62, 15.4, 5.04, 0.64, 2.48, 3.59, 3.27, 2.8, 0]
    gao = [58.48, 12.12, 4.6, 0.57, 3.73, 7.41, 2.57, 2.27, 7.36]

    # Set mixing endmembers
    endmembers = [this_study.sed this_study.volc this_study.plut]

    # Calculate percent of each endmember represented in each estimate
    mix_this_study = nnls(endmembers, this_study.bulk)
    mix_pease = nnls(endmembers, pease)
    mix_rudnickgao = nnls(endmembers, rudnickgao)
    mix_gao = nnls(endmembers, gao)

    # Calculate misfit
    misfit_this_study = sum((endmembers * mix_this_study - this_study.bulk).^2)
    misfit_pease = sum((endmembers * mix_pease - pease).^2)
    misfit_rudnickgao = sum((endmembers * mix_rudnickgao - rudnickgao).^2)
    misfit_gao = sum((endmembers * mix_gao - gao).^2)
    
    # Format for LaTeX-formatting Excel sheet (bulk, anhydrous, and other estimates)
    out = fill("", size(endmembers)[2] + 1)
    out = hcat(
        round.([mix_this_study.*100; misfit_this_study], sigdigits=3),
        round.([mix_rudnickgao.*100; misfit_rudnickgao], sigdigits=3),
        round.([mix_gao.*100; misfit_gao], sigdigits=3),
        round.([mix_pease.*100; misfit_pease], sigdigits=3),
    )
    println("sed / volc / plut mixing ratios for: bulk / R&G / Gao / Pease")
    for i in 1:size(out)[1]
        println(join(out[i,:], ";"))
    end
    

## --- Compute mixing ratios for our eroded material estimates 
    # Load data and calculate anhydrous composition
    fid = readdlm(erodedcomp_out)
    bulk_i = findfirst(x->x=="bulk", fid)[2]
    bulkcomp = float.(fid[2:length(majors)+1, bulk_i])

    anhydcomp = copy(bulkcomp)
    anhydcomp[end] = 0
    normalize!(anhydcomp)

    # Muller et al.
    muller = [65.1,18.7,5.67,0.8,2.1,2.9,1.1,2.7,0]

    # Compute endmember mixing and calculate misfit :)
    mix_bulk = nnls(endmembers, bulkcomp)
    mix_anhyd = nnls(endmembers, anhydcomp)
    mix_muller = nnls(endmembers, muller)

    misfit_bulk = sum((endmembers * mix_bulk - bulkcomp).^2)
    misfit_anhyd = sum((endmembers * mix_anhyd - anhydcomp).^2)
    misfit_muller = sum((endmembers * mix_muller - muller).^2)

    # Print to terminal 
    out = fill("", size(endmembers)[2] + 1)
    out = hcat(
        round.([mix_bulk.*100; misfit_bulk], sigdigits=3),
        round.([mix_anhyd.*100; misfit_anhyd], sigdigits=3),
        round.([mix_muller.*100; misfit_muller], sigdigits=3),
    )
    println("sed / volc / plut mixing ratios for: bulk / anhydrous / Muller")
    for i in 1:size(out)[1]
        println(join(out[i,:], ";"))
    end

## --- End of file