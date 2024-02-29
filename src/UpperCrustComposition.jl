## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using StatsBase
    using NNLS      # Unregistered package


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


## --- Compute and export composition of exposed crust!
    # Some elements just aren't measured, but a good whole rock geochemistry should
    # measure the major elements. If it's a NaN, it's probably just not there fr fr
    for i in majors
        zeronan!(mbulk[i])
    end

    # Save to a file 
    result = Array{Float64}(undef, (length(allelements), length(class)))
    result_err = similar(result)
    rows = string.(allelements)
    cols = hcat("element", reshape(string.(collect(keys(class))), 1, :))
    
    for i in eachindex(keys(class))
        result[:,i] .= [nanmean(mbulk[j][class[i]]) for j in allelements]
        result_err[:,i] .= [nanstd(mbulk[j][class[i]])./sqrt(npoints) for j in allelements]
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
        sem = round.([result_err[:,k][i]*2 for i in eachindex(majors)], sigdigits=1)
    ) for k in eachindex(keys(class))])

    # Pre-computed compositions
    target = (:volc, :plut, :sed, :bulk)
    out = fill("", length(majors)+1)
    for t in target
        for i in eachindex(majors) 
            out[i] *= "\$ $(comp[t].comp[i]) \\pm $(comp[t].sem[i]) \$; "
        end
        out[end] *= "$(round(sum(comp[t].comp), sigdigits=4)); "
    end

    # Anhydrous majors, doing The Most(™) so it's immune to changes in element order
    k = findfirst(x->x==:bulk, keys(class))
    index = collect(1:length(majors))[1:end .!= findfirst(x->x==:Volatiles, majors)]
    anhydrous_comp = [result[:,k][i] for i = index]
    anhydrous_sem = [result_err[:,k][i]*2 for i = index]
    
    # Normalize, without the function so I can make sure we're propagating error correctly
    # Also round at the same time
    sum_a = sum(anhydrous_comp)
    anhydrous_comp = round.(anhydrous_comp ./ sum_a .* 100, sigdigits=3)
    anhydrous_sem = round.(anhydrous_sem ./ sum_a .* 100, sigdigits=1)

    out_anh = fill("", length(majors)+1)
    for i = index
        out_anh[i] *= "\$ $(anhydrous_comp[i]) \\pm $(anhydrous_sem[i]) \$"
    end
    out_anh[end] = string(round(sum(anhydrous_comp), sigdigits=4))

    # Print to terminal 
    for i in eachindex(out)
        println("$(out[i] * out_anh[i])")
    end

    # # Alternate views
    # @info """Major element composition by selected lithologic class (± 2 s.e):
    #          $(join(rpad.(majors, 8), " "))
    # Bulk:    $(join(rpad.(comp.bulk.comp, 8), " "))
    #        ± $(join(rpad.(comp.bulk.sem, 8), " "))

    # Sed:     $(join(rpad.(comp.sed.comp, 8), " "))
    #        ± $(join(rpad.(comp.sed.sem, 8), " ")) 

    # Plut:    $(join(rpad.(comp.plut.comp, 8), " "))
    #        ± $(join(rpad.(comp.plut.sem, 8), " "))

    # Volc:    $(join(rpad.(comp.volc.comp, 8), " "))
    #        ± $(join(rpad.(comp.volc.sem, 8), " "))
    # """


## --- Compute mixing ratios for other estimates 
    # Hypothesis: other estimates can be approximated by mixing our estimates in different 
    # proportions. Separate into sed / plut / volc.
    # 
    # For comparison, we want to also show the mix for our bulk dataset, but because this 
    # is just a measure of surficial abundance, we want to normalize values by this number

    # In order: SiO2, Al2O3, FeOT, TiO2, MgO, CaO, Na2O, K2O, Volatiles
    this_study = NamedTuple{keys(class)}([result[:,k][i] for i in eachindex(majors)] 
        for k in eachindex(keys(class))
    )
    anhydrous = copy(this_study.bulk); anhydrous[findfirst(x->x==:Volatiles, majors)] = 0
    pease = [59.5, 17.3, 6.53, 0.901, 2.88, 5.39, 4.14, 3.24, 0.129]
    rudnickgao = [66.62, 15.4, 5.04, 0.64, 2.48, 3.59, 3.27, 2.8, 0]
    gao = [58.48, 12.12, 4.6, 0.57, 3.73, 7.41, 2.57, 2.27, 7.36]

    # Set mixing endmembers
    endmembers = [this_study.sed this_study.volc this_study.plut]

    # Calculate percent of each endmember represented in each estimate
    mix_this_study = nnls(endmembers, this_study.bulk)
    mix_anhydrous = nnls(endmembers, anhydrous)
    mix_pease = nnls(endmembers, pease)
    mix_rudnickgao = nnls(endmembers, rudnickgao)
    mix_gao = nnls(endmembers, gao)

    # Calculate misfit
    misfit_this_study = sum((endmembers * mix_this_study - this_study.bulk).^2)
    misfit_anhydrous = sum((endmembers * mix_anhydrous - anhydrous).^2)
    misfit_pease = sum((endmembers * mix_pease - pease).^2)
    misfit_rudnickgao = sum((endmembers * mix_rudnickgao - rudnickgao).^2)
    misfit_gao = sum((endmembers * mix_gao - gao).^2)

    # Normalize to surficial abundance of each endmember 
    mix_this_study_norm = mix_this_study ./ mix_this_study
    mix_anhydrous_norm = mix_anhydrous ./ mix_this_study
    mix_pease_norm = mix_pease ./ mix_this_study
    mix_rudnickgao_norm = mix_rudnickgao ./ mix_this_study
    mix_gao_norm = mix_gao ./ mix_this_study
    
    # Format for LaTeX-formatting Excel sheet (bulk, anhydrous, and other estimates)
    out = fill("", size(endmembers)[2] + 1)
    out = hcat(
        round.([mix_this_study; misfit_this_study], sigdigits=3),
        round.([mix_anhydrous; misfit_anhydrous], sigdigits=3),
        round.([mix_rudnickgao; misfit_rudnickgao], sigdigits=3),
        round.([mix_gao; misfit_gao], sigdigits=3),
        round.([mix_pease; misfit_pease], sigdigits=3),
    )
    for i in 1:size(out)[1]
        println(join(out[i,:], ";"))
    end
    
    
## --- End of file