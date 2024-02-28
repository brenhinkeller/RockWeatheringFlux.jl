## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using StatsBase


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
    for i in eachindex(out)
        println(out[i])
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
    for i in eachindex(out_anh)
        println(out_anh[i])
    end

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
    mix_this_study = endmembers \ this_study.bulk
    mix_anhydrous = endmembers \ anhydrous
    mix_pease = endmembers \ pease
    mix_rudnickgao = endmembers \ rudnickgao
    mix_gao = endmembers \ gao

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
    display(out)


## --- Distribution of lithologies exposed at the Earth's surface 
#     # Specifically, I want this to be directly comperable to the measurements of fractional
#     # contribution to total denudation.
#     # 
#     # Total denudation uses the match_cats lithologies, where metamorphic rocks have 
#     # been totally reassigned to their assumed protoliths. Therefore, while we want an 
#     # estimate of the surficial abundance of metamorphic rocks, we want the total 
#     # sedimentary and igneous abundances to add to 100%. Separate metamorphic into 
#     # metasedimentary and metaigneous for funsies.

#     # Major classes, where metamorphic is undifferentiated
#     include_minor!(macro_cats);
#     macro_cats.met .&= .!(macro_cats.sed .| macro_cats.ign);
#     majorclass = (:sed, :ign, :met) 
#     dist_major = NamedTuple{majorclass}(
#         normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in majorclass])
#     )

#     # Minor classes
#     exclude_minor!(macro_cats);
#     dist_minorsed = NamedTuple{minorsed}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorsed]).*(dist_major.sed/100))
#     dist_minorign = NamedTuple{minorign}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorign]).*(dist_major.ign/100))
#     dist_minorvolc = NamedTuple{minorvolc}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorvolc]).*(dist_minorign.volc/100))
#     dist_minorplut = NamedTuple{minorplut}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorplut]).*(dist_minorign.plut/100))

#     # Export to file 
#     sed = [[string.(collect(keys(dist_minorsed))) collect(values(dist_minorsed))];
#     ["total sedimentary" dist_major.sed]]
#     ign = [[string.(collect(keys(dist_minorign))) collect(values(dist_minorign))];
#         ["total igneous" dist_major.ign]]
#     volc = [[string.(collect(keys(dist_minorvolc))) collect(values(dist_minorvolc))];
#         ["total volcanic" dist_minorign.volc]]
#     plut = [[string.(collect(keys(dist_minorplut))) collect(values(dist_minorplut))];
#         ["total plutonic" dist_minorign.plut]]
#     met_total = ["undifferentiated metamorphic" dist_major.met]

#     collected_abundance = [sed; ign; volc; plut; met_total];
#     collected_abundance[:,2] = round.(collected_abundance[:,2], sigdigits=3)
#     writedlm(surficial_abundance_out, collected_abundance)


# ## ---
#     # (Option 1)
#     # We want to calculate absolute surficial abundances with two outputs: 
#     #
#     # 1) Include "undifferentiated" as a type. This will allow people to get a sense of 
#     # how many samples were randomly re-assigned during the matching process
#     #
#     # 2) Abundances of known types only. This will allow us to compare 
#     # 
#     # Absolute surficial abundances means that if an abundance is listed as 3%, it is 3%
#     # of all rocks exposed on Earth, not 3% of that subclass. The alternative to this is 
#     # to make the sum of the minor classes equal to 100, e.g. if a volcanic subclass is 
#     # 3%, that means 3% of all volcanic rocks are that subclass.

#     # Calculate the surficial abundance of metamorphic rocks, including metasedimentary
#     # and metaigneous rocks
#     met_undiff = ["total metamorphic" (count(macro_cats.met)/length(macro_cats.met))*100] 

#     # Calculate abundances of major classes, where metamorphic rocks refers to 
#     # undifferentiated metamorphic rocks
#     include_minor!(macro_cats);
#     macro_cats.met .&= .!(macro_cats.sed .| macro_cats.ign);
#     majorclass = (:sed, :ign, :met) 
#     dist_major = NamedTuple{majorclass}(
#         normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in majorclass])
#     )

#     # Calculate abundances of minor classes
#     exclude_minor!(macro_cats);
#     minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:end];
#     sed = (minorsed..., :sed);
#     ign = (minorign..., :ign);
#     volc = (minorvolc..., :volc);
#     plut = (minorplut..., :plut);

#     dist_minorsed = NamedTuple{sed}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in sed]).*(dist_major.sed/100))
#     dist_minorign = NamedTuple{ign}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in ign]).*(dist_major.ign/100))
#     dist_minorvolc = NamedTuple{volc}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in volc]).*(dist_minorign.volc/100))
#     dist_minorplut = NamedTuple{plut}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in plut]).*(dist_minorign.plut/100))

#     # Check 
#     # (sum(dist_minorsed) + sum(dist_minorvolc) + sum(dist_minorplut) + 
#     # dist_minorign.carbonatite + dist_minorign.ign + dist_major.met) == 100

#     # Export data to a file 
#     sed = [[string.(collect(keys(dist_minorsed))) collect(values(dist_minorsed))];
#         ["total sedimentary" dist_major.sed]]
#     sed[end-1, 1] = "undifferentiated sedimentary"
    
#     ign = [[string.(collect(keys(dist_minorign))) collect(values(dist_minorign))];
#         ["total igneous" dist_major.ign]]
#     ign[end-1, 1] = "undifferentiated igneous"
    
#     volc = [[string.(collect(keys(dist_minorvolc))) collect(values(dist_minorvolc))];
#         ["total volcanic" dist_minorign.volc]]
#     volc[end-1, 1] = "undifferentiated volcanic"

#     plut = [[string.(collect(keys(dist_minorplut))) collect(values(dist_minorplut))];
#         ["total plutonic" dist_minorign.plut]]
#     plut[end-1, 1] = "undifferentiated plutonic"

#     met_total = ["undifferentiated metamorphic" dist_major.met]

#     collected_abundance = [sed; ign; volc; plut; met_total; met_undiff];
#     collected_abundance[:,2] = round.(collected_abundance[:,2], sigdigits=3)
#     writedlm(surficial_abundance_total_out, collected_abundance)


# ## --- Distribution of lithologies exposed at the Earth's surface (Option 2)
#     # This option allows us to compare to the relative contribution of each class to 
#     # global sediment production 

#     # Major classes, where metamorphic is undifferentiated
#     include_minor!(macro_cats);
#     macro_cats.met .&= .!(macro_cats.sed .| macro_cats.ign);
#     majorclass = (:sed, :ign, :met) 
#     dist_major = NamedTuple{majorclass}(
#         normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in majorclass])
#     )

#     # Minor classes
#     exclude_minor!(macro_cats);
#     dist_minorsed = NamedTuple{minorsed}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorsed]).*(dist_major.sed/100))
#     dist_minorign = NamedTuple{minorign}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorign]).*(dist_major.ign/100))
#     dist_minorvolc = NamedTuple{minorvolc}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorvolc]).*(dist_minorign.volc/100))
#     dist_minorplut = NamedTuple{minorplut}(normalize!(
#         [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorplut]).*(dist_minorign.plut/100))

#     # Export to file 
#     sed = [[string.(collect(keys(dist_minorsed))) collect(values(dist_minorsed))];
#     ["total sedimentary" dist_major.sed]]
#     ign = [[string.(collect(keys(dist_minorign))) collect(values(dist_minorign))];
#         ["total igneous" dist_major.ign]]
#     volc = [[string.(collect(keys(dist_minorvolc))) collect(values(dist_minorvolc))];
#         ["total volcanic" dist_minorign.volc]]
#     plut = [[string.(collect(keys(dist_minorplut))) collect(values(dist_minorplut))];
#         ["total plutonic" dist_minorign.plut]]
#     met_total = ["undifferentiated metamorphic" dist_major.met]

#     collected_abundance = [sed; ign; volc; plut; met_total];
#     collected_abundance[:,2] = round.(collected_abundance[:,2], sigdigits=3)
#     writedlm(surficial_abundance_out, collected_abundance)

    
## --- End of file