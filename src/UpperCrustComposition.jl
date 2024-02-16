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


## --- Figure out how many geochemical samples explain 90% of the matches 
    # If this is changed, remember to change the values in CalculateFlux.jl!!
    c = countmap(mbulk.Sample_ID)
    npoints = count(<(percentile(values(c), 90)), values(c))


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


## --- Terminal printout to copy paste into the latex table formatting sheet
    comp = NamedTuple{keys(class)}([(
        comp = round.([result[:,k][i] for i in eachindex(majors)], sigdigits=3),
        sem = round.([result_err[:,k][i]*2 for i in eachindex(majors)], sigdigits=1)
    ) for k in eachindex(keys(class))])


    @info """Major element composition by selected lithologic class (± 2 s.e):
             $(join(rpad.(majors, 8), " "))
    Bulk:    $(join(rpad.(comp.bulk.comp, 8), " "))
           ± $(join(rpad.(comp.bulk.sem, 8), " "))

    Sed:     $(join(rpad.(comp.sed.comp, 8), " "))
           ± $(join(rpad.(comp.sed.sem, 8), " "))

    Shale:   $(join(rpad.(comp.shale.comp, 8), " "))
           ± $(join(rpad.(comp.shale.sem, 8), " "))

    Ign:     $(join(rpad.(comp.ign.comp, 8), " "))
           ± $(join(rpad.(comp.ign.sem, 8), " "))       

    Plut:    $(join(rpad.(comp.plut.comp, 8), " "))
           ± $(join(rpad.(comp.plut.sem, 8), " "))

    Granite: $(join(rpad.(comp.granite.comp, 8), " "))
           ± $(join(rpad.(comp.granite.sem, 8), " "))

    Volc:    $(join(rpad.(comp.volc.comp, 8), " "))
           ± $(join(rpad.(comp.volc.sem, 8), " "))

    Basalt:  $(join(rpad.(comp.basalt.comp, 8), " "))
           ± $(join(rpad.(comp.basalt.sem, 8), " "))
    
    """


## --- Distribution of lithologies exposed at the Earth's surface 
    # Specifically, I want this to be directly comperable to the measurements of fractional
    # contribution to total denudation.
    # 
    # Total denudation uses the match_cats lithologies, where metamorphic rocks have 
    # been totally reassigned to their assumed protoliths. Therefore, while we want an 
    # estimate of the surficial abundance of metamorphic rocks, we want the total 
    # sedimentary and igneous abundances to add to 100%. Separate metamorphic into 
    # metasedimentary and metaigneous for funsies.

    # Major classes, where metamorphic is undifferentiated
    include_minor!(macro_cats);
    macro_cats.met .&= .!(macro_cats.sed .| macro_cats.ign);
    majorclass = (:sed, :ign, :met) 
    dist_major = NamedTuple{majorclass}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in majorclass])
    )

    # Minor classes
    exclude_minor!(macro_cats);
    dist_minorsed = NamedTuple{minorsed}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorsed]).*(dist_major.sed/100))
    dist_minorign = NamedTuple{minorign}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorign]).*(dist_major.ign/100))
    dist_minorvolc = NamedTuple{minorvolc}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorvolc]).*(dist_minorign.volc/100))
    dist_minorplut = NamedTuple{minorplut}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorplut]).*(dist_minorign.plut/100))

    # Export to file 
    sed = [[string.(collect(keys(dist_minorsed))) collect(values(dist_minorsed))];
    ["total sedimentary" dist_major.sed]]
    ign = [[string.(collect(keys(dist_minorign))) collect(values(dist_minorign))];
        ["total igneous" dist_major.ign]]
    volc = [[string.(collect(keys(dist_minorvolc))) collect(values(dist_minorvolc))];
        ["total volcanic" dist_minorign.volc]]
    plut = [[string.(collect(keys(dist_minorplut))) collect(values(dist_minorplut))];
        ["total plutonic" dist_minorign.plut]]
    met_total = ["undifferentiated metamorphic" dist_major.met]

    collected_abundance = [sed; ign; volc; plut; met_total];
    collected_abundance[:,2] = round.(collected_abundance[:,2], sigdigits=3)
    writedlm(surficial_abundance_out, collected_abundance)


## ---
    # (Option 1)
    # We want to calculate absolute surficial abundances with two outputs: 
    #
    # 1) Include "undifferentiated" as a type. This will allow people to get a sense of 
    # how many samples were randomly re-assigned during the matching process
    #
    # 2) Abundances of known types only. This will allow us to compare 
    # 
    # Absolute surficial abundances means that if an abundance is listed as 3%, it is 3%
    # of all rocks exposed on Earth, not 3% of that subclass. The alternative to this is 
    # to make the sum of the minor classes equal to 100, e.g. if a volcanic subclass is 
    # 3%, that means 3% of all volcanic rocks are that subclass.

    # Calculate the surficial abundance of metamorphic rocks, including metasedimentary
    # and metaigneous rocks
    met_undiff = ["total metamorphic" (count(macro_cats.met)/length(macro_cats.met))*100] 

    # Calculate abundances of major classes, where metamorphic rocks refers to 
    # undifferentiated metamorphic rocks
    include_minor!(macro_cats);
    macro_cats.met .&= .!(macro_cats.sed .| macro_cats.ign);
    majorclass = (:sed, :ign, :met) 
    dist_major = NamedTuple{majorclass}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in majorclass])
    )

    # Calculate abundances of minor classes
    exclude_minor!(macro_cats);
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:end];
    sed = (minorsed..., :sed);
    ign = (minorign..., :ign);
    volc = (minorvolc..., :volc);
    plut = (minorplut..., :plut);

    dist_minorsed = NamedTuple{sed}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in sed]).*(dist_major.sed/100))
    dist_minorign = NamedTuple{ign}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in ign]).*(dist_major.ign/100))
    dist_minorvolc = NamedTuple{volc}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in volc]).*(dist_minorign.volc/100))
    dist_minorplut = NamedTuple{plut}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in plut]).*(dist_minorign.plut/100))

    # Check 
    # (sum(dist_minorsed) + sum(dist_minorvolc) + sum(dist_minorplut) + 
    # dist_minorign.carbonatite + dist_minorign.ign + dist_major.met) == 100

    # Export data to a file 
    sed = [[string.(collect(keys(dist_minorsed))) collect(values(dist_minorsed))];
        ["total sedimentary" dist_major.sed]]
    sed[end-1, 1] = "undifferentiated sedimentary"
    
    ign = [[string.(collect(keys(dist_minorign))) collect(values(dist_minorign))];
        ["total igneous" dist_major.ign]]
    ign[end-1, 1] = "undifferentiated igneous"
    
    volc = [[string.(collect(keys(dist_minorvolc))) collect(values(dist_minorvolc))];
        ["total volcanic" dist_minorign.volc]]
    volc[end-1, 1] = "undifferentiated volcanic"

    plut = [[string.(collect(keys(dist_minorplut))) collect(values(dist_minorplut))];
        ["total plutonic" dist_minorign.plut]]
    plut[end-1, 1] = "undifferentiated plutonic"

    met_total = ["undifferentiated metamorphic" dist_major.met]

    collected_abundance = [sed; ign; volc; plut; met_total; met_undiff];
    collected_abundance[:,2] = round.(collected_abundance[:,2], sigdigits=3)
    writedlm(surficial_abundance_total_out, collected_abundance)


## --- Distribution of lithologies exposed at the Earth's surface (Option 2)
    # This option allows us to compare to the relative contribution of each class to 
    # global sediment production 

    # Major classes, where metamorphic is undifferentiated
    include_minor!(macro_cats);
    macro_cats.met .&= .!(macro_cats.sed .| macro_cats.ign);
    majorclass = (:sed, :ign, :met) 
    dist_major = NamedTuple{majorclass}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in majorclass])
    )

    # Minor classes
    exclude_minor!(macro_cats);
    dist_minorsed = NamedTuple{minorsed}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorsed]).*(dist_major.sed/100))
    dist_minorign = NamedTuple{minorign}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorign]).*(dist_major.ign/100))
    dist_minorvolc = NamedTuple{minorvolc}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorvolc]).*(dist_minorign.volc/100))
    dist_minorplut = NamedTuple{minorplut}(normalize!(
        [count(macro_cats[k])/length(macro_cats[k])*100 for k in minorplut]).*(dist_minorign.plut/100))

    # Export to file 
    sed = [[string.(collect(keys(dist_minorsed))) collect(values(dist_minorsed))];
    ["total sedimentary" dist_major.sed]]
    ign = [[string.(collect(keys(dist_minorign))) collect(values(dist_minorign))];
        ["total igneous" dist_major.ign]]
    volc = [[string.(collect(keys(dist_minorvolc))) collect(values(dist_minorvolc))];
        ["total volcanic" dist_minorign.volc]]
    plut = [[string.(collect(keys(dist_minorplut))) collect(values(dist_minorplut))];
        ["total plutonic" dist_minorign.plut]]
    met_total = ["undifferentiated metamorphic" dist_major.met]

    collected_abundance = [sed; ign; volc; plut; met_total];
    collected_abundance[:,2] = round.(collected_abundance[:,2], sigdigits=3)
    writedlm(surficial_abundance_out, collected_abundance)

    
## --- End of file