## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using StatsBase


## --- Load data for the matched EarthChem samples
    # Indices and classes of matched samples
    fid = readdlm(matchedbulk_io)
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    match_cats = match_rocktype(string.(vec(fid[:,2]))[t]);
    include_minor!(match_cats)
    match_cats = delete_cover(match_cats)

    # Metamorphic samples
    fid = h5open("$macrostrat_io", "r")
    header = Tuple(Symbol.(read(fid["type"]["macro_cats_head"])))
    data = read(fid["type"]["metamorphic_cats"])
    data = @. data > 0
    metamorphic_cats = NamedTuple{header}([data[:,i][t] for i in eachindex(header)])
    include_minor!(metamorphic_cats)
    match_cats.met .|= (metamorphic_cats.sed .| metamorphic_cats.ign)

    # Lithology
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    include_minor!(macro_cats)
    macro_cats = delete_cover(macro_cats)
    close(fid)

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(fid)

    # Elements of interest
    majors, minors = get_elements()
    allelements = [majors; minors]


## --- Figure out how many geochemical samples explain 90% of the matches 
    c = countmap(mbulk.Sample_ID)
    npoints = count(>(percentile(values(c), 90)), values(c))
    # npoints = count(<(percentile(values(c), 95)), values(c))


## --- Compute and export composition of exposed crust!
    # Some elements just aren't measured, but a good whole rock geochemistry should
    # measure the major elements. If it's a NaN, it's probably just not there fr fr
    for i in majors
        zeronan!(mbulk[i])
    end

    # We want an option to filter for all samples 
    class = merge(match_cats, (bulk=trues(length(match_cats[1])),))

    # Save to a file 
    result = Array{Float64}(undef, (length(allelements), length(class)))
    result_err = similar(result)
    rows = string.(allelements)
    cols = hcat("element", reshape(string.(collect(keys(class))), 1, :))
    
    for i in eachindex(keys(class))
        result[:,i] .= [nanmean(mbulk[j][class[i]]) for j in allelements]
        result_err[:,i] .= [nanstd(mbulk[j][class[i]])./npoints for j in allelements]
    end
    writedlm("$ucc_out", vcat(cols, hcat(rows, result)))
    writedlm("$ucc_out_err", vcat(cols, hcat(rows, result_err)))

    # Terminal printout
    majorcomp = round.([result[:,end][i] for i in eachindex(majors)], digits=1)
    majorcomp_err = round.([result_err[:,end][i] for i in eachindex(majors)], sigdigits=1)

    @info """Bulk crustal composition ($geochem_fid | $macrostrat_io):
      $(join(rpad.(majors, 8), " "))
      $(join(rpad.(majorcomp, 8), " "))
    ± $(join(rpad.(majorcomp_err, 8), " "))

    Total (majors): $(round(nansum(result[:,end][1:length(majors)]), sigdigits=4))%
    Total (major + trace): $(round(nansum(result[:,end]), sigdigits=4))%
    """

    # A little comparative analysis, as a treat
    # [reshape(collect(keys(match_cats)), :, 1) collect(1:29)]  # What column do you want
    plutcomp = round.([result[:,end-4][i] for i in eachindex(majors)], digits=1);
    shalecomp = round.([result[:,2][i] for i in eachindex(majors)], digits=1);
    clastcomp = round.([result[:,1][i] for i in eachindex(majors)], digits=1);

    @info """Major element composition by selected lithologic class:
             $(join(rpad.(majors, 8), " "))
    Bulk:    $(join(rpad.(majorcomp, 8), " "))
    Plut:    $(join(rpad.(plutcomp, 8), " "))
    Shale:   $(join(rpad.(shalecomp, 8), " "))
    Clastic: $(join(rpad.(clastcomp, 8), " "))
    """

## --- Distribution of lithologies exposed at the Earth's surface 
    # We want to figure out the distributions where "undifferentiated" is considered a 
    # minor type. e.g. igneous might be 40% volcanic / 30% plutonic / 1% carbonatite /  
    # 29% undifferentiated
    include_minor!(macro_cats);
    majorclass = (:sed, :ign, :met) 
    dist_major = NamedTuple{majorclass}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in majorclass])
    )

    # Calculate "absolute surficial abundances", i.e. if an abundance of a volcanic 
    # subclass is 3%, that means 3% of all exposed bedrock is that subclass.
    # (sum(dist_minorsed) + sum(dist_minorvolc) + sum(dist_minorplut) + 
    # dist_minorign.carbonatite + dist_minorign.ign + dist_major.met) == 100

    # The alternative to this is to make the sum of the minor classes equal to 100, e.g.
    # if a volcanic subclass is 3%, that means 3% of all volcanic rocks are that subclass.

    macro_cats.met .&= .!(macro_cats.sed .| macro_cats.ign);
    exclude_minor!(macro_cats);
    minorsed, minorvolc, minorplut, minorign = get_rock_class()[2:end];
    sed = (minorsed..., :sed);
    ign = (minorign..., :ign);
    volc = (minorvolc..., :volc);
    plut = (minorplut..., :plut);

    dist_minorsed = NamedTuple{sed}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in sed]).*(dist_major.sed/100)
    )
    dist_minorign = NamedTuple{ign}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in ign]).*(dist_major.ign/100)
    )
    dist_minorvolc = NamedTuple{volc}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in volc]).*(dist_minorign.volc/100)
    )
    dist_minorplut = NamedTuple{plut}(
        normalize!([count(macro_cats[k])/length(macro_cats[k])*100 for k in plut]).*(dist_minorign.plut/100)
    )

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

    met_total = ["total metamorphic" dist_major.met]
    met_undiff = ["undifferentiated metamorphic" (count(macro_cats.met)/length(macro_cats.met))] 

    collected_abundance = [sed; ign; volc; plut; met_total; met_undiff];
    collected_abundance[:,2] = round.(collected_abundance[:,2], sigdigits=4)
    writedlm(surficial_abundance_out, collected_abundance)


## --- End of file