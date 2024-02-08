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
    # npoints = count(>(percentile(values(c), 95)), values(c))


## --- Compute and export composition of exposed crust!
    # Some elements just aren't measured, but a good whole rock geochemistry should
    # measure the major elements. If it's a NaN, it's probably just not there fr fr
    for i in majors
        zeronan!(mbulk[i])
    end

    # We don't want metamorphics, but we do want all samples 
    target = deleteat!(collect(keys(match_cats)), findall(x->x==:met, collect(keys(match_cats))))
    class = merge(
        NamedTuple{Tuple(target)}(match_cats[k] for k in target), 
        (bulk=trues(length(match_cats[1])),)
    )

    # Save to a file 
    result = Array{Float64}(undef, (length(allelements), length(class)))
    rows = string.(allelements)
    cols = hcat("element", reshape(string.(collect(keys(class))), 1, :))
    
    for i in eachindex(keys(class))
        result[:,i] .= [nanmean(mbulk[j][class[i]]) for j in allelements]
    end
    writedlm("$ucc_out", vcat(cols, hcat(rows, result)))

    # Terminal printout
    majorcomp = round.([result[:,end][i] for i in eachindex(majors)], digits=1)

    @info """Bulk crustal composition ($geochem_fid | $macrostrat_io):
    $(join(majors, " \t "))
    $(join(majorcomp, " \t "))

    Total (majors): $(round(nansum(result[:,end][1:length(majors)]), sigdigits=4))%
    Total (major + trace): $(round(nansum(result[:,end]), sigdigits=4))%
    """


## --- End of file