## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5


## --- Load data for the matched EarthChem samples
    # Indices and classes of matched samples
    fid = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0

    match_cats = match_rocktype(string.(vec(fid[:,2]))[t]);
    include_minor!(match_cats)

    # Matched geochemical data
    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][bulkidx[t]] for i in eachindex(header)])
    close(fid)

    # Elements of interest
    majors, minors = get_elements()
    allelements = [majors; minors]


## --- Compute and export composition of exposed crust!
    # Some elements just aren't measured, but a good whole rock geochemistry should
    # measure the major elements. If it's a NaN, it's probably just not there fr fr
    for i in majors
        zeronan!(bulk[i])
    end

    # Save to a file 
    rows = string.(allelements)
    cols = hcat("element", string.(reshape(collect(keys(match_cats)), 1, :)), "bulk")
    results = Array{Float64}(undef, (length(allelements), length(match_cats)+1))
    for i in eachindex(keys(match_cats))
        results[:,i] .= [nanmean(bulk[j][match_cats[i]]) for j in allelements]
    end

    bulkearth = [nanmean(bulk[i]) for i in allelements]
    results[:,end] .= bulkearth

    # Save to file
    writedlm("$ucc_out", vcat(cols, hcat(rows, results)))

    # Terminal printout
    majorcomp = round.([bulkearth[i] for i in eachindex(majors)], digits=1)

    @info """Bulk crustal composition ($geochem_fid | $macrostrat_io):
    $(join(majors, " \t "))
    $(join(majorcomp, " \t "))

    Total (majors): $(round(nansum(bulkearth[1:length(majors)]), sigdigits=4))%
    Total (major + trace): $(round(nansum(bulkearth), sigdigits=4))%
    """


## --- End of file