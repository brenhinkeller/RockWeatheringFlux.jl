## --- Set up
    # Packages
    using RockWeatheringFlux
    using DelimitedFiles, HDF5


## --- Load data for the matched EarthChem samples
    # Indices of matched EarthChem samples from SampleMatch.jl
    fid = readdlm("$matchedbulk_io")
    bulkidx = Int.(vec(fid[:,1]))
    t = @. bulkidx != 0
    
    # Macrostrat (for rock types)
    fid = h5open("$macrostrat_io", "r")
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)

    include_minor!(macro_cats)

    # Earthchem (geochemical data)
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

    UCC = (
        bulk = [nanmean(bulk[i]) for i in allelements],
        sed = [nanmean(bulk[i][macro_cats.sed]) for i in allelements],
        met = [nanmean(bulk[i][macro_cats.met]) for i in allelements],
        ign = [nanmean(bulk[i][macro_cats.ign]) for i in allelements],
        volc = [nanmean(bulk[i][macro_cats.volc]) for i in allelements],
        plut = [nanmean(bulk[i][macro_cats.plut]) for i in allelements],
        granite = [nanmean(bulk[i][macro_cats.granite]) for i in allelements],
        basalt = [nanmean(bulk[i][macro_cats.basalt]) for i in allelements],
    )

    # Save to file
    rows = string.(allelements)
    cols = hcat("element", string.(reshape(collect(keys(UCC)), (1, length(keys(UCC))))))
    results = Array{Float64}(undef, (length(allelements), length(UCC)))
    for i in eachindex(keys(UCC))
        results[:,i] = UCC[i]
    end
    writedlm("$ucc_out", vcat(cols, hcat(rows, results)))

    # Terminal printout
    majorcomp = round.([UCC.bulk[i] for i in eachindex(majors)], digits=1)

    @info """Bulk crustal composition ($geochem_fid | $macrostrat_io):
    $(join(majors, " \t "))
    $(join(majorcomp, " \t "))

    Total (majors): $(round(nansum(UCC.bulk[1:length(majors)]), sigdigits=4))%
    Total (major + trace): $(round(nansum(UCC.bulk), sigdigits=4))%
    """


## --- End of file