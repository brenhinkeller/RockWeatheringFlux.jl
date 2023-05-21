"""

Calculate the average wt. % and flux (kg/yr) by source of an element.
also returns global_flux

I want this to 
    - account for the fact that not all Macrostrat samples were matched, so some indices are 0
    - Be as generic as possible, i.e. don't do subindexing for a specific element (that should
        be done outside of the function)
    - return two NTuples, one with the wt.% and one with the flux by type
    - printout how many samples were not NaN
    - erosion, bulk_cats, crustal_area must contain the keys :siliciclast, :shale, :carb, :chert, 
    :evaporite, :coal, :sed,
    :volc, :plut, :ign,
    :metased, :metaign, :met,
    :cryst (others ok)

# Arguments
- `crustal_density::Number=2750`: Average crustal density (kg/m³).
- `elem::String=""`: Element being analyzed, for terminal printout
"""
function flux_source(bulk::AbstractArray, bulkidx::Vector{Int64}, erosion::NamedTuple, 
    bulk_cats::NamedTuple, crustal_area::NamedTuple; 
    crustal_density::Number=2750, elem::String="")

    # Preallocate
    allkeys = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal, :sed,
        :volc, :plut, :ign,
        :metased, :metaign, :met,
        :cryst
    )
    allinitvals = fill(NaN, length(allkeys))
    npoints = length(bulkidx)

    wt = Dict(zip(allkeys, allinitvals))
    flux = Dict(zip(allkeys, allinitvals))
    bulkdata = Array{Float64}(undef, npoints, 1)
    datacount = 0

    # Get EarthChem samples
    for i in eachindex(bulkidx)
        notzero = bulkidx[i] != 0
        if notzero
            bulkdata[i] = bulk[bulkidx[i]]
            datacount += 1
        else
            bulkdata[i] = NaN
        end
    end

    @info "$datacount of $npoints $elem samples ($(round(Int, datacount/npoints))%) are not NaN"

    # Calculate average wt.% for each rock type
    for i in eachindex(allkeys)
        wt[allkeys[i]] = nanmean(bulkdata[bulk_cats[i]])
    end
    wt = NamedTuple{Tuple(allkeys)}(values(wt))

    # Calculate provenance by rock type
    for i in eachindex(allkeys)
        flux[allkeys[i]] = erosion[i] * wt[i] * crustal_area[i] * crustal_density * 1e-6
    end
    flux = NamedTuple{Tuple(allkeys)}(values(flux))
    global_flux = flux.sed + flux.ign + flux.met

    return wt, flux, global_flux
end