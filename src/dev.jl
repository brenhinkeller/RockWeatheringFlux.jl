"""

Calculate the average wt. % and flux by source of an element.

I want this to 
    - account for the fact that not all Macrostrat samples were matched, so some indices are 0
    - Be as generic as possible, i.e. don't do subindexing for a specific element (that should
        be done outside of the function)
    - return two NTuples, one with the wt.% and one with the flux by type

# Arguments
- `crustal_density::Number=2750`: Average crustal density (kg/m³).
"""
function flux_source(bulk::AbstractArray, bulkidx::Vector{Int64}, erosion::NamedTuple, 
    crustal_area::NamedTuple; crustal_density::Number=2750)

    # Preallocate
    allkeys = (:siliciclast, :shale, :carb, :chert, :evaporite, :coal, :sed,
        :volc, :plut, :ign,
        :metased, :metaign, :met,
        :cryst
    )
    allinitvals = fill(NaN, length(allkeys))
    
    wt = Dict(zip(allkeys, allinitvals))
    flux = Dict(zip(allkeys, allinitvals))

    # Get EarthChem samples


    return wt, flux
end