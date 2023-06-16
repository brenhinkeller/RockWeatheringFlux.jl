## --- Set up
    # Packages
    using StatGeochem
    using LoopVectorization
    using HDF5
    using DelimitedFiles
    using Measurements
    using Plots

    # Local Utilities
    include("Utilities.jl")

    # Load results. Assume no change in element order during file save
    res = h5open("output/rwf_output6.h5", "r")
    rocks = read(res["bulkrockflux"]["rocktypes"])      # Names of rock subtypes
    elems = read(res["element_names"])                  # Names of element / element oxides

    # kg / Gt conversion. File is in units of kg/yr
    const kg_gt = 1000000000000 # (kg/Gt)


## --- Compute bulk flux without differentiating by element
    # Load data into computer memory
    bulkflux = read(res["bulkrockflux"]["val"]) .± read(res["bulkrockflux"]["std"])     # Currently no std.
    bulkflux = NamedTuple{Tuple(Symbol.(rocks))}(bulkflux)

    # Before passing GO... make sure the data was saved and loaded correctly
    minorsed_flux = 0.0
    for i in minorsed
        minorsed_flux += bulkflux[i]
    end
    @assert minorsed_flux < bulkflux.sed

    minorign_flux = 0.0
    for i in minorign
        minorign_flux += bulkflux[i]
    end
    @assert minorign_flux < bulkflux.ign

    minormet_flux = 0.0
    for i in minormet
        minormet_flux += bulkflux[i]
    end
    @assert minormet_flux < bulkflux.met

    # Total global flux (denundation)
    bulkglobalflux = (bulkflux.sed + bulkflux.ign + bulkflux.met) / kg_gt
    @info "Total global denundation: $(round(bulkglobalflux, digits=2)) Gt/yr"



## --- Total global flux (kg/yr) for each element
    fluxpath = res["elementflux"]

    globalelem = read(fluxpath["totalelemflux"]["val"]) .± read(fluxpath["totalelemflux"]["std"])
    globalelem = NamedTuple{Tuple(Symbol.(elems))}(globalelem)


## --- Flux (kg/yr) for each element by rock type
    # Preallocate
    blank = Dict(zip(Symbol.(elems), fill(NaN ± NaN, length(elems))))
    blanktuple = NamedTuple{Tuple(keys(blank))}(values(blank))
    flux_cats = Dict(zip(Symbol.(rocks), fill(blanktuple, length(rocks))))
    
    # Calculate and store as a NamedTuple
    for i in eachindex(rocks)
        path = fluxpath["byrocktype"][rocks[i]]
        flux = blank

        for j in eachindex(elems)
            flux[Symbol(elems[j])] = read(path["flux_val"])[j] ± read(path["flux_std"])[j]
        end

        flux_cats[Symbol(rocks[i])] = NamedTuple{Tuple(keys(flux))}(values(flux))
    end
    flux_cats = NamedTuple{Tuple(keys(flux_cats))}(values(flux_cats))

    # File will not be used anymore
    close(res)


## --- Relative flux (fraction) for each element by rock type
    # Preallocate
    relative_cats = Dict(zip(Symbol.(rocks), fill(blanktuple, length(rocks))))

    # Relative contribution is the absolute by rock type over the absolute global flux
    for i in eachindex(rocks)
        flux = blank
        rock = Symbol(rocks[i])     # Cannot assume element order!

        for j in eachindex(elems)
            element = Symbol(elems[j])
            flux[element] = flux_cats[rock][element] / globalelem[element]
        end

        relative_cats[Symbol(rocks[i])] = NamedTuple{Tuple(keys(flux))}(values(flux))
    end
    relative_cats = NamedTuple{Tuple(keys(relative_cats))}(values(relative_cats))


## --- 


## --- TABLE: Flux by rocktype by element (Gt), ignoring error for now
    bigmatrix = Array{Union{String, Float64}}(undef, length(elems), length(rocks))
    for i = 1:length(flux_cats)
        for j = 1:length(flux_cats[i])
            bigmatrix[j, i] = flux_cats[i][j].val / kg_gt
        end
    end

    # Get global flux for each element, making sure elements are in the right order
    globalflux = Array{Float64}(undef, length(elems), 1)
    keyorder = keys(flux_cats.sed)
    for i in eachindex(keyorder)
        globalflux[i] = netflux_elem[keyorder[i]].val / kg_gt
    end

    bigmatrix = vcat(reshape(string.(collect(keys(flux_cats))), 1, 14), bigmatrix)
    bigmatrix = hcat(bigmatrix, ["global"; globalflux])
    bigmatrix = hcat([""; string.(collect(keys(flux_cats.sed)))], bigmatrix)
    writedlm("output/flux_gt.csv", bigmatrix)

    # Same idea, but by relative contribution instead of total mass
    bigmatrix[2:end, 2:end] .= bigmatrix[2:end, 2:end] ./ globalflux
    writedlm("output/flux_relativecontrib.csv", bigmatrix)



    

    
## --- Calculate P provenance
    # planning to clean this code up later
    # path = res["elementflux"]["byrocktype"]

    # sed = read(path["sed"]["flux_val"])
    # ign = read(path["ign"]["flux_val"])
    # met = read(path["met"]["flux_val"])

    # totalp = sed[elems.P2O5] + ign[elems.P2O5] + met[elems.P2O5]
    # spv = sed[elems.P2O5]/totalp
    # ipv = ign[elems.P2O5]/totalp
    # mpv = met[elems.P2O5]/totalp

    # chunkval = [spv, ipv, mpv]

    # # Realizing I don't actually know how to deal with getting error from an absolute
    # # value to a value that makes sense for fractional contribution
    # sed_e = read(path["sed"]["flux_std"])[elems.P2O5]
    # ign_e = read(path["ign"]["flux_std"])[elems.P2O5]
    # met_e = read(path["met"]["flux_std"])[elems.P2O5]

    # chunkerr = [sed_e, ign_e, met_e]

    # # Make plot
    # x = [1, 2, 3]
    # b = plot(x, chunkval, seriestype=:bar, framestyle=:box, xticks=(x, ["Sed", "Ign", "Met"]), 
    #     label="", ylabel="Fractional contribution to P flux"
    # )
    # savefig(b, "prelim_P.pdf")

## --- End of File