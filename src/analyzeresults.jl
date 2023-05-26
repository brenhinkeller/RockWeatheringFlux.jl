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

    # Load results
    res = h5open("output/rwf_output3.h5", "r")

    # kg / Gt conversion
    const kg_gt = 1000000000000 # (kg/Gt)


## --- Get data from HDF5 file into computer memory
    # Organize rock types and element types for easy indexing
    rocks = read(res["bulkrockflux"]["rocktypes"])
    elems = read(res["element_names"])

    # Total flux for all rocks
    # Note that we don't have variance for global flux--only comes from wt.% rn
    netflux_bulk = read(res["bulkrockflux"]["val"]) .± read(res["bulkrockflux"]["std"])
    netflux_bulk = NamedTuple{Tuple(Symbol.(rocks))}(netflux_bulk)

    # Total global flux of each element
    path = res["elementflux"]["totalelemflux"]
    netflux_elem = read(path["val"]) .± read(path["std"])
    netflux_elem = NamedTuple{Tuple(Symbol.(elems))}(netflux_elem)

    # Flux and wt.% of each element by rocktype
    flux_cats = Dict{Symbol, NamedTuple}()
    wt_cats = Dict{Symbol, NamedTuple}()

    for i in eachindex(rocks)
        path = res["elementflux"]["byrocktype"][rocks[i]]
        fluxelem = Dict{Symbol, Measurement{Float64}}()
        wtelem = Dict{Symbol, Measurement{Float64}}()

        for j in eachindex(elems)
            fluxelem[Symbol(elems[j])] = read(path["flux_val"])[j] ± read(path["flux_std"])[j]
            wtelem[Symbol(elems[j])] = read(path["wtpct_val"])[j] ± read(path["wtpct_std"])[j]
        end
        flux_cats[Symbol(rocks[i])] = NamedTuple{Tuple(keys(fluxelem))}(values(fluxelem))
        wt_cats[Symbol(rocks[i])] = NamedTuple{Tuple(keys(wtelem))}(values(wtelem))
    end
    flux_cats = NamedTuple{Tuple(keys(flux_cats))}(values(flux_cats))
    wt_cats = NamedTuple{Tuple(keys(wt_cats))}(values(wt_cats))
    close(res)

## --- Start by making some tables...
    # Flux by rocktype by element (Gt), ignoring error for now
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


## --- Calculate total global denundation rate!
    bulkglobalflux = (netflux_bulk.sed + netflux_bulk.ign + netflux_bulk.met) / kg_gt

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