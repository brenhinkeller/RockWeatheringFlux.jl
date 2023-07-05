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
    res = h5open("output/erodedmaterial.h5", "r")
    rocks = read(res["bulkrockflux"]["rocktypes"])      # Names of rock subtypes
    elems = read(res["element_names"])                  # Names of element / element oxides

    # kg / Gt conversion. File is in units of kg/yr
    const kg_gt = 1000000000000 # (kg/Gt)


## --- Compute bulk flux without differentiating by element
    # Load data into computer memory
    bulkflux = read(res["bulkrockflux"]["val"]) .± read(res["bulkrockflux"]["std"])     # Currently no std.
    bulkflux = NamedTuple{Tuple(Symbol.(rocks))}(bulkflux)

    # Before passing GO... make sure the data was saved and loaded correctly
    # As a note to this and all the others... when running this as a script it won't work
    # because of minorsed_flux inside the loop is local. Haha
    # minorsed_flux = 0.0
    # for i in minorsed
    #     minorsed_flux += bulkflux[i]
    # end
    # @assert minorsed_flux < bulkflux.sed

    # minorign_flux = 0.0
    # for i in minorign
    #     minorign_flux += bulkflux[i]
    # end
    # @assert minorign_flux < bulkflux.ign

    # minormet_flux = 0.0
    # for i in minormet
    #     minormet_flux += bulkflux[i]
    # end
    # @assert minormet_flux < bulkflux.met

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


## --- Write absolute flux to a .csv file
    # Preallocate
    bigmatrix = Array{Union{Float64, String}}(undef, length(elems) + 1, length(rocks) + 2)

    # Convert data to a write-able matrix, converting units to Gt/yr
    for i in 1:length(flux_cats)
        for j in 1:length(flux_cats[i])
            bigmatrix[j+1, i+1] = flux_cats[i][j].val / kg_gt
        end
    end

    # Add global flux by element, converting units to Gt/yr
    row_order = collect(keys(flux_cats.sed))
    for i in eachindex(row_order)
        nextelement = row_order[i]
        bigmatrix[i+1, end] = globalelem[nextelement].val / kg_gt
    end

    # Add row names (elements) and column names (rock types)
    for i = 1:(size(bigmatrix)[1] - 1)
        bigmatrix[i+1, 1] = string.(row_order[i])
    end

    header = vcat("", string.(collect(keys(flux_cats))), "global")
    for i = 1:(size(bigmatrix)[2])
        bigmatrix[1, i] = string.(header[i])
    end

    # Write to file
    writedlm("output/flux_gt.csv", bigmatrix, delim=',')


## --- Write relative flux to a .csv file 
    # Preallocate
    bigmatrix = Array{Union{Float64, String}}(undef, length(elems) + 1, length(rocks) + 1)

    # Convert data to a write-able matrix
    for i in 1:length(relative_cats)
        for j in 1:length(relative_cats[i])
            bigmatrix[j+1, i+1] = relative_cats[i][j].val
        end
    end

    # Add row names (elements) and column names (rock types)
    for i = 1:(size(bigmatrix)[1] - 1)
        bigmatrix[i+1, 1] = string.(row_order[i])
    end

    header = vcat("", string.(collect(keys(relative_cats))), "global")
    for i = 1:(size(bigmatrix)[2])
        bigmatrix[1, i] = string.(header[i])
    end

    # Write to file
    writedlm("output/flux_relativecontrib.csv", bigmatrix, ',')


## --- End of File