# Visualize the change in silica distributions when the EarthChem dataset is restricted
# to only samples with between 84-100% silica.

## --- Set up
    # Packages
    using StatGeochem
    using MAT
    using HDF5
    using Plots

    # Local utilities
    using ProgressMeter, LoopVectorization, Measurements, Static
    include("../utilities/Utilities.jl")

    # Rock types we're interested in
    target = (:siliciclast, :shale, :carb, :sed,)


## --- Load unrestricted bulk data
    bulk = matread("data/bulk.mat")["bulk"];
    bulk_init = NamedTuple{Tuple(Symbol.(keys(bulk)))}(values(bulk));

    fid = h5open("output/bulk_unrestricted_types.h5", "r")
    header = read(fid["bulk_cats_head"])
    data = read(fid["bulk_cats"])
    data = @. data > 0
    bulk_cats_init = NamedTuple{Tuple(Symbol.(header))}(
        [data[:,i] for i in eachindex(header)])
    close(fid)

    # We only want to look at samples above sea level
    etopo = h5read("data/etopo/etopo1.h5", "vars/elevation")
    elev = find_etopoelev(etopo, bulk_init.Latitude, bulk_init.Longitude)
    abovesea = elev .> -140;


## --- Load intermediate file
    # Consider silica distribution for rocks without the assumed volatiles
    fid = h5open("output/intermediate_screen.h5", "r")
        bulkweight = read(fid["vars"]["bulkweight"])
        tᵢ = @. 84 <= bulkweight <= 104
    close(fid)


## --- Load restricted bulk data
    fid = h5open("output/bulk.h5", "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    bulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])

    header = read(fid["bulktypes"]["bulk_cats_head"])
    data = read(fid["bulktypes"]["bulk_cats"])
    data = @. data > 0
    bulk_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i] for i in eachindex(header)])


## --- Major rock types should include minor rock types
    # But we only really care about sedimentary rocks
    minorsed, minorign, minormet = get_minor_types()
    for type in minorsed
        bulk_cats_init.sed .|= bulk_cats_init[type]
        bulk_cats.sed .|= bulk_cats[type]
    end


## --- Resample (spatial weighting) unrestricted and restricted silica
    nsims = Int(1e6)

    # Unrestricted
    simout_init = NamedTuple{target}(Array{Float64}(undef, nsims) for _ in target)
    for t in target
        # Filter to samples of interest above sea level
        f = bulk_cats_init[t] .& abovesea;
        data = bulk_init.SiO2[f]
        uncert = fill(1.0, length(data))

        k = invweight_location(bulk_init.Latitude[f], bulk_init.Longitude[f])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        simout_init[t] .= bsresample(data, uncert, nsims, p)
    end

    # Restricted without assumed volatiles
    simout_novol = NamedTuple{target}(Array{Float64}(undef, nsims) for _ in target)
    for t in target
        f = bulk_cats_init[t] .& abovesea .& tᵢ;
        data = bulk_init.SiO2[f]
        uncert = fill(1.0, length(data))

        k = invweight_location(bulk_init.Latitude[f], bulk_init.Longitude[f])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        simout_novol[t] .= bsresample(data, uncert, nsims, p)
    end

    # Restricted
    simout = NamedTuple{target}(Array{Float64}(undef, nsims) for _ in target)
    for t in target
        # Filter to samples of interest (already filtered to above sea level)
        f = bulk_cats[t];
        data = bulk.SiO2[f]
        uncert = fill(1.0, length(data))

        k = invweight_location(bulk.Latitude[f], bulk.Longitude[f])
        p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        simout[t] .= bsresample(data, uncert, nsims, p)
    end


## --- Compare filtered and unfiltered silica distributions
    fig = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(fig)
        r = target[i]

        h = stephist(bulk_init.SiO2[bulk_cats_init[r] .& abovesea], bins=1:100, 
            normalize=:pdf, label="Initial Distribution", color=:black, linewidth=2)
        stephist!(simout_init[r], bins=1:100, normalize=:pdf, 
            label="Initial Distribution (Resampled)", color=:blue, linewidth=2)
        stephist!(simout_novol[r], bins=1:100, normalize=:pdf,
            label="Filtered without Assumed Volatiles", color=:red, linewidth=2)
        stephist!(h, simout[r], bins=1:100, normalize=:pdf, 
            label="Filtered", color=:green, linewidth=2,
            title="$r")
        fig[i] = h
    end

    nrows = ceil(Int, length(target)/2)
    h = plot(fig..., layout=(nrows,2), framestyle=:box, size=(1200, nrows*400), 
        titleloc=:left, titlefont = font(15),
        legendfontsize = 10, fg_color_legend=:white, legend=false,
        left_margin = (25,:px),
    )

    #  Make a legend
    leg = Plots.plot(yticks=:none, xticks=:none, framestyle=:none, 
        legendfontsize = 15, fg_color_legend=:white, legend=:inside)
    Plots.plot!(leg, [0],[0], label="Initial Distribution, Above Sea Level Only", 
        color=:black, linewidth=2)
    Plots.plot!(leg, [0],[0], label="Initial Distribution (Resampled)", color=:blue, 
        linewidth=2)
    Plots.plot!(leg, [0],[0], label="Filtered without Assumed Volatiles (Resampled)", 
        color=:red, linewidth=2)
    Plots.plot!(leg, [0],[0], label="Filtered Data (Resampled)", color=:green, 
        linewidth=2)

    # Make a plot layout
    l = @layout [
        a{0.9h} 
        b{0.1h} 
    ]
    h = Plots.plot(h, leg, layout = l, bottom_margin=(30,:px))
    display(h)


## --- End of File 