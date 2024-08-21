## --- Set up
    # Plot volcanic P/Alk and F-org over geologic time 

    # Load data and base packages
    include("Definitions.jl")

    # Make the volcanic color visible
    colorvolc = pal[5]

    # Definitions 
    xmin, xmax, nbins = 0, 3800, 38


## --- Load data
    # Carbon isotope data 
    fid = h5open("src/visualization/burial/resampled_carbon.h5", "r")
    head = keys(fid["vars"]["carb"])
    sim_carb = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["carb"]["$k"]) for k in head)
    head = keys(fid["vars"]["org"])
    sim_org = NamedTuple{Tuple(Symbol.(head))}(read(fid["vars"]["org"]["$k"]) for k in head)
    close(fid)

    # Resampled geochemical data (phosphorus / alkalinity ratios)
    fid = h5open("src/visualization/burial/resampled_geochem.h5", "r")
    head = read(fid["vars"]["ratio"]["head_ratio"])
    head_class = keys(fid["vars"]["ratio"]["data"])
    sim_ratio = NamedTuple{Tuple(Symbol.(head_class))}(NamedTuple{Tuple(Symbol.(head))}(
            read(fid["vars"]["ratio"]["data"]["$k"])[:,i] for i in eachindex(head)
        ) for k in head_class
    )
    close(fid)


## --- Calculate fraction of carbon buried as organic 
    mantle = -5.5
        
    c,m,e = binmeans(sim_carb.age, sim_carb.d13c_carb, xmin, xmax, nbins, relbinwidth=2)
    carbonate = m .± e

    c,m,e = binmeans(sim_org.age, sim_org.d13c_org_corrected, xmin, xmax, nbins)
    organic = m .± e

    frog = (mantle .- carbonate) ./ (organic .- carbonate)
    frog = (;
        val = Measurements.value.(frog),
        err = Measurements.uncertainty.(frog),
    )


## --- Plot!! 
    h = Plots.plot(
        xlabel="Age [Ma.]", ylabel="Phosphorus / Alkalinity [mol. ratio]",
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(20,:px), right_margin=(20,:px),

    );
    Plots.plot!(sim_ratio.volc.c, sim_ratio.volc.m,
        label="",
        color=colorvolc, lcolor=colorvolc, msc=:auto,
            y_foreground_color_border=colorvolc,
            y_foreground_color_text=colorvolc,
            y_foreground_color_axis=colorvolc,
            y_guidefontcolor=colorvolc,
        linestyle=:dot,
        linewidth=2,
    )
    Plots.plot!(sim_ratio.volc.c, sim_ratio.volc.m,
        yerror=(2*sim_ratio.volc.el, 2*sim_ratio.volc.eu),
        label="",
        color=colorvolc, lcolor=colorvolc, msc=:auto,
        markershape=:circle,
        seriestype=:scatter,
        linewidth=2,
    )
    Plots.plot!(twinx(), c, frog.val, 
        ribbon=2*frog.err,
        label="", ylabel="Fraction Buried as Organic Carbon",
        color=pal[9],
            y_foreground_color_border=pal[8],
            y_foreground_color_text=pal[8],
            y_foreground_color_axis=pal[8],
            y_guidefontcolor=pal[8],
        linewidth=2,
    )
    vline!([541, 717, 2500], label="", linestyle=:dash, color=:black)
    display(h)
    savefig(h, "$filepath/primaryproduction.pdf")
    savefig(h, "$filepath/primaryproduction.png")
    

## --- End of file