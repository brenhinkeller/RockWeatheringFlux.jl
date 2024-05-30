## --- Set up 
    # Look at the relative abundance of undifferentiated rocks over time 

    # Packages 
    using RockWeatheringFlux
    using DelimitedFiles, HDF5
    using Plots

    # Save figures to: 
    filepath = "results/figures/burial"

    # Definitions 
    xmin, xmax, nbins = 0, 3800, 38


## --- Load data 
    # Filter cats looks at rocks which are both mapped and matched as that rock type
    # Macro cats looks at rocks which are mapped as that rock type 

    # fid = h5open("src/burial/resampled_geochem.h5", "r")
    # head_data = Tuple(Symbol.(read(fid["vars"]["wt"]["head_data"])))
    # head_cats = Tuple(Symbol.(read(fid["vars"]["wt"]["head_cats"])))

    # data = read(fid["vars"]["wt"]["data"]["data"])
    # cats = read(fid["vars"]["wt"]["cats"]["macro_cats"])
    # cats = @. cats > 0

    # simba = NamedTuple{head_data}([data[:,i] for i in eachindex(head_data)])
    # macro_cats = NamedTuple{head_cats}([cats[:,i] for i in eachindex(head_cats)])
    # close(fid)

    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Macrostrat and mapped rock classes
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        rocktype = read(fid["vars"]["rocktype"])[t],
        rockdescrip = read(fid["vars"]["rockdescrip"])[t],
        rockname = read(fid["vars"]["rockname"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    include_minor!(macro_cats)
    macro_cats = delete_cover(macro_cats)
    close(fid)
    

## --- Propagate age uncertainty by resampling
    # Age and uncertainty (5% or 50 Ma.)
    sampleage = copy(macrostrat.age)
    ageuncert = nanadd.(macrostrat.agemax, .- macrostrat.agemin) ./ 2;
    for i in eachindex(ageuncert)
        ageuncert[i] = nanmaximum([ageuncert[i], sampleage[i] .* 0.05, 50])
    end

    # Rock class 
    classes = keys(macro_cats)
    macro_in = Array{Int64}(undef, length(macro_cats.sed), length(classes))
    for i in eachindex(classes)
        for j in eachindex(macro_cats[classes[i]])
            macro_in[j,i] = ifelse(macro_cats[classes[i]][j], 1, 0)
        end
    end

    # Resampling weights of 1
    p = ones(length(sampleage))

    # Resample! 
    data = hcat(sampleage, macro_in)
    uncert = hcat(ageuncert, zeros(size(macro_in)))
    resampled = bsresample(data, uncert, Int(1e5), p)

    # Re-parse 
    simage = resampled[:,1]
    sim_cats = delete_cover(get_cats(false, size(resampled)[1])[2]);
    for i = 2:size(resampled)[2]
        sim_cats[classes[i-1]] .= resampled[:,i] .> 0
    end


## --- Exclude undifferentiated rocks 
    exclude_minor!(sim_cats)
    undiff_cats = (
        sed = copy(sim_cats.sed),
        ign = copy(sim_cats.ign),
        volc = copy(sim_cats.volc),
        plut = copy(sim_cats.plut),
    )

    include_minor!(sim_cats)
    # macro_cats.sed .&= .!sed_undiff
    # macro_cats.ign .&= .!ign_undiff
    # macro_cats.volc .&= .!volc_undiff
    # macro_cats.plut .&= .!plut_undiff


## --- [PLOT] Fraction of rocks which are undifferentiated 
    # Base plot
    h = plot(
        xlabel="Age [Ma.]", 
        ylabel="Undifferentiated Fraction",
        # yaxis=:log10,
        ylims=(0,1),
        xlims=(xmin, xmax),
        framestyle=:box,
        fontfamily=:Helvetica,
        fg_color_legend=:white,
        legend=:topright,
        titleloc=:left,
        left_margin=(40,:px), right_margin=(25,:px), bottom_margin=(40,:px),
    );

    # Define targets and make plot
    target = (:sed, :ign, :volc, :plut)
    figs = Array{Plots.Plot{Plots.GRBackend}}(undef, length(target))
    for i in eachindex(target)
        k = target[i]
        
        # Samples per bin and fraction of undifferentiated samples 
        c, nₖ = bincounts(simage[sim_cats[k]], xmin, xmax, nbins*2)
        c, nᵤ = bincounts(simage[undiff_cats[k]], xmin, xmax, nbins*2)
        n = nᵤ ./ nₖ
        # n[isnan.(n)] .= 1
        # n[n .== 0] .= 1

        hᵢ = deepcopy(h)
        plot!(hᵢ, c, n,
            label="", 
            title="$k",
            color=colors[k], lcolor=colors[k], msc=:auto,
            # markershape=:circle,
            # linewidth=2,
            seriestype=:bar,
            barwidths=((xmax-xmin)/(nbins*2)),
        )
        figs[i] = hᵢ
    end

    h = plot(figs..., layout=(2,2), size=(1200,800))
    savefig(h, "$filepath/undifferentiated.pdf")


## --- End of file 