## --- Set up 
    # Visualize simulated crust compositions as a function of assumed volatiles

    # Load data and base packages
    include("Definitions.jl")

    # Define simulation filepaths
    # sim3 = "src/volatile_sensitivity/simout_prop2.h5"
    # sim3 = "src/volatile_sensitivity/simout_prop2_gard.h5"
    sim3 = "src/volatile_sensitivity/simout_prop2_combo.h5"


## --- Load simulation with proportional assumed volatiles (gypsum:dolomite:bassanite)
    fid = h5open(sim3, "r")
    allelements = read(fid["vars"]["elements"])
    sᵢ = findfirst(x->x=="SiO2", allelements)

    # Preallocate 
    simid = keys(fid["sims"]["UCC"])
    sim_prop2 = Array{Float64}(undef, length(simid), 1)
    silica_prop2 = Array{Float64}(undef, length(simid), 1)
    stdev_prop2 = Array{Float64}(undef, length(simid), 1)
    c = 0

    # Find bulk silica distribution for iterative experiments
    for i in eachindex(simid)
        silica_prop2[i] = read(fid["sims"]["UCC"][simid[i]])[sᵢ]
        stdev_prop2[i] = read(fid["sims"]["UCC_std"][simid[i]])[sᵢ]

        try 
            sim_prop2[i] = parse(Float64, split(simid[i],"_")[2]) 
        catch
            # Catch current experimental conditions
            c = i
            sim_prop2[i] = 47
        end
    end
    close(fid)

    # Convert s.d. to 2 s.e.
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)

    npoints = unique_sample(mbulk.Sample_ID, 90)
    stdev_prop2 ./= sqrt(npoints)

    # Sort data from smallest to largest permitted volatile (x) value
    p = sortperm(sim_prop2, dims=1)
    c = findfirst(x->x==c, p)

    sim_prop2 = sim_prop2[p]
    silica_prop2 = silica_prop2[p]
    stdev_prop2 = stdev_prop2[p]


## --- Plot simulations
    # Build plot 
    h = Plots.plot(
        framestyle=:box, 
        grid = false,
        fontfamily=:Helvetica,
        ylabel="Average Crustal Silica [wt.%]",
        xlabel="Maximum Permitted Sedimentary Volatiles [wt.%]",
        xlims=(10,100),
        ylims=(35,70),
        legend=:bottomleft,
        fg_color_legend=:white,
        legendfont=9
    );

    # Plot other estimates for comparison
    Plots.hline!([66.62], label="", 
        linewidth=2, linestyle=:dot, color=:dimgrey)
    Plots.annotate!([(97, 68.62, text("Rudnick and Gao, 2014", 9, :right, :top, :dimgrey))])
    Plots.hline!([59.5], label="", 
        linewidth=2, linestyle=:dot, color=:dimgrey)
    Plots.annotate!([(97, 61.5, text("Pease et al., 2023 [EarthChem prior]", 9, :right, :top, :dimgrey))])

    # Plot dolomite volatiles
    Plots.vline!([47], label="", linewidth=2, linestyle=:dash,
        color=colors_covariance.a)
    Plots.annotate!([(44, 40, text("Dolomite", 9, :left, :top, colors_covariance.a, rotation=90))])

    # Proportional volatiles (carbonate / evaporite / gypsum)
    Plots.plot!(sim_prop2, silica_prop2, ribbon=stdev_prop2, 
        markershape=:circle,
        color=colors_covariance.b, msc=:auto,
        fillalpha=0.15,
        # label="Sedimentary ≠ Evaporite ≠ Gypsum"
        label="Simulations"
    )

    # Current experimental conditions 
    Plots.scatter!([sim_prop2[c]], [silica_prop2[c]], # yerror=stdev[c], 
        markershape=:star5,
        markersize=10,
        linewidth=1,
        color=colors.rhyolite, linecolor=:black, msc=:black,
        label="This Study"
    )

    display(h)
    savefig(h, "$filepath/volatile_sensitivity.pdf")


## --- End of file