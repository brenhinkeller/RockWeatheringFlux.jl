## -- Set up
    # Model erosion as a function of slope for major rock types 

    # Packages 
    using RockWeatheringFlux
    using HDF5, DelimitedFiles
    using Plots 
    using Isoplot: yorkfit

    # using StatsBase
    # using CurveFit; using Isoplot


## --- Load data 
    # Matched samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0

    # Lithologic class 
    match_cats, metamorphic_cats, class, megaclass = get_lithologic_class();

    # # Matched geochemical data
    # fid = h5open(geochem_fid, "r")
    # header = read(fid["bulk"]["header"])
    # data = read(fid["bulk"]["data"])
    # mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    # close(fid)

    # Mapped data 
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
    )
    close(fid)
    
    # Sample age and uncertainty (5% or 50 Myr. uncertainty)
    # sampleage, ageuncert = resampling_age(
    #     mbulk.Age, mbulk.Age_Min, mbulk.Age_Max, 
    #     macrostrat.age, macrostrat.agemin, macrostrat.agemax, 
    #     uncert_rel=5, uncert_abs=50
    # )

    # Instead of using the geochemical age, we should use the mapped age since we're 
    # assuming a spatial relationship between age and erosion (e.g., that continental arcs 
    # are the primary erosion source)
    sampleage = copy(macrostrat.age)
    ageuncert = nanadd.(macrostrat.agemax, .- macrostrat.agemin) ./ 2;
    ageuncert[isnan.(ageuncert)] .= sampleage[isnan.(ageuncert)] .* 0.05

    # t = isnan.(sampleage);
    # sampleage[t] = mbulk.Age[t]
    # ageuncert[t] .= nanadd.(mbulk.Age_Max[t], .- mbulk.Age_Min[t]) ./ 2;

    # for i in eachindex(ageuncert)
    #     ageuncert[i] = nanmax(sampleage[i]*0.05, ageuncert[i], 50)
    # end



## --- Slope / erosion rate at each coordinate point
    # Load the slope variable from the SRTM15+ maxslope file
    srtm15_slope = h5read("output/srtm15plus_maxslope.h5", "vars/slope")
    
    srtm15_sf = h5read("output/srtm15plus_maxslope.h5", "vars/scalefactor")

    # Get slope at each coordinate point
    # Function returns the standard deviation of slope in each window, which we don't
    # actually care about propagating
    rockslope = movingwindow(srtm15_slope, macrostrat.rocklat, macrostrat.rocklon, srtm15_sf, n=5)
    rockslope = (;
        vals = Measurements.value.(rockslope),
        errs = Measurements.uncertainty.(rockslope),
    )
    
    # Calculate all erosion rates (mm/kyr) (propagate uncertainty)
    rock_ersn = emmkyr.(rockslope.vals);
    rock_ersn = (;
        vals = Measurements.value.(rock_ersn),
        errs = Measurements.uncertainty.(rock_ersn),
    )


## --- Resample (temporal) erosion / age relationship 
    # Definitions
    nsims = 10_000
    xmin, xmax, nbins = 0, 3800, 38
    c = cntr(xmin:(xmax-xmin)/nbins:xmax)
    
    # Preallocate
    target = (:sed, :volc, :plut)
    simout_pctile = NamedTuple{target}(Array{Float64}(undef, nbins, 3) for _ in target)
    simout_bulk = NamedTuple{target}(Array{Float64}(undef, nsims, 2) for _ in target)

    # Resample
    @showprogress for key in target
        s = .!isnan.(sampleage) .& match_cats[key]

        # Resampling weights 
        # k = invweight_age(sampleage[s])
        # p = 1.0 ./ ((k .* nanmedian(5.0 ./ k)) .+ 1.0)
        p = ones(count(s))

        # # Resample 
        # c,m,el,eu = bin_bsr(sampleage[s], rockslope.vals[s], xmin, xmax, nbins,
        #     x_sigma = ageuncert[s],
        #     y_sigma = rockslope.errs[s],
        #     nresamplings = nsims,
        #     sem = :pctile,
        #     p = p,
        # )
        # simout_pctile[key] .= [m eu el]

        simout_bulk[key] .= bsresample([sampleage[s] rockslope.vals[s]], 
            [ageuncert[s] rockslope.errs[s]], nsims, p
        )
    end


## --- Plot data 
    # Resampled and binned data
    h = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlims=(-50, 3850),
        grid=false,
        fg_color_legend=:white,
        xlabel="Age [Ma.]", ylabel="Slope",
        title="bin_bsr", titleloc=:left,
    );
    c = cntr(xmin:(xmax-xmin)/nbins:xmax)
    for k in keys(simout_pctile)
        plot!(h, c, simout_pctile[k][:,1], 
            # yerror=(simout_pctile.sed[:,2], simout_pctile.sed[:,3],),
            label="$k", color=colors[k], seriestype=:path, markershape=:circle
        )
    end
    display(h)

    # Resampled and unbinned data 
    h = Plots.plot(
        framestyle=:box,
        fontfamily=:Helvetica,
        xlims=(-50, 3850),
        grid=false,
        fg_color_legend=:white,
        xlabel="Age [Ma.]", ylabel="Slope",
        title="bsresample + binmeans", titleloc=:left,
    );

    for k in keys(simout_bulk)
        c, m, e = binmeans(simout_bulk[k][:,1], simout_bulk[k][:,2], xmin, xmax, nbins)
        plot!(h, c, m, # yerror=e, 
            label="$k", color=colors[k], seriestype=:path, markershape=:circle
        )
    end
    display(h)
    
    


## --- Model age / erosion as with an exponential decay function 
    # Fit a linear model in log space to account for uncertainty in both age and erosion 

    fobj = yorkfit(c, simout.sed[:,1], m, (simout.sed[:,2] .+ simout.sed[:,3])./2)

## --- Fit an exponential decay to the age / slope relationship
    x = 0:10:3800
    xₛ = xᵢ = xₘ = x

    # xₛ = 200:10:3800
    c,m,e = binmeans(simsed[:,c_age], simsed[:,c_slp], xₛ[1], xₛ[end], Int((xₛ[end] - xₛ[1])/100))
    sedexp = exp_fit(c, m)
    yₛ = @. sedexp[1] * exp(sedexp[2] * xₛ);

    # xᵢ = 0:10:2900
    c,m,e = binmeans(simign[:,c_age], simign[:,c_slp], xᵢ[1], xᵢ[end], Int((xₛ[end] - xₛ[1])/100))
    ignexp = exp_fit(c, m)
    yᵢ = @. ignexp[1] * exp(ignexp[2] * xᵢ);

    # xₘ = 0:10:2300
    c,m,e = binmeans(simmet[:,c_age], simmet[:,c_slp], xₘ[1], xₘ[end], Int((xₛ[end] - xₛ[1])/100))
    metexp = exp_fit(c, m)
    yₘ = @. metexp[1] * exp(metexp[2] * xₘ);




## --- Yorkfit in log space

    # Transform to log-space
    ersn_sed, ersn_sed_e = unmeasurementify(emmkyr.(simsed[:,c_slp]))
    ersn_ign, ersn_ign_e = unmeasurementify(emmkyr.(simign[:,c_slp]))
    ersn_met, ersn_met_e = unmeasurementify(emmkyr.(simmet[:,c_slp]))

    c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    xσ = ones(length(c))
    fobj = yorkfit(collect(c), xσ, m, e)
    esed(x) = x * (fobj.slope) + (fobj.intercept)

    c,m,e = binmeans(simign[:,c_age], ersn_ign, 0, 3800, 38)
    fobj = yorkfit(collect(c), xσ, m, e)
    eign(x) = x * (fobj.slope) + (fobj.intercept)

    c,m,e = binmeans(simmet[:,c_age], ersn_met, 0, 3800, 38)
    fobj = yorkfit(collect(c), xσ, m, e)
    emet(x) = x * (fobj.slope) + (fobj.intercept)


## --- Erosion rate over time in linear space
    x = collect(0:10:3800)
    c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    hₛ = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
        label="Sed", markershape=:circle,)
    Plots.plot!(x, unmeasurementify(esed.(x))[1], label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simign[:,c_age], simign[:,c_slp], 0, 3800, 38)
    hᵢ = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
        label="Ign", ylabel="Hillslope [m/km]", markershape=:circle,)
    Plots.plot!(x, unmeasurementify(eign.(x))[1], label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simmet[:,c_age], simmet[:,c_slp], 0, 3800, 38)
    hₘ = Plots.plot(c, m, yerror=e, color=:orange, lcolor=:orange, msc=:orange,
        label="Met", xlabel="Bedrock Age [Ma]", markershape=:circle,)
    Plots.plot!(x, unmeasurementify(emet.(x))[1], label="Model", color=:black, linewidth=2)

    h = Plots.plot(hₛ, hᵢ, hₘ, layout=(3,1), size=(600, 1200), left_margin=(30, :px), 
        framestyle=:box, legend=:topright
    )
    display(h)
    

## --- Exponential fit to erosion
    ersn_sed, ersn_sed_e = unmeasurementify(emmkyr.(simsed[:,c_slp]))
    ersn_ign, ersn_ign_e = unmeasurementify(emmkyr.(simign[:,c_slp]))
    ersn_met, ersn_met_e = unmeasurementify(emmkyr.(simmet[:,c_slp]))

    x = 0:10:3800

    c,m,e = binmeans(simsed[:,c_age], ersn_sed, 0, 3800, 38)
    sedexp = exp_fit(c, m)
    yₛ = @. sedexp[1] * exp(sedexp[2] * x);
    hₛ = Plots.plot(c, m, yerror=e, color=:blue, lcolor=:blue, msc=:blue,
        label="Sed", markershape=:circle,)
    Plots.plot!(x, yₛ, label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simign[:,c_age], ersn_ign, 0, 3800, 38)
    ignexp = exp_fit(c, m)
    yᵢ = @. ignexp[1] * exp(ignexp[2] * x);
    hᵢ = Plots.plot(c, m, yerror=e, color=:red, lcolor=:red, msc=:red,
        label="Ign", ylabel="Erosion [m/Myr]", markershape=:circle,)
    Plots.plot!(x, yᵢ, label="Model", color=:black, linewidth=2)

    c,m,e = binmeans(simmet[:,c_age], ersn_met,  0, 3800, 38)
    metexp = exp_fit(c, m)
    yₘ = @. metexp[1] * exp(metexp[2] * x);
    hₘ = Plots.plot(c, m, yerror=e, color=:orange, lcolor=:orange, msc=:orange,
        label="Met", xlabel="Bedrock Age [Ma]", markershape=:circle,)
    Plots.plot!(x, yₘ, label="Model", color=:black, linewidth=2)

    h = Plots.plot(hₛ, hᵢ, hₘ, layout=(3,1), size=(600, 1200), left_margin=(30, :px), 
        framestyle=:box, legend=:topright,
    )
    display(h)

    
## --- End of file 