#=
The purpose of this script is to allow local plotting. Slope and location data is supplied 
by the intermediate "parsed_responses.tsv" file generated in CalculateFlux.jl. Extraction / 
parsing of phosphorous data from bulk.mat is replicated from CalculateFlux.jl because that 
computation is relatively quick.

Plots are saved as .pdfs (for the write-up) and .pngs (for the presentation)
=#

# Packages and utilities
using StatGeochem
using DelimitedFiles
using StatsBase
using Plots
using CairoMakie
using GeoMakie
using ImageMagick           # Needed to save (Geo?)Makie plots
using MAT

include("Utilities.jl")     # Modify depending on which directory you're in


## --- Parse bulk EarthChem data [CalculateFlux.jl duplicate]
    # Read file
    bulk_file = matopen("C:/Users/rowan/Desktop/Research/Weathering/Bulk_EarthChem/bulk.mat")
    bulk_dict = read(bulk_file, "bulk")
    close(bulk_file)

    bulk = NamedTuple{Tuple(Symbol.(keys(bulk_dict)))}(values(bulk_dict))

    bulksed, bulkign, bulkmet = match_earthchem(bulk.Type)    # met excludes metaseds and metaigns
    bulkcryst = bulkign .| bulkmet                            # Crystalline (ign and met)
    bulk_not_matched = .!bulksed .& .!bulkign .& .!bulkmet    # Unmatched samples
    bulk_matched = .!bulk_not_matched                         # All matched samples

    # Calculate wt.% as only P2O5
    pwt_sed = nanmean(bulk.P2O5[bulksed])
    pwt_ign = nanmean(bulk.P2O5[bulkign])
    pwt_met = nanmean(bulk.P2O5[bulkmet])
    pwt_cryst = nanmean(bulk.P2O5[bulkcryst])


## --- Get slope and erosion rate data
ms_parsed = importdataset("output/responses_parsed.tsv", '\t', importas=:Tuple)

# Parse sed, met, ign, and cover into BitVectors
npoints = length(ms_parsed.sed)
sed = Bool.(ms_parsed.sed)
ign = Bool.(ms_parsed.ign)
met = Bool.(ms_parsed.met)
cover = Bool.(ms_parsed.cover)

# Get cryst and matched BitVectors
cryst = ign .| met
matched = sed .| ign .| met .| cover

# Lat, lon, and slope
rocklat = ms_parsed.lat
rocklon = ms_parsed.lon
rockslope = ms_parsed.slope

# Recalculate the erosion rate
rock_ersn = emmkyr.(rockslope)


## --- Calculate what % of total eroded P comes from each rock type [CalculateFlux.jl duplicate]
    # Constants
    const contl_area = 148940000 * 1000000    # Area of continents (m²)
    const crustal_density = 2750              # Average crustal density (kg/m³)

    # Fraction exposed crust that is each rock type - from our randomly selected samples 
    known = sed .| ign .| met   # Note that this excludes cover
    total_known = count(known)

    sed_area_frac = count(sed) / total_known
    ign_area_frac = count(ign) / total_known
    met_area_frac = count(met) / total_known
    cryst_area_frac = count(cryst) / total_known

    # Average erosion rates for each rock type
    ersn_sed_m = nanmean(rock_ersn[sed])
    ersn_ign_m = nanmean(rock_ersn[ign])
    ersn_met_m = nanmean(rock_ersn[met])
    ersn_cryst_m = nanmean(rock_ersn[cryst])

    # Crustal areas (m^2)
    sed_area = contl_area * sed_area_frac
    ign_area = contl_area * ign_area_frac
    met_area = contl_area * met_area_frac
    cryst_area = contl_area * cryst_area_frac

    # Calculate P flux by source contributions in kg/yr
    sed_contrib = ersn_sed_m * pwt_sed/100 * sed_area * crustal_density * 1e-6
    ign_contrib = ersn_ign_m * pwt_ign/100 * ign_area * crustal_density * 1e-6
    met_contrib = ersn_met_m * pwt_met/100 * met_area * crustal_density * 1e-6
    cryst_contrib = ersn_cryst_m * pwt_cryst/100 * cryst_area * crustal_density * 1e-6
    global_contrib = sed_contrib + ign_contrib + met_contrib


## --- Plot average erosion rate by rock type
    labels = ["$(round(pwt_sed, sigdigits=3)) wt.%", 
        "$(round(pwt_ign, sigdigits=3)) wt.%", 
        "$(round(pwt_met, sigdigits=3)) wt.%"]

    h = CairoMakie.barplot([1,2,3], [pwt_sed, pwt_ign, pwt_met], bar_labels = labels,
        flip_labels_at = 0 , label_color=:white, color = [:teal, :crimson, :darkorange], 
        axis = (ylabel = "Average Phosphorus [wt.%]", xticks = (1:3, ["Sedimentary", "Igneous", "Metamorphic"]),
        ylabelsize=25, xticklabelsize=25, yticklabelsize=20), label_size = 25,
    )

    display(h)
    save("output/figures_flux/pwt.pdf", h)
    save("output/figures_flux/pwt.png", h)


## -- Plot average phosphorous wt.% by rock type
    labels = ["$(round(ersn_sed_m, sigdigits=4)) m/Myr", 
        "$(round(ersn_ign_m, sigdigits=4)) m/Myr", 
        "$(round(ersn_met_m, sigdigits=4)) m/Myr"]

    h = CairoMakie.barplot([1,2,3], [ersn_sed_m, ersn_ign_m, ersn_met_m], bar_labels = labels,
    flip_labels_at = 0 , label_color=:white, color = [:teal, :crimson, :darkorange], 
    axis = (ylabel = "Average Erosion Rate [m/Myr]", xticks = (1:3, ["Sedimentary", "Igneous", "Metamorphic"]),
    ylabelsize=25, xticklabelsize=25, yticklabelsize=20), label_size = 25,
    )

    display(h)
    save("output/figures_flux/ersn-rate_avg.pdf", h)
    save("output/figures_flux/ersn-rate_avg.png", h)

## --- Plot crustal composition (area) by rock type
    labels = ["$(round(sed_area/contl_area*100, sigdigits=3))%", 
        "$(round(ign_area/contl_area*100, sigdigits=3))%", 
        "$(round(met_area/contl_area*100, sigdigits=3))%"]

    h = CairoMakie.barplot([1,2,3], [sed_area, ign_area, met_area], bar_labels = labels,
        flip_labels_at = 0 , label_color=:white, color = [:teal, :crimson, :darkorange], 
        axis = (ylabel = "Surface Composition [m²]", xticks = (1:3, ["Sedimentary", "Igneous", "Metamorphic"]),
        ylabelsize=25, xticklabelsize=25, yticklabelsize=20), label_size = 25,
    )

    display(h)
    save("output/figures_flux/crust_comp.pdf", h)
    save("output/figures_flux/crust_comp.png", h)


## --- Plot phosphorous flux contribution by rock type
    labels = ["$(round(sed_contrib/global_contrib*100, sigdigits=3))%", 
        "$(round(ign_contrib/global_contrib*100, sigdigits=3))%", 
        "$(round(met_contrib/global_contrib*100, sigdigits=3))%"]
    
    h = CairoMakie.barplot([1,2,3], [sed_contrib, ign_contrib, met_contrib], bar_labels = labels,
        flip_labels_at = 0 , label_color=:white, color = [:teal, :crimson, :darkorange], 
        axis = (ylabel = "P Flux [kg/year]", xticks = (1:3, ["Sedimentary", "Igneous", "Metamorphic"]),
        ylabelsize=25, xticklabelsize=25, yticklabelsize=20), label_size = 25,
    )
    display(h)
    save("output/figures_flux/pflux.pdf", h)
    save("output/figures_flux/pflux.png", h)


## --- Make a function! Everyone loves functions!
"""
```julia
global_visual_ersn(rocklon, rocklat, rock_ersn, name::String="all")
```
Global visualization of erosion rate `rock_ersn`. Saves .pdf and .png files as 
"ersnrate_`name`".
"""
function global_visual_ersn(rocklon, rocklat, rock_ersn, name::String="all")
    f = Figure(resolution = (1200, 600))
    ax = GeoAxis(f[1,1]; coastlines = true, dest = "+proj=wintri")
    h = CairoMakie.scatter!(ax, rocklon, rocklat, markersize = 5, 
        color = log.(rock_ersn), colormap=:thermal
    )
    Colorbar(f[1,2], h, label = "log Erosion Rate [log m/Myr]", height = Relative(0.9))
    display(f)
    save("output/figures_flux/ersnrate_$name.pdf", f)
    save("output/figures_flux/ersnrate_$name.png", f)
end


## --- Mass produce some erosion rate visualizations
    global_visual_ersn(rocklon, rocklat, rock_ersn, "all")  # Also a map of slope!
    global_visual_ersn(rocklon[sed], rocklat[sed], rock_ersn[sed], "sed")
    global_visual_ersn(rocklon[ign], rocklat[ign], rock_ersn[ign], "ign")
    global_visual_ersn(rocklon[met], rocklat[met], rock_ersn[met], "met")
    global_visual_ersn(rocklon[cover], rocklat[cover], rock_ersn[cover], "cover")


## --- Histograms of erosion rate and erosive flux
"""
```julia
plots_log_hist(data, barcolor::Symbol;
    type::Symbol=:rate,
    name::String, 
    anname::String="", 
    printxaxis::Bool=true
)
```
Returns and saves a `barcolor` colored log-scale histogram for a data set `data`. 
    
Keyword arguments and defaults:

    type

Dependent variable. `:rate` for erosion rate, `:flux` for erosive flux, or `:wt` for wt.%.
Default is `:rate`.

    name

Identifying string to append to saved file name.

    annname

String to annotate plot with. Defaults to empty.

    printxaxis

Include the x-axis label and ticks. Defaults to true.
"""
function plots_log_hist(data, barcolor::Symbol;
        type::Symbol=:rate,
        name::String, 
        anname::String="", 
        printxaxis::Bool=true
    )

    # Make plot of choice
    if type==:rate
        mylogbinedges = range(start=0.5, stop=6, length=26) # 25 bins
        counts = histcounts(log10.(data), mylogbinedges)    # Get count of each log-spaced bin
        filter = counts .> 0                                # Filter out bins with no data

        h = Plots.plot(cntr(mylogbinedges)[filter], (counts[filter]), seriestype=:bar, yaxis=:log,
            color=barcolor, linecolor=barcolor, xlims=(0.7, 5.5), 
            framestyle=:box, label="", ylabel="log Sample Count", tickfontsize=15, labelfontsize=15, 
            y_foreground_color_text=barcolor, y_guidefontcolor=barcolor
        )

        # Coordinate positions for annotation
        ypos = ylims(h)[2] - 0.5*ylims(h)[2]
        xpos = 3.5

    elseif type==:flux
        mylogbinedges = range(start=0.5, stop=6, length=26)           # 25 bins
        rates = cntr(mylogbinedges)                                   # Erosion rate for each bin (m/Myr)
        counts = histcounts(log10.(data), mylogbinedges)              # Get count of each log-spaced bin
        fluxes = counts / total_known * contl_area .* exp10.(rates)   # Erosive flux (m³/Myr)
        filter = fluxes .> 0                                          # Filter out bins with no data

        h = Plots.plot(cntr(mylogbinedges)[filter], (fluxes[filter]), seriestype=:bar, yaxis=:log,
            color=barcolor, linecolor=barcolor, xlims=(0.7, 5.5), ylims=(9e12, 12e14), fillrange=10e12,
            framestyle=:box, label="", ylabel="log Erosive Flux [m³/Myr]", tickfontsize=15, labelfontsize=15, 
            y_foreground_color_text=barcolor, y_guidefontcolor=barcolor
        )

        # Coordinate positions for annotation
        ypos = ylims(h)[2] - 0.25*ylims(h)[2]
        xpos = 3.5

    elseif type==:wt 
        mylogbinedges = range(start=0.01, stop=2, length=26)    # 25 bins
        counts = histcounts(log10.(data), mylogbinedges)        # Get count of each log-spaced bin
        filter = counts .> 0                                    # Filter out bins with no data

        h = Plots.plot(cntr(mylogbinedges)[filter], (counts[filter]), seriestype=:bar, yaxis=:log,
            color=barcolor, linecolor=barcolor, xlims=(-0.08,1.9),
            framestyle=:box, label="", ylabel="log Sample Count", tickfontsize=15, labelfontsize=15, 
            y_foreground_color_text=barcolor, y_guidefontcolor=barcolor
        )

        # Semi-log plot option. Keeping this because it's interesting but it's also easier to
        # see some of the differences between rock types in log-log space.

        # mylinbinedges = range(start = 0, stop=55, length=51)      # 50 bins, linspaced
        # counts = histcounts((data), mylinbinedges)                # count, linspaced

        # h = Plots.plot(cntr(mylinbinedges)[filter], (counts[filter]), seriestype=:bar, yaxis=:log,
        #     color=barcolor, linecolor=barcolor,
        #     framestyle=:box, label="", ylabel="log Sample Count", tickfontsize=15, labelfontsize=15, 
        #     y_foreground_color_text=barcolor, y_guidefontcolor=barcolor
        # )

        # Coordinate positions for annotation
        ypos = ylims(h)[2] - 0.4*ylims(h)[2]
        xpos = 1.6

    end

    # Annotate plot with the rock type
    Plots.annotate!(xpos, ypos, anname, annotationcolor=barcolor, annotationhalign=:left,
        annotationfontsize=18, 
    )   

    # Print the x-axis labels, or don't
    # x-ticks have to be done manually because Plots doesn't know we're using a log scale
    if !printxaxis 
        Plots.plot!(xticks=nothing)
    elseif type!=:wt
        Plots.plot!(xlabel="log Erosion Rate [m/Myr]", xticks=(1:5, ["10⁰", "10²", "10³", "10⁴", "10⁵"]))
    else
        Plots.plot!(xlabel="log P₂O₅ Content [wt.%]", xticks=(0:1, ["0", "10¹"]))
    end

    # Graphic design is my passion 🐸
    savefig(h, "output/figures_flux/log-ersn_$(string(type))_$name.pdf")
    savefig(h, "output/figures_flux/log-ersn_$(string(type))_$name.png")

    return h
end


## --- Generate histograms of erosion rate and erosive flux
    hsr = plots_log_hist(rock_ersn[sed], :teal, name="sed", anname="Sedimentary", printxaxis=false, type=:rate)
    hsf = plots_log_hist(rock_ersn[sed], :teal, name="sed", anname="Sedimentary", printxaxis=true, type=:flux)

    hir = plots_log_hist(rock_ersn[ign], :crimson, name="ign", anname="Igneous", printxaxis=false, type=:rate)
    hif = plots_log_hist(rock_ersn[ign], :crimson, name="ign", anname="Igneous", printxaxis=true, type=:flux)

    hmr = plots_log_hist(rock_ersn[met], :darkorange, name="met", anname="Metamorphic", printxaxis=false, type=:rate)
    hmf = plots_log_hist(rock_ersn[met], :darkorange, name="met", anname="Metamorphic", printxaxis=true, type=:flux)

    # # Option to make two plots, one for rate and one for flux
    # h = Plots.plot(hsr, hir, hmr, layout = (3, 1), size = (700,1200), left_margin=40Plots.px)
    # display(h)
    # savefig(h, "output/figures_flux/log-ersn_rate_combined.pdf")
    # savefig(h, "output/figures_flux/log-ersn_rate_combined.png")

    # h = Plots.plot(hsf, hif, hmf, layout = (3, 1), size = (700,1200), left_margin=40Plots.px)
    # display(h)
    # savefig(h, "output/figures_flux/log-ersn_flux_combined.pdf")
    # savefig(h, "output/figures_flux/log-ersn_flux_combined.png")

    # Option 2 for a even bigger plot with a different layout
    h = Plots.plot(hsr, hir, hmr, hsf, hif, hmf, layout = (2, 3), size = (2100,900), 
        left_margin=50Plots.px, bottom_margin=50Plots.px)
    display(h)
    savefig(h, "output/figures_flux/ersn_all.pdf")


## --- Plot log histogram distribution of phosphorus by rock type
    hsw = plots_log_hist(bulk.P2O5[bulksed], :teal, name="sed", anname="Sed.", printxaxis=false, type=:wt)
    hiw = plots_log_hist(bulk.P2O5[bulkign], :crimson, name="ign", anname="Ign.", printxaxis=false, type=:wt)
    hmw = plots_log_hist(bulk.P2O5[bulkmet], :darkorange, name="met", anname="Met.", printxaxis=true, type=:wt)

    h = Plots.plot(hsw, hiw, hmw, layout = (3, 1), size = (700,1200), left_margin=40Plots.px)
    display(h)
    savefig(h, "output/figures_flux/pwt_combined.pdf")
    savefig(h, "output/figures_flux/pwt_combined.png")
    

## --- Grab some rock type data for the entire Earth
    hartmann = readdlm("C:/Users/rowan/Desktop/Research/hartmann-moosdorf_2012/hartmann_2012.txt.asc", skipstart=6)
        
    # Use a coarser rock type than Hartmann does:
    # seds = 1, mets = 2, igns = 3, cover = 4, ice = 5, water / no data = NaN

    # I know this is a little cursed, but this is much easier to debug than a loop. Offsets
    # of 16 so nothing accidentaly gets reassigned and Australia (hypothetically) ends up
    # mostly igneous.
    hartmann[hartmann .== -9999] .= NaN    # Missing data (ocean)

    hartmann[hartmann .== 1]  .= 4 + 16    # Unconsolidated sediment   
    hartmann[hartmann .== 2]  .= 3 + 16    # Basic volcanics
    hartmann[hartmann .== 3]  .= 1 + 16    # Siliciclastic seds
    hartmann[hartmann .== 4]  .= 3 + 16    # Basic plutonics
    hartmann[hartmann .== 5]  .= 1 + 16    # Mixed seds
    hartmann[hartmann .== 6]  .= 1 + 16    # Carbonate seds
    hartmann[hartmann .== 7]  .= 3 + 16    # Acid volcanics
    hartmann[hartmann .== 8]  .= 2 + 16    # Metamorphics
    hartmann[hartmann .== 9]  .= 1 + 16    # Acid plutonics
    hartmann[hartmann .== 10] .= 1 + 16    # Intermediate volcanics
    hartmann[hartmann .== 11] .= NaN       # Water
    hartmann[hartmann .== 12] .= 3 + 16    # Pyroclastics
    hartmann[hartmann .== 13] .= 3 + 16    # Intermediate plutonics
    hartmann[hartmann .== 14] .= 1 + 16    # Evaporites
    hartmann[hartmann .== 15] .= NaN       # No data
    hartmann[hartmann .== 16] .= 5 + 16    # Ice and glaciers


## --- Geospatial rock type visualization
    # All rocks
    h = Plots.heatmap(reverse(hartmann, dims=1), size=(1200,600), framestyle=:box, 
        legend=nothing, xticks=nothing, yticks=nothing,
        c = [:teal, :darkorange, :crimson, :khaki, :lightcyan]      # The right colors are very important
    )
    display(h)
    savefig(h, "output/figures_flux/rocktype_all.pdf")
    savefig(h, "output/figures_flux/rocktype_all.png")

    # Seds
    sed_hm = copy(hartmann)
    sed_hm[sed_hm.!= 1 + 16] .= NaN
    h = Plots.heatmap(reverse(sed_hm, dims=1), size=(1200,600), framestyle=:box, c = :teal,
        legend=nothing, xticks=nothing, yticks=nothing,
    )
    display(h)
    savefig(h, "output/figures_flux/rocktype_sed.pdf")

    # Igns 
    ign_hm = copy(hartmann)
    ign_hm[ign_hm .!= 3 + 16] .= NaN
    h = Plots.heatmap(reverse(ign_hm, dims=1), size=(1200,600), framestyle=:box, c = :crimson,
        legend=nothing, xticks=nothing, yticks=nothing,
    )
    display(h)
    savefig(h, "output/figures_flux/rocktype_ign.pdf")

    # Mets 
    met_hm = copy(hartmann)
    met_hm[met_hm .!= 2 + 16] .= NaN
    h = Plots.heatmap(reverse(met_hm, dims=1), size=(1200,600), framestyle=:box, c = :darkorange,
        legend=nothing, xticks=nothing, yticks=nothing,
    )
    display(h)
    savefig(h, "output/figures_flux/rocktype_met.pdf")

    # Calculate relative abundance of rock type from this dataset
    hartmann_sed = nansum(sed_hm ./ 17)
    hartmann_ign = nansum(ign_hm ./ 19)
    hartmann_met = nansum(met_hm ./ 18)
    hartmann_total = hartmann_sed + hartmann_ign + hartmann_met

    hartmann_sed_rel = hartmann_sed / hartmann_total
    hartmann_ign_rel = hartmann_ign / hartmann_total
    hartmann_met_rel = hartmann_met / hartmann_total

## --- End of File