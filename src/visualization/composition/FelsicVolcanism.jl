## --- Set up
    # Why are igneous rocks more felsic at ~ 2 Ga and ~500 Ma?

    # Packages 
    using RockWeatheringFlux
    using StatsBase
    using DelimitedFiles, HDF5
    using Plots 


## --- Load data 
    # Heatmap
    suffix = RockWeatheringFlux.version * "_" * RockWeatheringFlux.tag
    fpath = "src/visualization/composition/shortcuts/SilicaAgeDistribution_SiO2_" * suffix * ".h5"

    target = (:ign, :plut, :volc)
    fid = h5open(fpath, "r")
    out_mbulk = NamedTuple{target}(read(fid["vars"]["matched"]["$key"]) for key in target)
    close(fid)

    # Matched geochemical samples 
    fid = readdlm(matchedbulk_io)
    gchem_ind = Int.(vec(fid[:,1]))
    t = @. gchem_ind != 0
    match_cats = get_lithologic_class()[1];

    fid = h5open(geochem_fid, "r")
    header = read(fid["bulk"]["header"])
    data = read(fid["bulk"]["data"])
    mbulk = NamedTuple{Tuple(Symbol.(header))}([data[:,i][gchem_ind[t]] for i in eachindex(header)])
    close(fid)
    
    # Macrostrat 
    fid = h5open("$macrostrat_io", "r")
    macrostrat = (
        rocklat = read(fid["vars"]["rocklat"])[t],
        rocklon = read(fid["vars"]["rocklon"])[t],
        age = read(fid["vars"]["age"])[t],
        agemax = read(fid["vars"]["agemax"])[t],
        agemin = read(fid["vars"]["agemin"])[t],
        scale = read(fid["vars"]["scale"])[t],
        rocktype = read(fid["vars"]["rocktype"])[t],
        rockname = read(fid["vars"]["rockname"])[t],
        rockdescrip = read(fid["vars"]["rockdescrip"])[t],
        reference = read(fid["vars"]["reference"])[t],
    )
    header = read(fid["type"]["macro_cats_head"])
    data = read(fid["type"]["macro_cats"])
    data = @. data > 0
    macro_cats = NamedTuple{Tuple(Symbol.(header))}([data[:,i][t] for i in eachindex(header)])
    close(fid)


## --- Pull heatmap: (scale = relative abundance)
    # Base heatmap for this target 
    h = plot(
        yflip=true,
        size = (600, 200),
        # colorbar=false, 
        grid=false,
        framestyle=:box,
    )

    # Set target and refine by looking at heatmap 
    ymin, ymax = 250, 550                    # Target 1 age 
    # ymin, ymax = 1500, 2500                  # Target 2 age 
    ybins = size(out_mbulk.volc)[1]
    xmin, xmax = 40, 80                     # SiO2
    xbins = size(out_mbulk.volc)[2]
    
    imin = Int(rescale_in_range(ymin, 0, 3800, 0, ybins))
    imax = Int(rescale_in_range(ymax, 0, 3800, 0, ybins))

    hsim = deepcopy(h);
    x = xmin:10:xmax 
    y = round.(Int, range(start=ymin, stop=ymax, length=5))
    heatmap!(hsim, out_mbulk.volc[imin:imax, :],
        xlims=(0.5,xbins+0.5),
        ylims=(0.5,(imax-imin)+0.5), 
        xticks=(rescale_in_range.(x, xmin, xmax, 0.5, xbins+0.5), x),
        yticks=(rescale_in_range.(y, ymin, ymax, 0.5, (imax-imin)+0.5), string.(y)),
    )


## --- Pull base data (scale = sample counts)
    # Since we don't resample, we need way coarser bins 
    ybinwidth = 10
    ybins = Int((ymax - ymin) / ybinwidth)
    yedges = ymin:(ymax-ymin)/ybins:ymax

    xbinwidth = 0.5
    xbins = Int((xmax - xmin) / xbinwidth)
    xedges = xmin:(xmax-xmin)/xbins:xmax

    imin = Int(rescale_in_range(ymin, ymin, ymax, 0, ybins))
    imax = Int(rescale_in_range(ymax, ymin, ymax, 0, ybins))

    intime = @. ymin .<= mbulk.Age .<= ymax;
    intarget = intime .& match_cats.volc
    out = zeros(ybins, xbins)
    for j = 1:ybins
        s = @. yedges[j] <= mbulk.Age[intarget] < yedges[j+1]
        c, n = bincounts(mbulk.SiO2[intarget][s], xmin, xmax, xbins)
        out[j,:] .= n
    end

    hobs = deepcopy(h);
    heatmap!(hobs, out,
        xlims=(0.5,xbins+0.5),
        ylims=(0.5,(imax-imin)+0.5), 
        xticks=(rescale_in_range.(x, xmin, xmax, 0.5, xbins+0.5), x),
        yticks=(rescale_in_range.(y, ymin, ymax, 0.5, (imax-imin)+0.5), string.(y)),
    )

    # Target 1:
    # It almost looks like it's only a few samples that are driving the felsic shift. 
    # Like something at ~400 Ma and ~72% silica has ~200 samples in it.... 

    # Target 2:
    # This is a little more spread out 


## --- Pull weird samples 
    # Find all the indices that are above the 99.9th percentile for sample counts 
    p = percentile(vec(out), 99.9)
    i = findall(>(p), out)

    # Translate to ages and SiO2 
    ages = zeros(length(i))
    silica = zeros(length(i))
    for j in eachindex(i)
        ages[j] = rescale_in_range(i[j][1], 0, ybins, ymin, ymax);
        silica[j] = rescale_in_range(i[j][2], 0, xbins, xmin, xmax);
        println("$(ages[j]) -> $(silica[j]) % (Ma, wt.%  SiO2) $(Int(out[i[j]])) samples")
    end


## --- [THIS CODE IS MEANT TO BE RUN LINE BY LINE]
    # Note for both targets we have a high silica mode and a low silica mode. Do some 
    # manual selection to just look at what's unexpected  
    # Grab weird shit

    # Target 1: two grid squares that we care about -- 400 Ma and 410 Ma, 71.5% SiO2
    target_age = @. (400 - ybinwidth) <= mbulk.Age < (410 + ybinwidth);
    target_SiO2 = @. (71.5 - xbinwidth) <= mbulk.SiO2 < (71.5 + xbinwidth);
    target = target_age .& target_SiO2 .& match_cats.volc;
    count(target)   # 341


    # # Target 2: 4-5 samples. For this, let's automate since they're not next to each other 
    # t = silica .> 70;
    # target = falses(length(mbulk.SiO2));
    # for j in 1:count(t)
    #     SiO2 = silica[t][j]
    #     age = ages[t][j]
    #     target .|= (
    #         (age-ybinwidth .<= mbulk.Age .< age+ybinwidth) .&
    #         (SiO2-xbinwidth .<= mbulk.SiO2 .< SiO2+xbinwidth) .&
    #         (match_cats.volc)
    #     )
    # end
    # count(target)   # 904

    # Look at the rocks we have (restrict to just those that represent a certain number
    # of samples, although you may not want to do this if you don't have a lot)
    c = countmap(macrostrat.rocktype[target])
    c = sort(collect(c), by=x->x[2], rev=true);
    p = percentile(last.(c), 50);
    q = last.(c) .> p;
    rocktype = hcat(last.(c)[q], first.(c)[q])  

    # Target 1: 267 of these 341 samples are "intermediate-felsic volcanic rocks"
        # damn.
    # Target 2: 519 of 904 samples
        # 210 are "sedimentary and volcaniclastic"
        # 195 are "medium-high grade orthogneiss/paragneiss/metavolcanic gneiss"
        # 111 are "metamorphic and undivided crystalline: orthogneiss"


## --- Look at the samples that are matching with those most-common types 
    # (everything that matches, not just those in our target)
    # That way we can see if there's any weird biases in how these super uninformative 
    # names are being treated

    # Maybe we want to grab a couple different names (or maybe we don't)
    r = falses(length(macrostrat.rocktype));
    # i = 1:3       # target 2
    i = 1:1         # target 1
    for j in i 
        println(rocktype[j,2])
        r .|= macrostrat.rocktype .== rocktype[j,2]
    end
    count(r)

    countmap(macrostrat.rockname[r])
    countmap(macrostrat.rockdescrip[r])
    countmap(macrostrat.age[r])

    # For both targets, we have NO useful information that we could use to get information 
    # from these rocks. 
    # Target 1 is just listed as volcanic and sedimentary rocks

    # Target 2 is just more of the same, but only for the undiff. metamorphics. A lot of them are 
    # listed as orthogneiss which I do belive we route to igneous, so there's that. But no actual info 
    # about composition. Empty strings in the rockdescrip
    # OOOhhhh so this is a pulse of orthogneiss maybe?

    # Where are they from?
    countmap(macrostrat.scale[r])
    countmap(macrostrat.reference[r])

    # Target 1: tiny resolution, generalized geologic map of the world
    # Target 2: small or tiny resolution, smattering of maps 

    mapplot(macrostrat.rocklon[r], macrostrat.rocklat[r], markersize=1, msw=0, label="")

    # What are they matching with?  
    histogram(mbulk.SiO2[r .& match_cats.volc])                
    unique_sample(mbulk.Sample_ID[r .& match_cats.volc], 90)   # Lots of different samples in here, no favorites

    # Target 1 is a mix of andesite and rhyolite        
    # Target 2 is a good bimodal but more mafic distribution 
    

## --- End of file 