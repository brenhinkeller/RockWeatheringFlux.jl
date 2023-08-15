# Functions for data analysis and meta-analysis

## --- Modal bin / sample functions

    """
    ```julia
    modal(data::AbstractArray)
    ```

    Find the mode and the frequency of the mode in `data`.

    # Example
    ```julia-repl
    julia> data = [1, 1, 1, 3, 4];

    julia> modal(data)
    (1, 3)
    ```
    """
    function modal(data::AbstractArray)
        values = unique(data)
        counts = [count(==(i), data) for i in values]

        f = findmax(counts)[1]              # Frequency of modal value
        v = values[findmax(counts)[2]]      # Modal value

        return v, f
    end

    """
    ```julia
    sameindex(type::Symbol, macro_cats, bulk, bulkidx; hist=:on)
    ```

    Count how many samples in the modal bin are from the same EarthChem sample. Optionally
    plot a histogram of the frequency of all selected indices. Print information to terminal
    and return the index of the EarthChem sample.

    ### Requirements and Assumptions
      * Data in `bulk` are unmatched samples in the original order
      * Indices in `bulkidx` are filtered so there are no indices of 0

    # Example
    ```julia-repl
    julia> iᵢ = sameindex(:ign, macro_cats, bulk, bulktext, bulkidx[t], bins=ignbin, hist=:off);
    ┌ Info: Type: ign
    │ Modal bin: 11 (50.0-51.0 wt.% SiO₂)
    │ Modal index count: 356 of 3746 (9.5%)
    │ Index: 413791
    │ 
    │ This sample represents 0.88% of all ign matches.
    │ ---
    │ Sample information:
    │ 
    │ Age: NaN
    │ Lat, Lon: 53.94, -148.53
    │ SiO₂: 50.98%
    │ 
    │ Name: basalt
    │ Type: volcanic
    └ Material: igneous
    ```
    """
    function sameindex(type::Symbol, macro_cats, bulk, bulktext, bulkidx; 
        bins, hist=:on)

        nzero = @. bulkidx == 0
        @assert count(nzero) == 0 "Filter bulkidx to remove indices of 0"

        filter = macro_cats[type]

        c, n, = bincounts(bulk.SiO2[bulkidx][filter], bins...,)

        # Find modal bin and all the points in it
        i = findmax(n)[2]
        s = step(c)/2
        tᵢ = @. c[i]-s <= bulk.SiO2[bulkidx][filter] <= c[i]+s

        # Find indices of those points, and get modal data
        ind = bulkidx[filter][tᵢ]
        j, f = modal(ind)

        # What percent of total indices are the modal index?
        totalcount = count(==(j), bulkidx[filter])
        totalindex = length(bulkidx[filter])

        # Terminal printout
        @info """
        Type: $type
        Modal bin: $i ($(c[i]-s)-$(c[i]+s) wt.% SiO₂)
        Modal index count: $f of $(length(ind)) ($(round(f/length(ind)*100, sigdigits=3))%)
        Index: $(j)

        This sample represents $(round(totalcount/totalindex*100, digits=3))% of all $type matches.
        ---
        Sample information:

        Age: $(bulk.Age[j])
        Lat, Lon: $(bulk.Latitude[j]), $(bulk.Longitude[j])
        SiO₂: $(round(bulk.SiO2[j], sigdigits=3))%

        Name: $(bulktext.Rock_Name[j])
        Type: $(bulktext.Type[j])
        Material: $(bulktext.Material[j])
        """

        # Histogram
        if hist==:on
            c = sort!(unique(bulkidx[filter]))
            n = [count(==(i), bulkidx[filter]) for i in c]
            h = plot(c, n, seriestype=:bar, label="$type", xlabel="Index [all bins]", 
                ylabel="Frequency", framestyle=:box, ylims=(0, maximum(n)+100), color=:black,
            )
            display(h)
        end

        return j
    end


    """
    ```julia
    get_matched_samples(i::Int64, bulkidx::AbstractArray{<:Int64}, macrostrat::NamedTuple;
        [savefile::Bool], [filepath::String], 
        [desc::String],
        [filter::BitArray]
    )
    ```

    Get Macrostrat sample information for all samples matched with index `i`. 
        
    ### Optional kwargs:
      * `savefile`: save a .csv file to `filepath` with the data for all matched samples.
        Boolean; `false` by default.
      * `filter`: filter `macrostrat` samples.
      * `desc`: Additional information to identify the terminal printout.
    """
    function get_matched_samples(i::Int64, bulkidx::AbstractArray{<:Int64}, macrostrat::NamedTuple;
        savefile::Bool=false, filepath::String="", desc::String="\b", 
        filter::BitArray{1}=trues(length(bulkidx))
    )

        # Check inputs
        @assert i > 0 "Index i must be greater than 0."
        @assert length(bulkidx) == length(macrostrat[1]) == length(filter) """
            Lengths of bulkidx, filter, and Macrostrat sample arrays must be equal."""
        if savefile
            @assert filepath != "" "Specify file path to save file."
            @assert containsi(filepath, ".csv") "Filepath must specify .csv extension."
        end
        
        required = ifelse(savefile, (:rocktype, :rockname, :rockdescrip, :rocklat, :rocklon, :age),
            (:rocktype, :rockname, :rockdescrip, :age)
        )
        present = keys(macrostrat)
        for k in required
            @assert k in present "Macrostrat key $k required."
        end

        # Get matched samples
        t = @. bulkidx == i && bulkidx != 0 && filter
        @assert count(t) > 0 "No samples found"

        rocktype = macrostrat.rocktype[t]
        rockname = macrostrat.rockname[t]
        rockdescrip = macrostrat.rockdescrip[t]
        age = macrostrat.age[t]
        if savefile
            rocklat = macrostrat.rocklat[t]
            rocklon = macrostrat.rocklon[t]
        end

        # Get modal data
        agem = (modal(age)..., length(age))
        typem = (modal(rocktype)..., length(rocktype))
        namem = (modal(rockname)..., length(rockname))
        descm = (modal(rockdescrip)..., length(rockdescrip))

        # Print to terminal
        @info """
        $desc Macrostrat samples matched with index $i 
        
        Age [Ma]:
        Mean: $(round(nanmean(age), digits=2))
        Standard deviation: $(round(nanstd(age), digits=2))
        Median: $(nanmedian(age))
        Mode: $(agem[1]) (n = $(agem[2]) of $(agem[3]))
        ---

        Modal sample description:
        Rock type: $(typem[1]) (n = $(typem[2]) of $(typem[3]))
        Rock name: $(namem[1]) (n = $(namem[2]) of $(namem[3]))
        Rock description: $(descm[1]) (n = $(descm[2]) of $(descm[3]))
        """

        # Save to file, if directed
        if savefile
            header = ["age" "lat" "lon" "type" "name" "desc"]
            writedlm("$filepath", vcat(header, hcat(age, rocklat, rocklon, rocktype, rockname,
                rockdescrip))
            )
        end
    end


## --- End of file