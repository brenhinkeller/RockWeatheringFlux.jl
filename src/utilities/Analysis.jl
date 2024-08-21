# Functions for data analysis and meta-analysis

## --- Format an array as a LaTeX table 
    """
    ```julia
    format_latex(A)
    ```

    Format `A` as a LaTeX table, to be copy-pasted into a LaTeX document.
    """
    function format_latex(A::AbstractArray)
        row = size(A)[1]
        for r = 1:row 
            println(join(A[r,:], " & ") * " \\" * "\\")
        end
    end
    export format_latex
    

## --- Visualize rock class tuples
    """
    ```julia
    terminal_cats(cats1, cats2; [catkeys, name1, name2])
    ```

    Print the number of `true` values in each key of `cat1` and `cat2`. Optionally specify
    specific keys to print as `catkeys`. 
    """
    function terminal_cats(cats1::NamedTuple, cats2::NamedTuple; 
            catkeys::AbstractArray{Symbol} = collect(keys(cats1)), 
            name1::String="cats1", name2::String="cats2"
        )
        @assert keys(cats1) == keys(cats2)
        # catkeys = collect(keys(cats1))
        allkeys = collect(keys(cats1))

        # Get pad lengths
        keymax = length(string(allkeys[sortperm(string.(allkeys), by = length)[end]]))
        max1 = max(length(name1), 6)
        max2 = max(length(name2), 6)

        # Terminal printout
        println("$(rpad("Rock Class", keymax)) \t $(lpad(name1, max1)) $(lpad(name2, max1)) \t Maximum")
        for k in catkeys
            c1 = count(cats1[k])
            c2 = count(cats2[k])

            larger = ifelse(c1>c2, rpad(name1, max1), rpad(name2, max1))
            diff = max(c1, c2) - min(c1, c2)

            println("$(rpad(k, keymax)) \t $(lpad(c1, max1)) $(lpad(c2, max1)) \t $larger (+$diff)")
        end
    end

    # I need this code to work, but if I don't have a kwarg in the previous method, this
    # overwrites it and I don't understand why :(
    # function terminal_cats(cats1::NamedTuple, cats2::NamedTuple; catkeys::AbstractArray{Symbol})
    #     allkeys = collect(keys(cats1))
    #     strmax = length(string(allkeys[sortperm(string.(allkeys), by = length)[end]]))

    #     for k in catkeys
    #         c1 = count(cats1[k])
    #         c2 = count(cats2[k])
    #         println("$(rpad(k, strmax)) \t $(lpad(c1, 6)) $(lpad(c2, 6))")
    #     end
    # end
    export terminal_cats


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
    export modal

    """
    ```julia
    sameindex(bulkidx::AbstractArray, geochemkeys::AbstractArray{Symbol}, bins::Tuple,
        bulk, bulktext
    )
    ```

    Count how many matched samples are the same EarthChem sample, in both the modal silica
    bin and the entire set of matched indices.

    Count how many samples in the modal bin are from the same EarthChem sample.
    
    Print information to terminal and return the index of the sample.

    ### Arguments

        \tbulkidx

    Matched indices to test.

        \tgeochemkeys

    Elements to be retrived from EarthChem for modal indices.

        \tbins

    Bins (xmin, xmax, nbins) for `bincounts`. Required to determine modal bin.

        \tbulk

    EarthChem samples; as originally loaded for the sample matching algorithm.

        \tbulktext

    EarthChem sample metadata.    

    # Example
    ```julia-repl
    julia> sample_cats = match_rocktype(string.(sampletypes));

    julia> j = sameindex(matches[sample_cats.shale], geochemkeys, (0,100,100), bulk, bulktext)
    ┌ Info: Frequency of the modal sample in:
    │ Modal bin: 28.6%       i = 74106
    │ 
    │ Age: 385.5
    │ Lat, Lon: 69.0, 89.0
    │ 
    │ Name: argillite
    │ Type: 
    │ Material: sedimentary
    │ 
    │ --- --- --- --- ---
    │ All matches: 0.85%     i = 142971
    │ 
    │ Age: 1.305
    │ Lat, Lon: 43.0, 11.0
    │ 
    │ Name: shale
    │ Type: 
    │ Material: sedimentary
    │ 
    │ --- --- --- --- ---
    │ Geochemistry [wt.%] for modal sample in:
    │ 
    │             SiO2       Al2O3   FeOT    TiO2    MgO     CaO     Na2O    K2O     H2O     CO2
    │ Modal bin:   66.0      21.0    4.2     1.5     2.7     0.24    0.89    3.3     0.0     NaN
    └ All samples: 62.0      21.0    6.1     0.87    2.1     4.5     0.75    2.9     0.0     NaN
    ```
    """
    function sameindex(bulkidx::AbstractArray, geochemkeys::AbstractArray{Symbol}, bins::Tuple,
        bulk, bulktext
        )
        # Interpret user input
        t = @. bulkidx != 0
        bulkidx = bulkidx[t]

        # Modal index of modal bin, and of all samples
        c, n, = bincounts(bulk.SiO2[bulkidx], bins...,)     # Get bins
        i = findmax(n)[2]                                   # Index of largest bin
        s = step(c)/2                                       # Dist. bin centers to edges
        tᵢ = @. c[i]-s <= bulk.SiO2[bulkidx] <= c[i]+s      # Filter samples in that bin
        ind = bulkidx[tᵢ]                                   # Indices of the bin
        j, f = modal(ind)                                   # Value, frequency of modal index
        j₁, f₁ = modal(bulkidx)                             # Modal index of all samples

        # Major element composition of modal indices
        majel = round.([bulk[i][j] for i in geochemkeys], sigdigits=2)
        majel₁ = round.([bulk[i][j₁] for i in geochemkeys], sigdigits=2)

        # Terminal printout
        @info """Frequency of the modal sample in:
        Modal bin: $(round(f/length(ind)*100, sigdigits=3))% \t i = $j
        
        Age [Ma]: $(bulk.Age[j])
        Lat, Lon: $(bulk.Latitude[j]), $(bulk.Longitude[j])
        Name:     $(bulktext.Rock_Name[j])
        Type:     $(bulktext.Type[j])
        Material: $(bulktext.Material[j])

        --- --- --- --- ---
        All matches: $(round(f₁/length(bulkidx)*100, sigdigits=3))% \t i = $j₁

        Age [Ma]: $(bulk.Age[j₁])
        Lat, Lon: $(bulk.Latitude[j₁]), $(bulk.Longitude[j₁])
        Name:     $(bulktext.Rock_Name[j₁])
        Type:     $(bulktext.Type[j₁])
        Material: $(bulktext.Material[j₁])

        --- --- --- --- ---
        Geochemistry [wt.%] for modal sample in:

                     $(join(geochemkeys, " \t "))
        Modal bin:   $(join(majel, " \t "))
        All samples: $(join(majel₁, " \t "))
        """

        return j
    end
    export sameindex


    """
    ```julia
    get_matched_samples(i::Int64, bulkidx::AbstractArray{<:Int64}, macrostrat::NamedTuple;
        [desc::String],
        [filter::BitArray]
    )
    ```

    Get Macrostrat sample information for all samples matched with index `i`. 
        
    ### Optional kwargs:
      * `filter`: filter `macrostrat` samples.
      * `desc`: Additional information to identify the terminal printout.

    # Example
    ```julia
    i = sameindex(bulkidx[t][macro_cats.met], majors, (25,100,75), bulk, bulktext);
    get_matched_samples(i, bulkidx[t], macrostrat, filter=macro_cats.met, desc="Metamorphic")
    ```
    """
    function get_matched_samples(i::Int64, bulkidx::AbstractArray{<:Int64}, 
        macrostrat::NamedTuple;
        desc::String="\b", 
        filter::BitArray{1}=trues(length(bulkidx))
    )

        # Check inputs
        @assert i > 0 "Index i must be greater than 0."
        @assert length(bulkidx) == length(macrostrat[1]) == length(filter) """
            Lengths of bulkidx, filter, and Macrostrat sample arrays must be equal."""

        
        required = (:rocktype, :rockname, :rockdescrip, :age)
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

        # Get modal data
        agem = (modal(age)..., length(age))
        typem = (modal(rocktype)..., length(rocktype))
        namem = (modal(rockname)..., length(rockname))
        descm = (modal(rockdescrip)..., length(rockdescrip))

        # Print to terminal
        @info """
        $desc Macrostrat samples matched with index $i 
        
        Age [Ma]:
        \tMean: $(round(nanmean(age), digits=2))
        \tStandard deviation: $(round(nanstd(age), digits=2))
        \tMedian: $(nanmedian(age))
        \tMode: $(agem[1]) (n = $(agem[2]) of $(agem[3]))
        
        Modal sample description:
        \tRock type: $(typem[1]) (n = $(typem[2]) of $(typem[3]))
        \tRock name: $(namem[1]) (n = $(namem[2]) of $(namem[3]))
        \tRock description: $(descm[1]) (n = $(descm[2]) of $(descm[3]))
        """
    end
    export get_matched_samples

## --- End of file