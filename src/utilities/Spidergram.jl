## Fully independent (does not need to call other files) function for plotting REE spider
# diagrams

## --- Dependencies
    using Plots

## --- Definitions
    # Chondrite values from Taylor and McLennan (1985).
    chondrite_taylormclennan = (
        La = 0.367,
        Ce = 0.957,
        Pr = 0.137,
        Nd = 0.711,
        Sm = 0.231,
        Eu = 0.087,
        Gd = 0.306,
        Tb = 0.058,
        Dy = 0.381,
        Ho = 0.085,
        Er = 0.249,
        Tm = 0.036,
        Yb = 0.248,
        Lu = 0.038,
    )

    


## --- Functions
    """
    ```julia
    spidergram(data::Dict, [chondrite::NamedTuple]; 
        markershape::Symbol=:auto,
        seriescolor::Symbol=:auto,
        label::AbstractString=""
    )
    ```

    Construct a `chondrite` normalized multi-element diagrams (spider diagram) from the rare 
    earth elements in `data`. 

    Use `spidergram` to create a new plot object, and `spidergram!` to add to an existing 
    one:
    ```julia
    spidergram(args...; kw...)              # Create a new spider diagram
    spidergram!(plotobj, args...; kw...)    # Add to the plot `plotobj`
    ```

    ## Data Formatting
    REE `data` should be in a dictionary organized by element, in units of ppm (mg/g):
    ```julia
    Dict{String, Float64} with 14 entries:
    "La" => 0.001883
    "Ce" => 0.00323
    "Pr" => 0.00021
    "Nd" => 0.00142
    ⋮    => ⋮
    ```

    (Optional) Chondrite normalized values should be in a `Dictonary` organised by element.
    If no values are specified, the plot will be normalized to chondrite values from 
    Taylor and McLennan, 1985:
    ```julia
    NamedTuple with 14 elements:
    La  = Float64 0.367
    Ce  = Float64 0.957
    Pr  = Float64 0.137
    ⋮   = ⋮
    ```
    """
    function spidergram(data::Dict, chondrite::NamedTuple=chondrite_taylormclennan; 
            markershape::Symbol=:auto,
            seriescolor::Symbol=:auto,
            label::AbstractString=""
        )
        # Get REEs in chondrite normalized space
        REEᵢ = NamedTuple{Tuple(keys(chondrite))}(data[string(e)] / chondrite[e] 
            for e in keys(chondrite)
        )

        # X-Axis should skip a space for Pm
        REEs = [:La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu,]
        i = findfirst(x->x==:Pm, REEs)
        x = collect([1:i-1; i+1:length(REEs)])
        REEs = string.(REEs)
        REEs[i] = ""

        # Build plot
        h = Plots.plot(
            ylabel="Chondrite Normalized",
            fg_color_legend=:white,
            framestyle=:box,
            grid=false,
            yaxis=:log10,
            ylims=(10^0, 10^3),
            yticks=(10.0.^(0:3), ("1", "10", "100", "1000")),
            xticks=(1:length(REEs), string.(REEs)),
            yminorticks=log.(1:10),
        )
        Plots.plot!(h, x, collect(values(REEᵢ)),
            markershape=markershape, 
            seriescolor=seriescolor, msc=seriescolor,
            label=label,
        )

        return h
    end

    """
    ```julia
    spidergram!(h, data::Dict, [chondrite::NamedTuple]; 
    markershape::Symbol=:auto,
    seriescolor::Symbol=:auto,
    label::AbstractString="")
    ```

    Add a new `spidergram` to the plot object `h`.
    """
    function spidergram!(h::Plots.Plot{Plots.GRBackend}, 
            data::Dict, chondrite::NamedTuple=chondrite_taylormclennan; 
            markershape::Symbol=:auto,
            seriescolor::Symbol=:auto,
            label::AbstractString=""
        )

        # Get REEs in chondrite normalized space
        REEᵢ = NamedTuple{Tuple(keys(chondrite))}(data[string(e)] / chondrite[e] 
            for e in keys(chondrite)
        )

        # X-Axis should skip a space for Pm
        REEs = [:La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu,]
        i = findfirst(x->x==:Pm, REEs)
        x = collect([1:i-1; i+1:length(REEs)])
        REEs = string.(REEs)
        REEs[i] = ""

        Plots.plot!(h, x, collect(values(REEᵢ)),
            markershape=markershape, 
            seriescolor=seriescolor, msc=seriescolor,
            label=label,
        )

        return h
    end


## --- End of File