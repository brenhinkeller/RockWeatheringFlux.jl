# Convert a .csv of geologic time boundaries to a set of NamedTuples

## --- Main function 
    """
    ```julia
    get_GTS_boundaries()
    ```

    Return NamedTuples with the upper and lower ages (with 1σ uncertainty) for the 
    geologic time scale eons, eras, periods, epochs, and stages.

    Phanerozoic stage boundaries are from Green et al., 2022 (10.1073/pnas.2120441119). 
    Eons, eras, Precambrian periods and boundaries were added from the GTS 2012 
    (https://doi.org/10.1016/C2011-1-08249-8).

    # Example
    ```julia-repl
    eon, era, period, epoch, stage = get_GTS_boundaries()
    ```
    """
    function get_GTS_boundaries()
        # Load stage boundaries
        # Note ages are for the base of the named age / stage
        boundaries = importdataset("data/boundaries_green2022.csv", ',', importas=:Tuple)

        return (
            GTS_to_age(boundaries, lowercase.(boundaries.Eon)),
            GTS_to_age(boundaries, lowercase.(boundaries.Era)),
            GTS_to_age(boundaries, lowercase.(boundaries.Period)),
            GTS_to_age(boundaries, lowercase.(boundaries.Epoch)),
            GTS_to_age(boundaries, lowercase.(boundaries.Age_Stage_above_boundary)),
        )
    end
    export get_GTS_boundaries


## --- Supporting (non-exported) functions
    # Translate the .csv file to a NamedTuple
    function GTS_to_age(boundaries::NamedTuple, division::AbstractArray{String})
        # Get list of unique names
        unique_names = unique(division)

        # Preallocate
        result = [[NaN ± NaN, NaN ± NaN] for _ in unique_names]
        
        # Get results
        for i in eachindex(unique_names)
            result[i] = get_boundaries(unique_names[i], boundaries, division)
        end

        return NamedTuple{Tuple(Symbol.(unique_names))}(Tuple.(result))
    end


    # Get the upper and lower bounds, with uncertainty, of the GTS time scale name
    function get_boundaries(name::String, boundaries::NamedTuple, division::AbstractArray)
        age = [NaN ± NaN, NaN ± NaN]    

        # Find edge cases 
        i₀, i₁ = firstindex(division), lastindex(division)

        # Find first
        for i in eachindex(division)
            if name == division[i]
                if i == i₀
                    age[1] = boundaries.Age_Ma[i] ± boundaries.Age_sigma_Ma[i]
                    break
                else
                    age[1] = boundaries.Age_Ma[i-1] ± boundaries.Age_sigma_Ma[i-1]
                    break
                end
            end
        end

        # Find last
        for i in eachindex(division)
            if name == division[i] && i == i₁
                age[2] = boundaries.Age_Ma[i] ± boundaries.Age_sigma_Ma[i]
                break
            elseif name==division[i] && name != division[i+1]
                age[2] = boundaries.Age_Ma[i] ± boundaries.Age_sigma_Ma[i]
                break
            end
        end
        
        return age
    end


## --- End of File