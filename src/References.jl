## --- Set up
    # Packages 
    using RockWeatheringFlux
    using HDF5, DelimitedFiles


## --- Macrostrat / Burwell 
    refs = h5read(macrostrat_io, "vars/reference")
    refs = sort(unique(refs))
    writedlm("results/references/macrostrat_ref.txt", refs)


## --- OCTOPUS
    octopusdata = importdataset("output/octopusdata.tsv",'\t', importas=:Tuple)
    refs = Array{String}(undef, length(octopusdata.pubyear))
    for i in eachindex(refs)
        refs[i] = join([octopusdata.auth[i], Int(octopusdata.pubyear[i]), 
            octopusdata.refdoi[i]], " | "
        )
    end
    refs = sort(unique(refs))
    writedlm("results/references/octopus_ref.txt", refs)
    

## --- Gard et al., 2019
    refs = h5read(geochem_fid, "bulktext/Reference")
    refs = sort(unique(refs))
    writedlm("results/references/gard_ref.txt", refs)
    

## --- End of File 