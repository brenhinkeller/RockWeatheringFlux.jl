# `/visualization/chapter1_composition`
Scripts to visualize results from the lithology / geochemistry method and associated estimation of the composition of the uppermost continental crust.

## Organization
Each figure (or "family" of related figures) has its own file. Common packages and data are loaded in `Definitions.jl`. Unique packages and data are loaded by file.

## Contents
 * `ArcheanRocks.jl`: Experiments with distributions of rock class and silica in Archean (> 2500 Ma.) samples.
 * `GeochemLoc.jl`: Age and location of samples in the geochemical dataset.
 * `Phosphorous.jl`: Changes in the concentration of phosphorus over geologic time.
 * `SilicaAgeDistribution.jl`: 2D histograms of EarthChem and matched samples as a function of age and silica content, after Fig. 6.9 from Keller (2016).
 * `SilicaDistribution.jl`: Silica distribution by rock class.
 * `SlopeErosion.jl`: Erosion as a function of slope.
 * `Spidergram.jl`: REE patterns in my estimation of the bulk continental crust compared to those from Rudnick and Gao, 2014.
 * `VolatileSensitivity.jl`: Results from sensitivity testing the impact of wt.% assumed volatiles in data screening.
