# RockWeatheringFlux.jl
Estimate the composition of the continental crust, and correlate this estimate with erosion rate to estimate volume and composition of eroded material.

In addition to performing computations and defining its own functions, RockWeatheringFlux reexports (and internally depends on):

* [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl)
* [StatGeochem.jl](https://github.com/brenhinkeller/StatGeochem.jl)
* [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl): `@showprogress`, `Progress`, `next!`

## Contents
### `/data`

`/octopus`: data from the OCTOPUS compilation (https://doi.org10.5194/essd-10-2123-2018); cosmogenic isotopes used for erosion rates. 

### `/src`

All source code.

`/utilities`: functions specific to this project.

`/visualization`: figure generation and visualization of results.

`/volatile_sensitivity`: sensitivity testing for assummptions of volatile content in sedimentary rocks.

### `/test`

Unit tests.

## Data Availability
Datasets are too large to be stored in GitHub, even with Git LFS. All data will be stored in a non-GitHub repository when this research is published. Email rowan.m.gregoire.23 (at) dartmouth.edu with questions.

### Current Data Access

Erosion rates from the OCTOPUS compilation (https://doi.org/10.5194/essd-10-2123-2018) can be downloaded from the source, or from the `.gz` file in `/data/octopus`.

Lithologic data can be requested from Macrostrat (https://macrostrat.org/burwell) using the scripts in `src/ParseMacrotrat.jl`.

SRTM15+ elevation data (https://doi.org/10.1029/2019EA000658) can be downloaded from the source, or using the `get_srtm15plus` function from [StatGeochem.jl](https://github.com/brenhinkeller/StatGeochem.jl).

EarthChem geochemical data is freely available at http://portal.earthchem.org/.