# RockWeatheringFlux.jl
Sample exposed rock types and correlate with local erosion rates as a function of slope to calculate total erosive flux

## Contents
### /data

/octopus: data from the OCTOPUS compilation (https://doi.org10.5194/essd-10-2123-2018); cosmogenic isotopes used for erosion rates. 

Current datasets are generally too large to be stored in GitHub, even with Git LFS. All data will be stored in a repository after this researh is published. Email rowan.m.gregoire.23 (at) dartmouth.edu with questions.

Lithologic data can be requested from Macrostrat (https://macrostrat.org/burwell) using the scripts in `src/ParseMacrotrat.jl`.

SRTM15+ elevation data (https://doi.org/10.1029/2019EA000658) can be downloaded from the source, or using the `get_srtm15plus` function from [StatGeochem.jl](https://github.com/brenhinkeller/StatGeochem.jl).

EarthChem geochemical data is freely available at http://portal.earthchem.org/.

### /src

All source code.

/utilities: functions specific to this project.

/visualization: figure generation and visualization of results.

### /test

Unit tests.

/volatile_sensitivity: sensitivity testing for assummptions of volatile content in sedimentary rocks.