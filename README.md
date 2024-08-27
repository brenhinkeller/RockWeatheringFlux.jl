# RockWeatheringFlux.jl
Estimate the composition of the continental crust, and correlate this estimate with erosion rate to estimate volume and composition of eroded material.

In addition to performing computations and defining its own functions, RockWeatheringFlux reexports:

* [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl)
* [StatGeochem.jl](https://github.com/brenhinkeller/StatGeochem.jl)
* [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl): `@showprogress`, `Progress`, `next!`


## Running the Code
Some setup and definitions file (e.g. `src/utilities/Definitions.jl`) call bash commands in order to create directories, which causes errors on Windows machines. The easiest way to run this code on these machines is to use [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install).

## Contents
### `/data`

`/octopus`: data from the OCTOPUS compilation (https://doi.org10.5194/essd-10-2123-2018); cosmogenic isotopes used for erosion rates. 

### `/src`

All source code. See `run_endtoend.sh` for the recommended code run order (does not include figures).

`/utilities`: functions specific to this project.

* See `src/utilities/ScreenOutliers.jl` for criteria used to remove outliers from geochemical data.

`/visualization`: figure generation and visualization of results.

`/volatile_sensitivity`: sensitivity testing for assummptions of volatile content in sedimentary rocks.

### `/test`

Unit tests.
