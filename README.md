A modern Fortran implementation of the Jacchia-Roberts atmospheric density model.

## Status

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![CI Status](https://github.com/jacobwilliams/jacchia-roberts-fortran/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/jacchia-roberts-fortran/actions)
[![GitHub release](https://img.shields.io/github/release/jacobwilliams/jacchia-roberts-fortran.svg?style=plastic)](https://github.com/jacobwilliams/jacchia-roberts-fortran/releases/latest)

## Overview

The Jacchia-Roberts atmosphere model computes atmospheric density for altitudes between 100 km and 2500 km. This model is widely used in satellite orbit determination and prediction applications. The model automatically retrieves F10.7 solar flux and Kp geomagnetic indices from a user-specified CSSI Space Weather file.

## Building the Code

This project uses the Fortran Package Manager (fpm) for building:

```bash
# activate the environment
pixi shell

# Build the library and test
fpm build

# Run the test program
fpm test

# Build with specific compiler flags
fpm build --profile release
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following processor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To use this as a dependency in another fpm project:
```toml
[dependencies]
jacchia-roberts = { git = "https://github.com/jacobwilliams/jacchia-roberts-fortran.git" }
```

## Example

The following example program shows how to use the model:

```fortran
program example

use jacchia_roberts_module, only: jacchia_roberts_type
use jacchia_roberts_kinds,  only: ip, dp

implicit none

type(jacchia_roberts_type) :: jr
integer(ip) :: status
real(dp) :: density

! example inputs:
real(dp),parameter :: rad_earth = 6356.766_dp ! Earth polar radius (km)
character(len=*),parameter :: sw_file = 'data/SpaceWeather-All-v1.2.txt' ! space weather file to load
real(dp),parameter :: utc_mjd = 59215.5_dp ! MJD: Jan 1, 2021 12:00 UTC
real(dp),parameter :: alt_km = 200.0_dp ! geodetic altitude (km)
real(dp),parameter :: lat = 20.0_dp ! geodetic latitude (deg)
real(dp),dimension(3),parameter :: r = [7000.0_dp, 0.0_dp, 0.0_dp] ! spacecraft position vector (km)
real(dp),dimension(3),parameter :: sun = [1.0_dp, 0.0_dp, 0.0_dp] ! sun direction unit vector

! initialize the model:
call jr%initialize(rad_earth, sw_file, status)

! compute the density:
density = jr%density(alt_km, r, sun_vector, lat, utc_mjd)

end program example
```

## Documentation

The latest API documentation can be found [here](https://jacobwilliams.github.io/jacchia-roberts-fortran/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## References

1. Roberts, C. E., Jr., "[An Analytic Model for Upper Atmosphere Densities Based Upon Jacchia's 1970 Models](https://ui.adsabs.harvard.edu/scan/manifest/1971CeMec...4..368R)", *Celestial Mechanics*, Vol. 4, pp. 368-377, 1971.
2. Jacchia, L. G., "[New Static Models of the Thermosphere and Exosphere with Empirical Temperature Profiles](https://articles.adsabs.harvard.edu/pdf/1970SAOSR.313.....J)", *Smithsonian Astrophysical Observatory Special Report No. 313*, 1970.
3. Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th Edition, Microcosm Press, 2013.
4. [GMAT](https://github.com/nasa/GMAT) -- The fortran code is based on the C++ code from here
5. [CelesTrak](https://celestrak.org/SpaceData/) -- Source of space weather files.