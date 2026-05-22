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

To use this as a dependency in another fpm project:
```toml
[dependencies]
jacchia-roberts = { git = "https://github.com/jacobwilliams/jacchia-roberts-fortran.git" }
```

## Documentation

The latest API documentation can be found [here](https://jacobwilliams.github.io/jacchia-roberts-fortran/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## References

1. Roberts, C. E., Jr., "An Analytic Model for Upper Atmosphere Densities Based Upon Jacchia's 1970 Models", *Celestial Mechanics*, Vol. 4, pp. 368-377, 1971.

2. Jacchia, L. G., "New Static Models of the Thermosphere and Exosphere with Empirical Temperature Profiles", *Smithsonian Astrophysical Observatory Special Report No. 313*, 1970.

3. Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th Edition, Microcosm Press, 2013.

4. GMAT: https://github.com/nasa/GMAT

5. CelesTrak: https://celestrak.org/SpaceData/