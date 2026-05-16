# Jacchia-Roberts Atmosphere Model - Fortran Implementation

A modern Fortran implementation of the Jacchia-Roberts atmospheric density model, ported from the GMAT (General Mission Analysis Tool) C++ implementation.

## Overview

The Jacchia-Roberts atmosphere model computes atmospheric density for altitudes between 100 km and 2500 km. This model is widely used in satellite orbit determination and prediction applications.

### Reference

Roberts, C. E., Jr., "An Analytic Model for Upper Atmosphere Densities Based Upon Jacchia's 1970 Models", *Celestial Mechanics*, Vol. 4, pp. 368-377, 1971.

## Features

- **Modern Fortran 2008**: Clean, readable code using modern Fortran standards
- **Standalone Module**: No external dependencies beyond standard Fortran libraries
- **Altitude Regimes**: 
  - 90-100 km (constant density region)
  - 100-125 km (transition region)
  - 125-2500 km (high altitude region)
- **Physical Effects Modeled**:
  - Exospheric temperature variations
  - Geomagnetic activity corrections
  - Semiannual density variations
  - Seasonal-latitudinal variations
  - Solar activity effects (via F10.7 and Kp indices)

## Files

```
├── src/
│   └── jacchia_roberts_module.f90  # Main module with atmosphere model
├── test_jacchia_roberts.f90        # Test/example program
├── Makefile                         # Build configuration
└── README.md                        # This file
```

## Building the Code

### Prerequisites

- Fortran compiler (gfortran, ifort, etc.)
- Fortran Package Manager (fpm) - recommended
- Pixi package manager (optional, for managed dependencies)

### Using FPM (Recommended)

This project uses the Fortran Package Manager (fpm) for building:

```bash
# Build the library and executable with pixi
pixi run fpm build

# Run the test program
pixi run fpm run test_jacchia_roberts

# Build and run in one command
pixi run fpm run

# Build with specific compiler flags
pixi run fpm build --profile release

# Run tests (if you add them later)
pixi run fpm test
```

Or if you have fpm installed directly:
```bash
fpm build
fpm run test_jacchia_roberts
```

### Using Make (Alternative)

```bash
# Build the test program
make

# Build and run the test
make test

# Build with debug flags
make debug

# Clean build files
make clean
```

### Manual Compilation (Not Recommended)

If you prefer not to use fpm or make:
```bash
# Create directories
mkdir -p build bin

# Compile the module
pixi run gfortran -O2 -std=f2008 -c src/jacchia_roberts_module.f90 -o build/jacchia_roberts_module.o -Jbuild

# Compile the test program
pixi run gfortran -O2 -std=f2008 -c test_jacchia_roberts.f90 -o build/test_jacchia_roberts.o -Jbuild -Ibuild

# Link the executable
pixi run gfortran build/jacchia_roberts_module.o build/test_jacchia_roberts.o -o bin/test_jacchia_roberts

# Run
./bin/test_jacchia_roberts
```

Or without pixi (if gfortran is in your PATH):
```bash
# Create directories
mkdir -p build bin

# Compile the module
gfortran -O2 -std=f2008 -c src/jacchia_roberts_module.f90 -o build/jacchia_roberts_module.o -Jbuild

# Compile the test program
gfortran -O2 -std=f2008 -c test_jacchia_roberts.f90 -o build/test_jacchia_roberts.o -Jbuild -Ibuild

# Link the executable
gfortran build/jacchia_roberts_module.o build/test_jacchia_roberts.o -o bin/test_jacchia_roberts

# Run
./bin/test_jacchia_roberts
```

## Usage

### Basic Example

```fortran
program example
   use jacchia_roberts_module
   implicit none
   
   real(8) :: density, height, position(3), sun_vector(3)
   real(8) :: geo_lat, sun_dec, utc_mjd
   type(geoparms_type) :: geo
   
   ! Initialize module with Earth's polar radius (km)
   call jr_init(6356.766d0)
   
   ! Set up parameters
   height = 400.0d0                    ! Altitude in km
   position = [7000.0d0, 0.0d0, 0.0d0] ! Position vector (km)
   sun_vector = [1.0d0, 0.0d0, 0.0d0]  ! Unit vector to Sun
   geo_lat = 0.0d0                     ! Geodetic latitude (degrees)
   sun_dec = 0.0d0                     ! Sun declination (radians)
   utc_mjd = 51544.5d0                 ! Modified Julian Date
   
   ! Geomagnetic parameters
   geo%tkp = 3.0d0                     ! Kp index
   geo%xtemp = 865.0d0                 ! Exospheric temperature (K)
   
   ! Calculate density (returns kg/m^3)
   density = jacchia_roberts_density(height, position, sun_vector, &
                                     geo_lat, sun_dec, geo, utc_mjd)
   
   write(*,*) 'Density at', height, 'km:', density, 'kg/m^3'
   
end program example
```

### Calculating Exospheric Temperature

The exospheric temperature is calculated from solar flux indices:

```fortran
! Given F10.7 and F10.7a (81-day average)
real(8) :: f107, f107a, xtemp

f107 = 150.0d0   ! Daily F10.7 flux
f107a = 150.0d0  ! 81-day average

xtemp = 379.0d0 + 3.24d0 * f107a + 1.3d0 * (f107 - f107a)
geo%xtemp = xtemp
```

## Module Interface

### Public Subroutines and Functions

#### `jr_init(cb_polar_radius)`
Initialize the module with the central body's polar radius.

**Arguments:**
- `cb_polar_radius` (real(8), in): Polar radius in km (6356.766 for Earth)

#### `jacchia_roberts_density(...)`
Calculate atmospheric density.

**Arguments:**
- `height` (real(8), in): Altitude above reference ellipsoid (km)
- `position(3)` (real(8), in): Position vector in TOD/GCI frame (km)
- `sun_vector(3)` (real(8), in): Unit vector to Sun in TOD/GCI frame
- `geo_lat` (real(8), in): Geodetic latitude (degrees)
- `sun_dec` (real(8), in): Sun declination (radians)
- `geo` (geoparms_type, in): Geomagnetic parameters structure
- `utc_mjd` (real(8), in): UTC Modified Julian Date

**Returns:**
- `density` (real(8)): Atmospheric density in kg/m³

### Types

#### `geoparms_type`
Structure containing geomagnetic parameters:
- `tkp` (real(8)): Geomagnetic Kp index
- `xtemp` (real(8)): Exospheric temperature (K)

## Placeholder Functions

The following functions are currently implemented as placeholders and should be replaced with actual implementations for production use:

### `calculate_geodetics_placeholder`
Converts Cartesian position to geodetic coordinates. Currently uses a simplified spherical approximation.

**TODO:** Implement proper geodetic coordinate transformation accounting for Earth's oblate spheroid.

### `get_solar_flux_placeholder`
Retrieves solar flux (F10.7) and geomagnetic (Kp) indices. Currently returns constant nominal values.

**TODO:** Implement reading from solar flux data files (e.g., CSSI Space Weather files) or connect to a space weather data service.

## Model Validity and Limitations

### Valid Range
- **Altitude:** 100 km to 2500 km
- Below 100 km: Model returns zero density (not valid)
- Above 2500 km: Model returns zero density (exosphere)

### Required Inputs
- Accurate solar flux indices (F10.7, F10.7a)
- Geomagnetic activity index (Kp)
- Precise spacecraft position
- Sun position

### Coordinate Systems
- Position vectors should be in True-of-Date (TOD) or Geocentric Celestial Inertial (GCI) frame
- Sun vectors should be in the same frame as position vectors

## Physical Constants

The model uses the following physical constants (as per the original GMAT implementation):

- **RHO_ZERO:** 3.46×10⁻⁹ g/cm³ (low altitude density)
- **TZERO:** 183 K (temperature at 90 km)
- **G_ZERO:** 9.80665 m/s² (gravitational acceleration)
- **GAS_CON:** 8.31432 J/(K·mol) (gas constant)
- **AVOGADRO:** 6.022045×10²³ (Avogadro's number)

Note: These values are tuned for consistency with the GMAT implementation and may differ slightly from standard references.

## Atmospheric Constituents

The model accounts for six atmospheric species:
1. Nitrogen (N₂)
2. Argon (Ar)
3. Helium (He)
4. Molecular Oxygen (O₂)
5. Atomic Oxygen (O)
6. Hydrogen (H) - only above 500 km

## Testing

The included test program (`test_jacchia_roberts.f90`) demonstrates:
- Module initialization
- Density calculation at multiple altitudes
- Generation of a density profile from 100-1000 km

Run the test with fpm:
```bash
pixi run fpm run test_jacchia_roberts
```

Or with make:
```bash
make test
```

This will generate a file `density_profile.dat` with density values at 10 km intervals.

## FPM Package Structure

This project uses the Fortran Package Manager (fpm) for streamlined building and distribution:

- **fpm.toml** - Package manifest with metadata and build configuration
- **src/** - Library source files (the main module)
- **test_jacchia_roberts.f90** - Example executable in the root directory

To use this as a dependency in another fpm project:
```toml
[dependencies]
jacchia-roberts = { git = "https://github.com/yourusername/jacchia-roberts-fortran.git" }
```

## Original Source

This Fortran implementation is ported from:
- **Project:** GMAT (General Mission Analysis Tool)
- **Copyright:** 2002-2026 United States Government (NASA)
- **License:** Apache License 2.0
- **Original Authors:** Waka A. Waktola, Darrel J. Conway, and GMAT development team

## License

This Fortran implementation maintains the Apache License 2.0 of the original GMAT code.

## Further Development

### Recommended Enhancements

1. **Solar Flux Data Integration:**
   - Implement reading from CSSI Space Weather files
   - Add support for Schatten prediction model
   - Implement automatic flux data interpolation

2. **Geodetic Transformations:**
   - Implement proper geodetic coordinate calculations
   - Add support for different reference ellipsoids
   - Include coordinate frame transformations

3. **Performance Optimizations:**
   - Vectorize density calculations for multiple points
   - Cache frequently used calculations
   - Implement OpenMP parallelization for large-scale calculations

4. **Extended Functionality:**
   - Add temperature and pressure outputs
   - Implement scale height calculations
   - Add partial density outputs for individual species

5. **Validation and Testing:**
   - Compare against reference data
   - Add unit tests for all functions
   - Validate against GMAT implementation

## Contact and Support

For questions, bug reports, or contributions, please refer to the project repository.

## References

1. Roberts, C. E., Jr., "An Analytic Model for Upper Atmosphere Densities Based Upon Jacchia's 1970 Models", *Celestial Mechanics*, Vol. 4, pp. 368-377, 1971.

2. Jacchia, L. G., "New Static Models of the Thermosphere and Exosphere with Empirical Temperature Profiles", *Smithsonian Astrophysical Observatory Special Report No. 313*, 1970.

3. Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th Edition, Microcosm Press, 2013.

4. GMAT Documentation: [https://software.nasa.gov/software/GSC-17177-1](https://software.nasa.gov/software/GSC-17177-1)
