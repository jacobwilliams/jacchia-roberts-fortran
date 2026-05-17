# Jacchia-Roberts Atmosphere Model - Fortran Implementation

A modern Fortran implementation of the Jacchia-Roberts atmospheric density model, ported from the GMAT (General Mission Analysis Tool) C++ implementation.

## Overview

The Jacchia-Roberts atmosphere model computes atmospheric density for altitudes between 100 km and 2500 km. This model is widely used in satellite orbit determination and prediction applications.

### Reference

Roberts, C. E., Jr., "An Analytic Model for Upper Atmosphere Densities Based Upon Jacchia's 1970 Models", *Celestial Mechanics*, Vol. 4, pp. 368-377, 1971.

## Features

- **Modern Fortran 2008**: Clean, readable code using modern Fortran standards
- **Space Weather Integration**: Reads CSSI Space Weather files from CelesTrak
- **High Performance**: O(1) direct indexing for daily data, efficient monthly data lookup
- **Altitude Regimes**:
  - 90-100 km (constant density region)
  - 100-125 km (transition region)
  - 125-2500 km (high altitude region)
- **Physical Effects Modeled**:
  - Exospheric temperature variations (from F10.7 solar flux)
  - Geomagnetic activity corrections (from Kp indices)
  - Semiannual density variations
  - Seasonal-latitudinal variations
  - Solar activity effects automatically retrieved from space weather data

## Files

```
├── src/
│   ├── jacchia_roberts_module.f90   # Main module with atmosphere model
│   └── space_weather_module.f90     # Space weather data reader (CSSI format)
├── test/
│   └── test_jacchia_roberts.f90     # Test/example program
├── data/
│   └── SpaceWeather-All-v1.2.txt    # CSSI Space Weather data file
├── original_cpp_code/               # Reference C++ implementation from GMAT
│   ├── JacchiaRobertsAtmosphere.cpp
│   ├── JacchiaRobertsAtmosphere.hpp
│   ├── SolarFluxReader.cpp
│   └── SolarFluxReader.hpp
├── fpm.toml                         # FPM package manifest
├── pixi.toml                        # Pixi environment configuration
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
# activate the environment
pixi shell

# Build the library and test
fpm build

# Run the test program
fpm test

# Build with specific compiler flags
fpm build --profile release
```

## Usage

```fortran
integer :: sw_status

   ! Initialize module with Earth's polar radius (km)
   call jr_init(6356.766d0)

   ! Load space weather data from CSSI file
   call jr_load_space_weather('data/SpaceWeather-All-v1.2.txt', sw_status)
   if (sw_status /= 0) then
      write(*,*) 'Warning: Could not load space weather data'
      write(*,*) '         Using nominal values'
   end if

   ! Set up parameters
   height = 400.0d0                    ! Altitude in km
   position = [7000.0d0, 0.0d0, 0.0d0] ! Position vector (km)
   sun_vector = [1.0d0, 0.0d0, 0.0d0]  ! Unit vector to Sun
   geo_lat = 0.0d0                     ! Geodetic latitude (degrees)
   sun_dec = 0.0d0                     ! Sun declination (radians)
   utc_mjd = 59215.5d0                 ! Modified Julian Date (Jan 1, 2021)

   ! Calculate density (returns kg/m^3)
   ! Space weather data (F10.7, Kp) is automatically retrieved for the given date
   density = jacchia_roberts_density(height, position, sun_vector, &
                                     geo_lat, sun_dec, utc_mjd)

   write(*,*) 'Density at', height, 'km:', density, 'kg/m^3'

   ! Clean up when done
   call jr_cleanup()e

   ! Geomagnetic parameters
   geo%tkp = 3.0d0                     ! Kp index
   geo%xtemp = 865.0d0                 ! Exospheric temperature (K)

   ! Calculate density (returns kg/m^3)
   density = jacchia_roberts_density(height, position, sun_vector, &
    Space Weather Data
```

The model automatically retrieves F10.7 solar flux and Kp geomagnetic indices from the loaded space weather file. The exospheric temperature is calculated internally using the formula:

```
T_exo = 379.0 + 3.24 * F10.7a + 1.3 * (F10.7 - F10.7a)
```

where:
- F10.7 is the observed daily solar flux (10⁻²² W/m²/Hz)
- F10.7a is the 81-day centered average

The space weather data file contains:
- **Observed data**: Historical measurements from 1957-10-01 to present
- **Daily predicted**: 45-day forecasts from NOAA SWPC
- **Monthly predicted**: Long-term predictions extending to ~2109

Data coverage: ~25,545 records spanning October 1957 through 2109.l(8) :: f107, f107a, xtemp

f107 = 150.0d0   ! Daily F10.7 flux
f107a = 150.0d0Jacchia-Roberts module with the central body's polar radius.

**Arguments:**
- `cb_polar_radius` (real(8), in): Polar radius in km (6356.766 for Earth)

#### `jr_load_space_weather(filename, status)`
Load space weather data from a CSSI Space Weather file.

**Arguments:**
- `filename` (character, in): Path to CSSI Space Weather file
- `status` (integer, out): Status code (0=success, non-zero=error)

**Notes:**
- Must be called after `jr_init()` and before `jacchia_roberts_density()`
- Loads observed, daily predicted, and monthly predicted data
- Uses the same algorithm as GMAT's SolarFluxReader for compatibility

#### `jacchia_roberts_density(...)`
Calculate atmospheric density at a given position and time.

**Arguments:**
- `height` (real(8), in): Altitude above reference ellipsoid (km)
- `position(3)` (real(8), in): Position vector in TOD/GCI frame (km)
- `sun_vector(3)` (real(8), in): Unit vector to Sun in TOD/GCI frame
- `geo_lat` (real(8), in): Geodetic latitude (degrees)
- `sun_dec` (real(8), in): Sun declination (radians)
- `utc_mjd` (real(8), in): UTC Modified Julian Date
CSSI Space Weather File Format

The module reads CSSI Space Weather files from CelesTrak in the legacy fixed-width format:

**Format:** `FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,F4.1,I2,I4,F6.1,I2,5F6.1)`

**File Structure:**
- Header with metadata and column descriptions
- `BEGIN OBSERVED` / `END OBSERVED` - Historical measurements
- `BEGIN DAILY_PREDICTED` / `END DAILY_PREDICTED` - Near-term forecasts (45 days)
- `BEGIN MONTHLY_PREDICTED` / `END MONTHLY_PREDICTED` - Long-term predictions

**Data Fields:**
- Date (year, month, day)
- Bartels Solar Rotation Number and day number
- 8 three-hourly Kp indices (multiplied by 10)
- 8 three-hourly Ap indices
- Daily average Ap
- F10.7 solar flux (observed and adjusted)
- 81-day centered averages

**Download:** Space weather files are available from CelesTrak at:
https://celestrak.org/SpaceData/

The included `data/SpaceWeather-All-v1.2.txt` file provides coverage from 1957 to 2109
### Space Weather Module

The `space_weather_module` provides the underlying space weather data functionality:

#### `sw_init(filename, status)`
Low-level initialization of space weather data (called by `jr_load_space_weather`)

#### `sw_get_flux_data(mjd, flux_data, status)`
Retrieve space weather data for a specific Modified Julian Date

#### `flux_data_type`
Structure containing space weather parameters:
- `mjd` (real(8)): Modified Julian Date
- `f107_obs` (real(8)): Observed F10.7 solar flux
- `f107_adj` (real(8)): Adjusted F10.7 solar flux
- `f107a_obs_ctr` (real(8)): Observed F10.7 81-day centered average
- `f107a_adj_ctr` (real(8)): Adjusted F10.7 81-day centered average
- `kp(8)` (real(8)): Kp indices for 8 3-hour periods
- `ap_avg` (real(8)): Daily average Ap indexters structure
- `utc_mjd` (real(8), in): UTC Modified Julian Date

**Returns:**
- `density` (real(8)): Atmospheric density in kg/m³

### Types

#### `geoparms_type`
Structure containing geomagnetic parameters:
- `tkp` (real(8)): Geomagnetic Kp index
- `xtemp` (real(8)): Exospheric temperature (K)
/test_jacchia_roberts.f90`) demonstrates:
- Module initialization with Earth's polar radius
- Loading space weather data from CSSI file
- Density calculation at multiple altitudes (150, 250, 400, 600, 1000 km)
- Generation of a density profile from 100-1000 km

Run the test with fpm:
```bash
fpm test
```

The test will:
1. Load 25,545 space weather records from `data/SpaceWeather-All-v1.2.txt`
2. Calculate densities for Jan 1, 2021 (MJD 59215.5)
3. Generate `density_profile.dat` with density values at 10 km intervals

**Expected Output:**
```
Loading space weather data...
Space weather data loaded: 25545 records
  Date range (MJD):     36112.50 to     86882.50

Test Scenario:
  Position magnitude:        7000.00 km
  UTC MJD:                   59215.jacchia_roberts_module, space_weather_module)
- **test/** - Test programs and examples
- **data/** - Space weather data files
```

To use this as a dependency in another fpm project:
```toml
[dependencies]
jacchia-roberts = { git = "https://github.com/jacobwilliams/jacchia-roberts-fortran.git" }
```

**Note:** When using as a dependency, you'll need to provide your own CSSI Space Weather data file.

## Testing

The included test program (`test_jacchia_roberts.f90`) demonstrates:
- Module initialization
- Density calculation at multiple altitudes
- Generation of a density profile from 100-1000 km

Run the test with fpm:
```bash
fpm test
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

## References

1. Roberts, C. E., Jr., "An Analytic Model for Upper Atmosphere Densities Based Upon Jacchia's 1970 Models", *Celestial Mechanics*, Vol. 4, pp. 368-377, 1971.

2. Jacchia, L. G., "New Static Models of the Thermosphere and Exosphere with Empirical Temperature Profiles", *Smithsonian Astrophysical Observatory Special Report No. 313*, 1970.

3. Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th Edition, Microcosm Press, 2013.

4. GMAT Documentation: [https://software.nasa.gov/software/GSC-17177-1](https://software.nasa.gov/software/GSC-17177-1)
