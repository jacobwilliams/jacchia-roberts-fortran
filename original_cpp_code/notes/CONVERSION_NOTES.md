# C++ to Fortran Conversion Notes

## Overview

This document describes the conversion of the Jacchia-Roberts atmosphere model from C++ (GMAT implementation) to modern Fortran.

## Conversion Summary

### Language Features Mapping

| C++ Feature | Fortran Equivalent | Notes |
|-------------|-------------------|-------|
| `class JacchiaRobertsAtmosphere` | `module jacchia_roberts_module` | Module with contained procedures |
| Private member variables | `module` private variables | Module-level state in `jr_state` |
| Static const arrays | `parameter` arrays | Compile-time constants |
| Member functions | Module procedures | Contained in module |
| Real (double) | `real(dp)` with `dp = real64` | Explicit kind parameter |
| Integer | `integer(ip)` with `ip = int32` | Explicit kind parameter |
| Inheritance | Not used | Standalone module |

### Key Structural Changes

#### 1. Class to Module Conversion
**C++:**
```cpp
class JacchiaRobertsAtmosphere : public AtmosphereModel {
private:
    Real root1, root2, x_root, y_root;
    // ...
};
```

**Fortran:**
```fortran
module jacchia_roberts_module
   type :: jr_state_type
      real(dp) :: root1, root2, x_root, y_root
      ! ...
   end type jr_state_type
   
   type(jr_state_type), save :: jr_state
end module
```

#### 2. Array Indexing
- **C++:** 0-based indexing (`array[0]` to `array[n-1]`)
- **Fortran:** 1-based indexing (`array(1)` to `array(n)`)
- All array accesses were adjusted accordingly

#### 3. Array Storage Order
- **C++:** Row-major order (C-style)
- **Fortran:** Column-major order
- Multi-dimensional array `CON_DEN[5][7]` converted using `reshape` with proper ordering

#### 4. Loop Constructs
**C++:**
```cpp
for (i = 4; i >= 0; i--) {
    sum = CON_C[i] + sum * height;
}
```

**Fortran:**
```fortran
do i = 4, 1, -1
   sum = CON_C(i) + sum * height
end do
```

### Function Conversions

#### 1. JacchiaRoberts → jacchia_roberts_density
Main density calculation function:
- Simplified interface (removed count parameter for multiple spacecraft)
- Added explicit height limit checking
- Returns single density value

#### 2. exotherm
Temperature calculation:
- Converted from void function to function with return value
- Module state variables used instead of object members
- Boundary checking for trigonometric functions enhanced

#### 3. rho_100, rho_125, rho_high, rho_cor
Density calculation helpers:
- Direct 1:1 conversion
- Array indexing adjusted for Fortran
- Module state variables accessed via `jr_state%`

#### 4. roots
Polynomial root finder:
- Complex number handling using 2-element real arrays
- Array passed with explicit shape
- Convergence tolerance as parameter

#### 5. deflate_polynomial
Polynomial deflation:
- Direct conversion
- Array slicing handled explicitly

### Constants and Parameters

All physical constants converted to Fortran `parameter` with explicit `_dp` kind:

```fortran
real(dp), parameter :: RHO_ZERO = 3.46e-9_dp
real(dp), parameter :: TZERO = 183.0_dp
real(dp), parameter :: G_ZERO = 9.80665_dp
```

Array constants defined with explicit initialization:
```fortran
real(dp), parameter :: CON_C(5) = [ &
   -89284375.0_dp, &
      3542400.0_dp, &
      -52687.5_dp, &
         340.5_dp, &
          -0.8_dp ]
```

### Mathematical Functions

| C++ | Fortran | Notes |
|-----|---------|-------|
| `sqrt()` | `sqrt()` | Direct mapping |
| `exp()` | `exp()` | Direct mapping |
| `log()` | `log()` | Natural log |
| `pow(x, y)` | `x**y` | Fortran power operator |
| `fabs()` | `abs()` | Intrinsic function |
| `atan2()` | `atan2()` | Direct mapping |
| `sin()`, `cos()` | `sin()`, `cos()` | Direct mapping |

### Memory Management

- **C++:** Dynamic allocation with constructors/destructors
- **Fortran:** Module-level state variables, no dynamic allocation needed
- Initialization handled by `jr_init()` subroutine

### Removed Dependencies

The following GMAT-specific dependencies were removed or converted to placeholders:

1. **AtmosphereModel base class** → Standalone module
2. **MessageInterface** → Standard Fortran `write()` statements
3. **SolarFluxReader** → Placeholder `get_solar_flux_placeholder()`
4. **CalculateGeodetics** → Placeholder `calculate_geodetics_placeholder()`
5. **GmatConstants** → Direct constant definitions
6. **CelestialBody** → Polar radius passed to `jr_init()`

### Error Handling

**C++:**
```cpp
throw AtmosphereException("Error message");
```

**Fortran:**
```fortran
write(*,*) 'ERROR: Error message'
density = 0.0_dp
return
```

Note: Exception handling replaced with error messages and safe returns.

### Precision

Both implementations use double precision:
- **C++:** `Real` type (typically `double`)
- **Fortran:** `real(dp)` with `dp = real64`

Explicit `_dp` suffix added to all floating-point literals in Fortran.

## Notable Implementation Differences

### 1. State Management

**C++:** Object-oriented with member variables
```cpp
class JacchiaRobertsAtmosphere {
    Real t_infinity;
    Real tx;
    Real sum;
};
```

**Fortran:** Module-level state structure
```fortran
type(jr_state_type), save :: jr_state
! Access: jr_state%t_infinity
```

### 2. Initialization

**C++:** Constructor
```cpp
JacchiaRobertsAtmosphere::JacchiaRobertsAtmosphere(const std::string &name) :
   cbPolarRadius(6356.766) { }
```

**Fortran:** Explicit initialization subroutine
```fortran
subroutine jr_init(cb_polar_radius)
   real(dp), intent(in) :: cb_polar_radius
   jr_state%cb_polar_radius = cb_polar_radius
   ! ...
end subroutine
```

### 3. Array Parameters

**C++:** Pointer-based
```cpp
Real JacchiaRoberts(Real height, Real space_craft[3], Real sun[3], ...
```

**Fortran:** Assumed-shape arrays
```fortran
function jacchia_roberts_density(height, position, sun_vector, ...
   real(dp), intent(in) :: position(3)
   real(dp), intent(in) :: sun_vector(3)
```

### 4. Complex Numbers in Root Finding

**C++:** Two separate real values
```cpp
Real croots[][2];  // [real, imaginary]
```

**Fortran:** Same approach (not using native complex type)
```fortran
real(dp) :: croots(:,:)  ! (i, 1:2) where 1=real, 2=imag
```

## Testing and Validation

### Validation Points

To ensure correctness of the conversion:

1. **Constant values:** All physical constants verified against C++ code
2. **Array dimensions:** Multi-dimensional arrays carefully verified
3. **Index offsets:** All loops checked for proper 1-based indexing
4. **Mathematical expressions:** Complex expressions verified term-by-term
5. **Algorithm flow:** Control flow verified against C++ implementation

### Known Differences

1. **Error handling:** Fortran uses print statements instead of exceptions
2. **External dependencies:** Placeholder functions require implementation
3. **Multi-spacecraft:** Single spacecraft only (simplified interface)

## Placeholder Functions

Two placeholder functions need real implementations:

### 1. calculate_geodetics_placeholder
**Purpose:** Convert Cartesian position to geodetic coordinates

**Current:** Simple spherical approximation
```fortran
r = sqrt(sum(position**2))
height = r - polar_radius
latitude = atan2(z, sqrt(x^2 + y^2))
```

**Needed:** Proper ellipsoidal geodetic conversion accounting for:
- Earth's oblate spheroid shape
- Iterative latitude computation
- Proper altitude above reference ellipsoid

### 2. get_solar_flux_placeholder
**Purpose:** Get F10.7, F10.7a, and Kp values

**Current:** Returns constant values
```fortran
f107 = 150.0_dp
f107a = 150.0_dp
kp = 3.0_dp
```

**Needed:** 
- Read CSSI Space Weather files
- Interpolate values for requested epoch
- Handle historic vs. predicted data
- Support Schatten prediction model

## Build System

Two build systems provided:

### Make
```bash
make           # Build
make test      # Build and run
make clean     # Clean
```

### CMake
```bash
mkdir build && cd build
cmake ..
make
ctest         # Run tests
```

## Performance Considerations

1. **No vectorization:** Current implementation is scalar
2. **Cache efficiency:** Module state may benefit from better organization
3. **OpenMP opportunities:** Parallelization possible for multiple points
4. **Array operations:** Could use Fortran array syntax more extensively

## Future Enhancements

1. **Vectorization:** Handle multiple positions in one call
2. **Real solar flux integration:** Read actual space weather data
3. **Proper geodetics:** Implement ellipsoidal coordinate transformations
4. **Output options:** Add temperature, pressure, scale height outputs
5. **Extended species:** Output individual constituent densities
6. **Performance:** Optimize hot paths, add OpenMP
7. **Validation:** Compare against reference datasets

## Compatibility Notes

- **Fortran Standard:** F2008 (some compilers may need F2018 for full support)
- **Tested Compilers:** gfortran 9+, Intel Fortran 19+
- **Real Kind:** Uses `iso_fortran_env` for portable real kinds
- **Integer Kind:** Uses explicit kinds for portability

## References

1. Original C++ implementation: GMAT (NASA/GSFC)
2. Fortran 2008 standard: ISO/IEC 1539-1:2010
3. Numerical Recipes (for root finding algorithm reference)
4. Roberts (1971) - Original model paper

---

**Conversion Date:** 2026  
**Original C++ Source:** GMAT (General Mission Analysis Tool)  
**Fortran Version:** 2008  
**Conversion Status:** Complete, placeholder functions need implementation
