# Bugs Fixed in Fortran Port of Jacchia-Roberts Atmosphere Model

This document describes critical bugs discovered in the original GMAT C++ implementation of the Jacchia-Roberts atmosphere model that have been corrected in this Fortran port.

## Bug 1: Division by Zero in Helium Correction (Equinox Case)

### Location
- **C++ Code**: `JacchiaRobertsAtmosphere.cpp`, line 1277-1279
- **Fortran Code**: `jacchia_roberts_module.f90`, `rho_high` function

### Description
When computing atmospheric density for altitudes above 125 km, the model applies a helium correction factor that depends on solar declination:

```cpp
// Original C++ code (BUGGY)
if (i == 2)  // Helium
{
   exp1 -= 0.38;
   f = 4.9914 * fabs(sun_dec) * (pow(sin(0.25 * GmatMathConstants::PI - 0.5 * geo_lat *
                                         sun_dec / fabs(sun_dec)), 3)
                                 - 0.35355) / GmatMathConstants::PI;
   f = pow(10.0, f);
}
```

The expression `sun_dec / fabs(sun_dec)` is used to extract the sign of the sun declination. However, when `sun_dec = 0` (which occurs at the equinoxes when the sun is over the equator), this becomes **0 / 0 = NaN**.

### Impact
- **Severity**: Critical
- **Conditions**: Occurs whenever sun declination is zero (equinoxes)
- **Result**: All density calculations return NaN, making the model completely unusable for equinox conditions
- **Propagation**: The NaN value propagates through all subsequent calculations, rendering the entire density computation invalid

### Fix
The Fortran implementation detects this edge case and handles it properly:

```fortran
! Fixed Fortran code
if (i == 3) then  ! Helium
   exp1 = exp1 - 0.38_dp
   ! Handle sun_dec = 0 case to avoid NaN from 0/0
   if (abs(sun_dec) < 1.0e-15_dp) then
      ! When sun_dec = 0, the formula simplifies to f = 10^0 = 1
      f = 1.0_dp
   else
      f = 4.9914_dp * abs(sun_dec) * &
          (sin(0.25_dp * PI - 0.5_dp * geo_lat * sun_dec / abs(sun_dec))**3 - &
           0.35355_dp) / PI
      f = 10.0_dp**f
   end if
end if
```

When `sun_dec ≈ 0`, the mathematically correct value is `f = 1.0` (since `4.9914 * 0 * ... = 0` and `10^0 = 1`).

---

## Bug 2: Overly Restrictive Collinear Geometry Check

### Location
- **C++ Code**: `JacchiaRobertsAtmosphere.cpp`, lines 814-820
- **Fortran Code**: `jacchia_roberts_module.f90`, `exotherm` function

### Description
The original C++ code throws an exception when the spacecraft position and sun vectors are collinear in the xy-plane:

```cpp
// Original C++ code (OVERLY RESTRICTIVE)
sun_denom = sqrt(sun[0]*sun[0] + sun[1]*sun[1]);
cross_denom = fabs(sun[0] * space_craft[1] - sun[1] * space_craft[0]);
cos_denom = sqrt(space_craft[0]*space_craft[0] +
      space_craft[1]*space_craft[1]);

if ((cross_denom < GmatRealConstants::REAL_TOL) || (cos_denom < GmatRealConstants::REAL_TOL))
   throw AtmosphereException("Numerical precision error in "
         "JacchiaRobertsAtmosphere::exotherm; denominator is too close to 0.0");
```

This check prevents the model from running when the cross product is near zero, which occurs when the sun and spacecraft are collinear (aligned or opposite) in the equatorial plane.

### Impact
- **Severity**: High
- **Conditions**: Occurs when spacecraft is directly sunward or anti-sunward from Earth in the equatorial plane
- **Result**: Model throws an exception and terminates
- **Physical Reality**: This is a **physically valid geometry** that occurs regularly in orbital mechanics

### Why This is Overly Restrictive
The cross product is used to determine the **sign** of the hour angle. However, when vectors are collinear:
- If they point in the same direction (cosAlpha ≈ 1): hour angle = 0
- If they point in opposite directions (cosAlpha ≈ -1): hour angle = ±π

In both cases, the hour angle can be determined from the dot product alone, without needing the cross product.

### Fix
The Fortran implementation handles collinear geometry gracefully:

```fortran
! Fixed Fortran code
! Calculate cosine of angle between vectors
cosAlpha = (sun_vector(1) * position(1) + sun_vector(2) * position(2)) / &
           (sun_denom * cos_denom)

! Compute cross product for sign determination
cross_denom = abs(sun_vector(1) * position(2) - sun_vector(2) * position(1))

! Calculate hour angle with boundary checking
if (cosAlpha >= 1.0_dp - error_tolerance) then
   ! Vectors aligned: hour angle = 0
   hour_angle = 0.0_dp
else if (cosAlpha <= -1.0_dp + error_tolerance) then
   ! Vectors opposite: hour angle = ±π
   if (cross_denom < real_tol) then
      hour_angle = PI
   else
      hour_angle = sign(1.0_dp, sun_vector(1) * position(2) - &
                       sun_vector(2) * position(1)) * PI
   end if
else
   ! General case: use cross product to determine sign
   if (cross_denom < real_tol) then
      hour_angle = 0.0_dp
   else
      hour_angle = sign(1.0_dp, sun_vector(1) * position(2) - &
                       sun_vector(2) * position(1)) * acos(cosAlpha)
   end if
end if
```

This allows the model to handle all physically valid geometries, including collinear cases.

---

## Additional Fix: Array Layout Correction

While not a bug in the C++ code itself, the Fortran port required careful attention to array layout differences and indexing:

### Issue 1: CON_DEN Array Layout
The `CON_DEN` array stores polynomial coefficients for atmospheric constituent densities. C++ uses row-major order, while Fortran uses column-major order by default.

### Fix
Added explicit ordering to the `reshape()` function:

```fortran
real(dp), parameter :: CON_DEN(5,7) = reshape([ &
   ! ... data ...
   ], [5, 7], order=[2, 1])  ! Specify row-major order to match C++
```

Without this correction, the polynomial evaluation produced incorrect values (underflow to zero or overflow to infinity).

### Issue 2: Array Index Offset
C++ uses 0-based array indexing while Fortran uses 1-based indexing. During the port, one instance was missed:

**Bug**: In the `rho_125` function, the code used `CON_C(4)` instead of `CON_C(5)`:
```fortran
! INCORRECT (uses 4th element = 340.5):
factor_k = -1500625.0_dp * G_ZERO * jr_state%cb_polar_squared / &
           (GAS_CON * CON_C(4) * (jr_state%tx - TZERO))

! CORRECT (uses 5th element = -0.8):
factor_k = -1500625.0_dp * G_ZERO * jr_state%cb_polar_squared / &
           (GAS_CON * CON_C(5) * (jr_state%tx - TZERO))
```

The C++ code uses `CON_C[4]` which refers to the 5th element (value: -0.8), so the Fortran code must use `CON_C(5)` to access the same element.

All other array indices were correctly converted from 0-based to 1-based indexing. See [ARRAY_INDEX_VERIFICATION.md](ARRAY_INDEX_VERIFICATION.md) for a complete verification of all array accesses.

---

## Testing

Both bugs were discovered and fixed during testing with:
- **Sun declination**: 0 radians (equinox)
- **Spacecraft position**: [7000, 0, 0] km (collinear with x-axis)
- **Sun vector**: [1, 0, 0] (collinear with spacecraft position)

The original C++ code would fail under these conditions, while the Fortran implementation produces physically correct atmospheric density values.

---

## Recommendations

These bugs should be fixed in the original GMAT C++ implementation to improve robustness and handle all physically valid input conditions.
