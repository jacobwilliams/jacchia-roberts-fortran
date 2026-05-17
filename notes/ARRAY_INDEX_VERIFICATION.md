# Array Index Verification: C++ (0-based) to Fortran (1-based)

## Array Definitions

| Array | C++ Size | C++ Indices | Fortran Size | Fortran Indices | Status |
|-------|----------|-------------|--------------|-----------------|--------|
| CON_C | [5] | 0-4 | (5) | 1-5 | ✓ |
| CON_L | [5] | 0-4 | (5) | 1-5 | ✓ |
| M_CON | [7] | 0-6 | (7) | 1-7 | ✓ |
| S_CON | [6] | 0-5 | (6) | 1-6 | ✓ |
| S_BETA | [6] | 0-5 | (6) | 1-6 | ✓ |
| ZETA_CON | [7] | 0-6 | (7) | 1-7 | ✓ |
| MOL_MASS | [6] | 0-5 | (6) | 1-6 | ✓ |
| NUM_DENS | [5] | 0-4 | (5) | 1-5 | ✓ |
| CON_DEN | [5][7] | [0-4][0-6] | (5,7) | (1-5,1-7) | ✓ |

## Critical Index Mappings

### CON_C Array
- **C++**: `CON_C[4]` = -0.8 (5th element)
- **Fortran**: `CON_C(5)` = -0.8 (5th element) ✓ **FIXED**
- **Location**: `rho_125` function, factor_k calculation
- **Original Bug**: Used `CON_C(4)` = 340.5 instead of `CON_C(5)` = -0.8

### CON_L Array
- **C++**: `CON_L[4]` = 0.2462708e-9 (5th element)
- **Fortran**: `CON_L(5)` = 0.2462708e-9 (5th element) ✓

### M_CON Array
- **C++**: `M_CON[6]` = -0.000000697444 (7th element)
- **Fortran**: `M_CON(7)` = -0.000000697444 (7th element) ✓

### ZETA_CON Array
- **C++**: `ZETA_CON[6]` = 0.7026942e-32 (7th element)
- **Fortran**: `ZETA_CON(7)` = 0.7026942e-32 (7th element) ✓

### CON_DEN Array (2D)
- **C++**: `CON_DEN[i][6]` (7th column)
- **Fortran**: `CON_DEN(i,7)` (7th column) ✓
- **Note**: Also required `order=[2,1]` in reshape() to match row-major layout

## Loop Index Conversions

### rho_125 Function
**C++**:
```cpp
for (rho_sum = 0.0, i=0; i<=4; i++) {
    rhoi = MOL_MASS[i] * NUM_DENS[i] * ...
    if (i == 2) {  // Helium (3rd element)
        rhoi *= pow(t_100/temperature, -0.38);
    }
}
```

**Fortran**:
```fortran
do i = 1, 5
    rhoi = MOL_MASS(i) * NUM_DENS(i) * ...
    if (i == 3) then  ! Helium (3rd element)
        rhoi = rhoi * (t_100 / temperature)**(-0.38_dp)
    end if
end do
```
✓ Indices correctly shifted: C++ `i=0..4` → Fortran `i=1..5`, C++ `i==2` → Fortran `i==3`

### rho_high Function
**C++**:
```cpp
for (rho_out = 0.0, i = 0; i <= 5; i++) {
    if (i <= 4) { ... }  // First 5 constituents
    if (i == 2) { ... }  // Helium (3rd element)
    if (height > 500.0 && i == 5) {  // Hydrogen (6th element)
        r = MOL_MASS[5] * ...
    }
    else if (i <= 4) {
        r = ... * MOL_MASS[i] * ...
    }
}
```

**Fortran**:
```fortran
do i = 1, 6
    if (i <= 5) then ... end if  ! First 5 constituents
    if (i == 3) then ... end if  ! Helium (3rd element)
    if (height > 500.0_dp .and. i == 6) then  ! Hydrogen (6th element)
        r = MOL_MASS(6) * ...
    else if (i <= 5) then
        r = ... * MOL_MASS(i) * ...
    end if
end do
```
✓ Indices correctly shifted: C++ `i=0..5` → Fortran `i=1..6`, C++ `i==2` → Fortran `i==3`, C++ `i==5` → Fortran `i==6`

### Polynomial Evaluations

#### CON_C Polynomial
**C++**: `for (sum = CON_C[4], i = 3; i >= 0; i--)` → accesses CON_C[4,3,2,1,0]
**Fortran**: `jr_state%sum = CON_C(5); do i = 4, 1, -1` → accesses CON_C(5,4,3,2,1)
✓ Same 5 elements in same order

#### CON_L Polynomial
**C++**: `for (sum = CON_L[4], i = 3; i >= 0; i--)` → accesses CON_L[4,3,2,1,0]
**Fortran**: `jr_state%sum = CON_L(5); do i = 4, 1, -1` → accesses CON_L(5,4,3,2,1)
✓ Same 5 elements in same order

#### M_CON Polynomial
**C++**: `for (m_poly = M_CON[6], i = 5; i >= 0; i--)` → accesses M_CON[6,5,4,3,2,1,0]
**Fortran**: `m_poly = M_CON(7); do i = 6, 1, -1` → accesses M_CON(7,6,5,4,3,2,1)
✓ Same 7 elements in same order

#### ZETA_CON Polynomial
**C++**: `for (rho_prime = ZETA_CON[6], i=5; i>=0; i--)` → accesses ZETA_CON[6,5,4,3,2,1,0]
**Fortran**: `rho_prime = ZETA_CON(7); do i = 6, 1, -1` → accesses ZETA_CON(7,6,5,4,3,2,1)
✓ Same 7 elements in same order

#### S_CON/S_BETA Access
**C++**: `for (i = 0; i <= 5; i++) { b[i] = S_CON[i] + S_BETA[i] * ... }`
**Fortran**: `do i = 1, 6; b(i) = S_CON(i) + S_BETA(i) * ...; end do`
✓ Same 6 elements

#### CON_DEN Polynomial
**C++**: `for (log_di=CON_DEN[i][6], j=5; j>=0; j--) { log_di = log_di * t_infinity + CON_DEN[i][j]; }`
**Fortran**: `log_di = CON_DEN(i,7); do j = 6, 1, -1; log_di = log_di * t_infinity + CON_DEN(i,j); end do`
✓ Accesses all 7 columns in correct order

## Summary

All array indices have been verified to correctly map from C++ 0-based indexing to Fortran 1-based indexing. The critical bug in `CON_C(4)` → `CON_C(5)` has been fixed.

**Total Arrays Checked**: 9
**Issues Found**: 1 (CON_C indexing in rho_125)
**Issues Fixed**: 1
**Status**: ✓ ALL VERIFIED
