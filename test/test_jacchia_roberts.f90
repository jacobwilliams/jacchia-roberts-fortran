!------------------------------------------------------------------------------
! Program: test_jacchia_roberts
!------------------------------------------------------------------------------
! Description:
!   Test program for the Jacchia-Roberts atmosphere model
!
! This program demonstrates how to:
!   1. Initialize the module
!   2. Set up input parameters
!   3. Call the density calculation
!   4. Display results
!------------------------------------------------------------------------------

program test_jacchia_roberts
   use jacchia_roberts_module
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none

   ! Kind parameter
   integer, parameter :: dp = real64

   ! Constants
   real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp
   real(dp), parameter :: RAD_PER_DEG = PI / 180.0_dp

   ! Variables
   real(dp) :: density
   real(dp) :: height
   real(dp) :: position(3)
   real(dp) :: sun_vector(3)
   real(dp) :: geo_lat
   real(dp) :: sun_dec
   real(dp) :: utc_mjd
   type(geoparms_type) :: geo

   ! Earth polar radius (km)
   real(dp), parameter :: EARTH_POLAR_RADIUS = 6356.766_dp

   ! Test parameters
   integer :: i
   real(dp) :: test_heights(5)

   ! Banner
   write(*,'(A)') ''
   write(*,'(A)') '======================================================'
   write(*,'(A)') '  Jacchia-Roberts Atmosphere Model - Test Program'
   write(*,'(A)') '======================================================'
   write(*,'(A)') ''

   ! Initialize the module
   call jr_init(EARTH_POLAR_RADIUS)
   write(*,'(A)') 'Module initialized with Earth polar radius = 6356.766 km'
   write(*,'(A)') ''

   ! Set up test scenario
   ! Position: 7000 km from Earth center (approximately 650 km altitude)
   position(1) = 7000.0_dp
   position(2) = 0.0_dp
   position(3) = 0.0_dp

   ! Sun vector (unit vector pointing to sun)
   ! Test collinear case: sun and spacecraft both on x-axis
   sun_vector(1) = 1.0_dp
   sun_vector(2) = 0.0_dp
   sun_vector(3) = 0.0_dp

   ! Geodetic latitude (degrees)
   geo_lat = 0.0_dp

   ! Sun declination (radians) - equinox
   sun_dec = 0.0_dp

   ! UTC Modified Julian Date (arbitrary date)
   utc_mjd = 51544.5_dp  ! Jan 1, 2000 12:00 UTC

   ! Geomagnetic parameters
   geo%tkp = 3.0_dp      ! Moderate geomagnetic activity
   ! Calculate exospheric temperature from F10.7 and F10.7a
   ! xtemp = 379 + 3.24 * F10.7a + 1.3 * (F10.7 - F10.7a)
   ! Using nominal values: F10.7 = 150, F10.7a = 150
   geo%xtemp = 379.0_dp + 3.24_dp * 150.0_dp + 1.3_dp * (150.0_dp - 150.0_dp)
   geo%xtemp = 865.0_dp  ! Simplified calculation

   write(*,'(A)') 'Test Scenario:'
   write(*,'(A,F10.2,A)') '  Position magnitude:     ', sqrt(sum(position**2)), ' km'
   write(*,'(A,F10.4)') '  Geodetic latitude:      ', geo_lat
   write(*,'(A,F10.4)') '  Sun declination:        ', sun_dec
   write(*,'(A,F12.2)') '  UTC MJD:                ', utc_mjd
   write(*,'(A,F8.2)') '  Geomagnetic Kp:         ', geo%tkp
   write(*,'(A,F8.2)') '  Exospheric temp:        ', geo%xtemp
   write(*,'(A)') ''

   ! Test at various altitudes
   test_heights = [150.0_dp, 250.0_dp, 400.0_dp, 600.0_dp, 1000.0_dp]

   write(*,'(A)') 'Density Calculations:'
   write(*,'(A)') '------------------------------------------------------'
   write(*,'(A)') '  Altitude (km)    Density (kg/m^3)    Density (g/cm^3)'
   write(*,'(A)') '------------------------------------------------------'

   do i = 1, size(test_heights)
      height = test_heights(i)

      ! Calculate density
      density = jacchia_roberts_density(height, position, sun_vector, &
                                       geo_lat, sun_dec, geo, utc_mjd)

      ! Display results
      write(*,'(F12.2,8X,ES14.6,4X,ES14.6)') height, density, density * 1.0e-3_dp
   end do

   write(*,'(A)') '------------------------------------------------------'
   write(*,'(A)') ''

   ! Additional test: Density profile from 100 to 1000 km
   write(*,'(A)') 'Generating density profile (100-1000 km)...'
   write(*,'(A)') ''

   open(unit=10, file='density_profile.dat', status='replace', action='write')
   write(10,'(A)') '# Altitude (km), Density (kg/m^3), Density (g/cm^3)'

   do i = 100, 1000, 10
      height = real(i, dp)
      density = jacchia_roberts_density(height, position, sun_vector, &
                                       geo_lat, sun_dec, geo, utc_mjd)
      write(10,'(F12.2,2(2X,ES16.8))') height, density, density * 1.0e-3_dp
   end do

   close(10)

   write(*,'(A)') 'Density profile written to: density_profile.dat'
   write(*,'(A)') ''

   write(*,'(A)') '======================================================'
   write(*,'(A)') '  Test completed successfully!'
   write(*,'(A)') '======================================================'
   write(*,'(A)') ''

end program test_jacchia_roberts
