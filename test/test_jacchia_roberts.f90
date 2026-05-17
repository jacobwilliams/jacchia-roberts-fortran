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
   use jacchia_roberts_module, only: jacchia_roberts_type
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   ! Kind parameter
   integer, parameter :: dp = real64
   integer, parameter :: ip = int32

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
   integer(ip) :: sw_status
   type(jacchia_roberts_type) :: atm

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
   call atm%initialize(EARTH_POLAR_RADIUS, 'data/SpaceWeather-All-v1.2.txt', sw_status)
   write(*,'(A)') 'Module initialized with Earth polar radius = 6356.766 km'
   if (sw_status /= 0) then
      write(*,'(A)') 'WARNING: Could not load space weather file'
      write(*,'(A)') '         Using nominal values (F10.7=150, Kp=3)'
   end if
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

   ! UTC Modified Julian Date (Jan 1, 2021 12:00 UTC)
   ! This date is in the space weather file
   utc_mjd = 59215.5_dp

   write(*,'(A)') 'Test Scenario:'
   write(*,'(A,F10.2,A)') '  Position magnitude:     ', sqrt(sum(position**2)), ' km'
   write(*,'(A,F10.4)') '  Geodetic latitude:      ', geo_lat
   write(*,'(A,F10.4)') '  Sun declination:        ', sun_dec
   write(*,'(A,F12.2)') '  UTC MJD:                ', utc_mjd
   write(*,'(A)') '  (Space weather retrieved from data file)'
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
      density = atm%density(height, position, sun_vector, &
                            geo_lat, sun_dec, utc_mjd)

      ! Display results
      write(*,'(F12.2,8X,ES14.6,4X,ES14.6)') height, density, density * 1.0e-3_dp
   end do

   write(*,'(A)') '------------------------------------------------------'
   write(*,'(A)') ''

   ! test different sun vectors (general, collinear, orthogonal, anti-collinear)
   write(*,'(A)') 'Testing different sun vector orientations...'

   ! Collinear case (already tested above)
   write(*,'(A)') '  Collinear case: sun and spacecraft both on x-axis (already tested)'

   ! Orthogonal case: sun on y-axis, spacecraft on x-axis
   sun_vector(1) = 0.0_dp
   sun_vector(2) = 1.0_dp
   sun_vector(3) = 0.0_dp
   density = atm%density(400.0_dp, position, sun_vector, &
                           geo_lat, sun_dec, utc_mjd)
   write(*,'(A)') '  Orthogonal case: sun on y-axis, spacecraft on x-axis'
   write(*,'(A,F12.2,8X,ES14.6,4X,ES14.6)') '    Altitude: 400.00 km', density, density * 1.0e-3_dp

   ! Anti-collinear case: sun on negative x-axis, spacecraft on positive x-axis
   sun_vector(1) = -1.0_dp
   sun_vector(2) = 0.0_dp
   sun_vector(3) = 0.0_dp
   density = atm%density(400.0_dp, position, sun_vector, &
                           geo_lat, sun_dec, utc_mjd)
   write(*,'(A)') '  Anti-collinear case: sun on negative x-axis, spacecraft on positive x-axis'
   write(*,'(A,F12.2,8X,ES14.6,4X,ES14.6)') '    Altitude: 400.00 km', density, density * 1.0e-3_dp

   ! General case: sun vector at 45 degrees in x-y plane
   sun_vector(1) = 1.0_dp / sqrt(2.0_dp)
   sun_vector(2) = 1.0_dp / sqrt(2.0_dp)
   sun_vector(3) = 0.0_dp
   density = atm%density(400.0_dp, position, sun_vector, &
                           geo_lat, sun_dec, utc_mjd)
   write(*,'(A)') '  General case: sun vector at 45 degrees in x-y plane'
   write(*,'(A,F12.2,8X,ES 14.6,4X,ES14.6)') '    Altitude: 400.00 km', density, density * 1.0e-3_dp



   write(*,'(A)') '------------------------------------------------------'
   write(*,'(A)') ''

   ! Additional test: Density profile from 90 to 2510 km
   write(*,'(A)') 'Generating density profile (90-2510 km)...'
   write(*,'(A)') ''

   open(unit=10, file='density_profile.dat', status='replace', action='write')
   write(10,'(A)') '# Altitude (km), Density (kg/m^3), Density (g/cm^3)'

   do i = 90, 2510, 5
      height = real(i, dp)
      density = atm%density(height, position, sun_vector, &
                            geo_lat, sun_dec, utc_mjd)
      write(10,'(F12.2,2(2X,ES16.8))') height, density, density * 1.0e-3_dp
   end do

   close(10)

   write(*,'(A)') 'Density profile written to: density_profile.dat'
   write(*,'(A)') ''

   ! Clean up
   call atm%destroy()

   write(*,'(A)') '======================================================'
   write(*,'(A)') '  Test completed successfully!'
   write(*,'(A)') '======================================================'
   write(*,'(A)') ''

end program test_jacchia_roberts
