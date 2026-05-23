!------------------------------------------------------------------------------
!>
!  Unit test for Jacchia-Roberts density calculation using override values
!  for space weather parameters instead of reading from a space weather file.
!
!  This demonstrates the use of f107_override, f107a_override, and kp_override
!  to specify constant solar activity conditions.

program test_density_with_overrides
   use jacchia_roberts_kinds, only: ip, dp
   use jacchia_roberts_module, only: jacchia_roberts_type
   implicit none

   type(jacchia_roberts_type) :: jr_model
   integer(ip) :: status, i, j
   real(dp) :: density
   real(dp) :: height, geo_lat, utc_mjd
   real(dp) :: position(3), sun_vector(3)
   real(dp) :: f107, f107a, kp

   ! Test altitudes (km)
   real(dp), parameter :: test_heights(7) = [100.0_dp, 150.0_dp, 200.0_dp, &
                                              300.0_dp, 500.0_dp, 800.0_dp, 1200.0_dp]

   print *, ''
   print *, '=========================================='
   print *, 'JR Density Test with Override Values'
   print *, '=========================================='
   print *, ''
   print *, 'Testing Jacchia-Roberts atmospheric density'
   print *, 'using f107_override, f107a_override, kp_override'
   print *, 'instead of space weather file.'
   print *, ''

   ! Fixed test geometry
   geo_lat = 0.0_dp      ! Equator
   utc_mjd = 55781.5_dp  ! Jan 1, 2012 00:00 UTC
   sun_vector = [1.0_dp, 0.0_dp, 0.0_dp]  ! Sun along X-axis

   ! ===== TEST CASE 1: Solar Minimum Conditions =====
   print *, '=========================================='
   print *, 'TEST CASE 1: Solar Minimum'
   print *, '=========================================='
   print *, ''

   f107 = 70.0_dp
   f107a = 70.0_dp
   kp = 1.0_dp

   write(*, '(A, F6.1)') ' F10.7  daily:   ', f107
   write(*, '(A, F6.1)') ' F10.7a average: ', f107a
   write(*, '(A, F4.1)') ' Kp index:       ', kp
   print *, ''

   ! Initialize with override values
   call jr_model%initialize(6356.766_dp, status, &
                           f107_override=f107, &
                           f107a_override=f107a, &
                           kp_override=kp)

   if (status /= 0) then
      error stop 'ERROR: Failed to initialize JR model'
   end if

   write(*, '(A)') ' Altitude      Density        log10(ρ)'
   write(*, '(A)') ' --------      -------        ---------'

   do i = 1, size(test_heights)
      height = test_heights(i)
      position = [6378.137_dp + height, 0.0_dp, 0.0_dp]

      density = jr_model%density(height, position, sun_vector, geo_lat, utc_mjd)

      write(*, '(F7.1, A, ES14.6, A, F8.4)') &
         height, ' km   ', density, ' kg/m³   ', log10(density)
   end do

   call jr_model%destroy()
   print *, ''

   ! ===== TEST CASE 2: Solar Maximum Conditions =====
   print *, '=========================================='
   print *, 'TEST CASE 2: Solar Maximum'
   print *, '=========================================='
   print *, ''

   f107 = 250.0_dp
   f107a = 250.0_dp
   kp = 5.0_dp

   write(*, '(A, F6.1)') ' F10.7  daily:   ', f107
   write(*, '(A, F6.1)') ' F10.7a average: ', f107a
   write(*, '(A, F4.1)') ' Kp index:       ', kp
   print *, ''

   call jr_model%initialize(6356.766_dp, status, &
                           f107_override=f107, &
                           f107a_override=f107a, &
                           kp_override=kp)

   if (status /= 0) then
      error stop 'ERROR: Failed to initialize JR model'
   end if

   write(*, '(A)') ' Altitude      Density        log10(ρ)'
   write(*, '(A)') ' --------      -------        ---------'

   do i = 1, size(test_heights)
      height = test_heights(i)
      position = [6378.137_dp + height, 0.0_dp, 0.0_dp]

      density = jr_model%density(height, position, sun_vector, geo_lat, utc_mjd)

      write(*, '(F7.1, A, ES14.6, A, F8.4)') &
         height, ' km   ', density, ' kg/m³   ', log10(density)
   end do

   call jr_model%destroy()
   print *, ''

   ! ===== TEST CASE 3: Moderate Conditions =====
   print *, '=========================================='
   print *, 'TEST CASE 3: Moderate Solar Activity'
   print *, '=========================================='
   print *, ''

   f107 = 150.0_dp
   f107a = 150.0_dp
   kp = 3.0_dp

   write(*, '(A, F6.1)') ' F10.7  daily:   ', f107
   write(*, '(A, F6.1)') ' F10.7a average: ', f107a
   write(*, '(A, F4.1)') ' Kp index:       ', kp
   print *, ''

   call jr_model%initialize(6356.766_dp, status, &
                           f107_override=f107, &
                           f107a_override=f107a, &
                           kp_override=kp)

   if (status /= 0) then
      error stop 'ERROR: Failed to initialize JR model'
   end if

   write(*, '(A)') ' Altitude      Density        log10(ρ)'
   write(*, '(A)') ' --------      -------        ---------'

   do i = 1, size(test_heights)
      height = test_heights(i)
      position = [6378.137_dp + height, 0.0_dp, 0.0_dp]

      density = jr_model%density(height, position, sun_vector, geo_lat, utc_mjd)

      write(*, '(F7.1, A, ES14.6, A, F8.4)') &
         height, ' km   ', density, ' kg/m³   ', log10(density)
   end do

   call jr_model%destroy()
   print *, ''

   ! ===== TEST CASE 4: Geomagnetic Storm (High Kp) =====
   print *, '=========================================='
   print *, 'TEST CASE 4: Geomagnetic Storm'
   print *, '=========================================='
   print *, ''

   f107 = 150.0_dp
   f107a = 150.0_dp
   kp = 8.0_dp  ! Severe storm

   write(*, '(A, F6.1)') ' F10.7  daily:   ', f107
   write(*, '(A, F6.1)') ' F10.7a average: ', f107a
   write(*, '(A, F4.1)') ' Kp index:       ', kp
   print *, ''

   call jr_model%initialize(6356.766_dp, status, &
                           f107_override=f107, &
                           f107a_override=f107a, &
                           kp_override=kp)

   if (status /= 0) then
      error stop 'ERROR: Failed to initialize JR model'
   end if

   write(*, '(A)') ' Altitude      Density        log10(ρ)'
   write(*, '(A)') ' --------      -------        ---------'

   do i = 1, size(test_heights)
      height = test_heights(i)
      position = [6378.137_dp + height, 0.0_dp, 0.0_dp]

      density = jr_model%density(height, position, sun_vector, geo_lat, utc_mjd)

      write(*, '(F7.1, A, ES14.6, A, F8.4)') &
         height, ' km   ', density, ' kg/m³   ', log10(density)
   end do

   call jr_model%destroy()

   print *, ''
   print *, 'TEST PASSED'
   print *, ''

end program test_density_with_overrides
