!------------------------------------------------------------------------------
! Program: test_space_weather
!------------------------------------------------------------------------------
! Description:
!   Test program for the space weather data module
!
!   Tests:
!   1. Reading CSSI space weather file
!   2. Retrieving data from daily section
!   3. Retrieving data from monthly section (lines 311-322 in space_weather_module)
!   4. Boundary conditions
!------------------------------------------------------------------------------

program test_space_weather
   use space_weather_module, only: sw_data_type, flux_data_type
   use jacchia_roberts_kinds, only: ip, dp
   implicit none

   type(sw_data_type) :: sw_data
   type(flux_data_type) :: flux_data
   integer(ip) :: init_status
   logical :: status
   integer :: test_count, pass_count
   real(dp) :: mjd_test
   logical :: test_passed

   ! Initialize counters
   test_count = 0
   pass_count = 0

   write(*,'(A)') ''
   write(*,'(A)') '=========================================='
   write(*,'(A)') '  Space Weather Module - Test Suite'
   write(*,'(A)') '=========================================='
   write(*,'(A)') ''

   ! Test 1: Initialize with space weather file
   write(*,'(A)') 'Test 1: Initialize with space weather file'
   call sw_data%initialize('data/SpaceWeather-All-v1.2.txt', init_status)
   test_count = test_count + 1
   if (init_status == 0) then
      write(*,'(A)') '  ✓ Initialization successful'
      pass_count = pass_count + 1
   else
      write(*,'(A)') '  ✗ Initialization failed'
      stop 1
   end if
   write(*,'(A)') ''

   ! Test 2: Retrieve data from daily section (Jan 1, 2021)
   write(*,'(A)') 'Test 2: Retrieve data from daily section'
   mjd_test = 59215.5_dp  ! Jan 1, 2021 noon
   call sw_data%get_flux_data(mjd_test, flux_data, status)
   test_count = test_count + 1
   test_passed = (status .and. flux_data%f107_obs > 0.0_dp)
   if (test_passed) then
      write(*,'(A,F10.2)') '  ✓ Retrieved daily data for MJD ', mjd_test
      write(*,'(A,F8.2)') '    F10.7 obs: ', flux_data%f107_obs
      pass_count = pass_count + 1
   else
      write(*,'(A)') '  ✗ Failed to retrieve daily data'
   end if
   write(*,'(A)') ''

   ! Test 3: Retrieve data from monthly section (Aug 1, 2025)
   write(*,'(A)') 'Test 3: Retrieve data from monthly section'
   write(*,'(A)') '  (This tests lines 311-322 in space_weather_module.f90)'
   ! Aug 1, 2025: daily data ends May 4, 2025 (≈MJD 60800.5)
   ! Aug 1 is ~88 days after May 4
   mjd_test = 60888.5_dp
   call sw_data%get_flux_data(mjd_test, flux_data, status)
   test_count = test_count + 1
   test_passed = (status .and. flux_data%f107_obs > 0.0_dp)
   if (test_passed) then
      write(*,'(A,F10.2)') '  ✓ Retrieved monthly data for MJD ', mjd_test
      write(*,'(A,F8.2)') '    F10.7 obs: ', flux_data%f107_obs
      write(*,'(A,F8.2)') '    F10.7a ctr: ', flux_data%f107a_obs_ctr
      pass_count = pass_count + 1
   else
      write(*,'(A)') '  ✗ Failed to retrieve monthly data'
   end if
   write(*,'(A)') ''

   ! Test 4: Retrieve data between monthly points (Aug 15, 2025)
   write(*,'(A)') 'Test 4: Retrieve data between monthly points'
   mjd_test = 60902.5_dp  ! Aug 15, 2025 noon (between Aug 1 and Sep 1)
   call sw_data%get_flux_data(mjd_test, flux_data, status)
   test_count = test_count + 1
   test_passed = (status .and. flux_data%f107_obs > 0.0_dp)
   if (test_passed) then
      write(*,'(A,F10.2)') '  ✓ Retrieved interpolated monthly data for MJD ', mjd_test
      write(*,'(A,F8.2)') '    F10.7 obs: ', flux_data%f107_obs
      pass_count = pass_count + 1
   else
      write(*,'(A)') '  ✗ Failed to retrieve interpolated monthly data'
   end if
   write(*,'(A)') ''

   ! Test 5: Retrieve data near end of monthly section (Dec 2026)
   write(*,'(A)') 'Test 5: Retrieve data near end of monthly section'
   mjd_test = 61376.5_dp  ! Dec 1, 2026 noon
   call sw_data%get_flux_data(mjd_test, flux_data, status)
   test_count = test_count + 1
   test_passed = (status .and. flux_data%f107_obs > 0.0_dp)
   if (test_passed) then
      write(*,'(A,F10.2)') '  ✓ Retrieved data near end of monthly section for MJD ', mjd_test
      write(*,'(A,F8.2)') '    F10.7 obs: ', flux_data%f107_obs
      pass_count = pass_count + 1
   else
      write(*,'(A)') '  ✗ Failed to retrieve data near end'
   end if
   write(*,'(A)') ''

   ! Test 6: Retrieve data before file start (should use first record)
   write(*,'(A)') 'Test 6: Retrieve data before file start'
   mjd_test = 30000.0_dp  ! Way before file starts (1957)
   call sw_data%get_flux_data(mjd_test, flux_data, status)
   test_count = test_count + 1
   if (status) then
      write(*,'(A,F10.2)') '  ✓ Retrieved first record for MJD ', mjd_test
      pass_count = pass_count + 1
   else
      write(*,'(A)') '  ✗ Failed to retrieve data before start'
   end if
   write(*,'(A)') ''

   ! Test 7: Retrieve data after file end (should use last record)
   write(*,'(A)') 'Test 7: Retrieve data after file end'
   mjd_test = 90000.0_dp  ! Way after file ends
   call sw_data%get_flux_data(mjd_test, flux_data, status)
   test_count = test_count + 1
   if (status) then
      write(*,'(A,F10.2)') '  ✓ Retrieved last record for MJD ', mjd_test
      pass_count = pass_count + 1
   else
      write(*,'(A)') '  ✗ Failed to retrieve data after end'
   end if
   write(*,'(A)') ''

   ! Cleanup
   call sw_data%destroy()

   ! Print summary
   write(*,'(A)') ''
   write(*,'(A)') '=========================================='
   write(*,'(A)') 'Test Summary:'
   write(*,'(A)') '=========================================='
   write(*,'(A,I0)') '  Total tests:  ', test_count
   write(*,'(A,I0)') '  Passed:       ', pass_count
   write(*,'(A,I0)') '  Failed:       ', test_count - pass_count
   write(*,'(A)') ''

   if (pass_count == test_count) then
      write(*,'(A)') 'All tests PASSED!'
      stop 0
   else
      write(*,'(A)') 'Some tests FAILED!'
      stop 1
   end if

end program test_space_weather
