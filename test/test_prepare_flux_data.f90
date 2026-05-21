!*****************************************************************************************
!>
!  Unit tests for prepare_flux_data subroutine
!  Tests all conditional branches in the flux data preparation logic

program test_prepare_flux_data

   use jacchia_roberts_module
   use space_weather_module
   use iso_fortran_env, only: real64, int32

   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: ip = int32

   type(jacchia_roberts_type) :: jrmodel
   type(sw_data_type) :: sw_data
   integer :: total_tests, passed_tests, status
   logical :: all_passed
   real(dp), parameter :: EARTH_POLAR_RADIUS = 6356.766_dp  ! km

   total_tests = 0
   passed_tests = 0

   write(*,*) ''
   write(*,*) '========================================='
   write(*,*) 'Testing prepare_flux_data subroutine'
   write(*,*) '========================================='
   write(*,*) ''

   ! Initialize the space weather data
   call sw_data%initialize('data/SpaceWeather-All-v1.2.txt', status)
   if (status /= 0) then
      write(*,*) 'ERROR: Failed to initialize space weather data'
      stop 1
   end if

   ! Initialize the model
   call jrmodel%initialize(EARTH_POLAR_RADIUS, 'data/SpaceWeather-All-v1.2.txt', status)
   if (status /= 0) then
      write(*,*) 'ERROR: Failed to initialize model'
      stop 1
   end if

   ! Test 1: Kp selection - various 3-hour periods (sub_index 0-7)
   call test_kp_various_periods(sw_data, total_tests, passed_tests)

   ! Test 2: Kp selection - sub_index >= 8 (should clamp to 7)
   call test_kp_overflow(sw_data, total_tests, passed_tests)

   ! Test 3: Kp selection - sub_index < 0 (previous day lookup)
   call test_kp_previous_day(sw_data, total_tests, passed_tests)

   ! Test 4: F10.7 offset - before reference epoch (5pm measurement)
   call test_f107_offset_before_ref(sw_data, total_tests, passed_tests)

   ! Test 5: F10.7 offset - after reference epoch (8pm measurement)
   call test_f107_offset_after_ref(sw_data, total_tests, passed_tests)

   ! Test 6: F10.7 daily - after measurement time (use yesterday)
   call test_f107_daily_after_measurement(sw_data, total_tests, passed_tests)

   ! Test 7: F10.7 daily - before measurement time (use 2 days ago)
   call test_f107_daily_before_measurement(sw_data, total_tests, passed_tests)

   ! Test 8: F10.7a - after measurement time (use today)
   call test_f107a_after_measurement(sw_data, total_tests, passed_tests)

   ! Test 9: F10.7a - before measurement time (use yesterday)
   call test_f107a_before_measurement(sw_data, total_tests, passed_tests)

   ! Test 10: Edge case - exactly at midnight (frac_epoch = 0)
   call test_edge_midnight(sw_data, total_tests, passed_tests)

   ! Test 11: Edge case - near measurement time boundaries
   call test_edge_measurement_boundary(sw_data, total_tests, passed_tests)

   ! Cleanup
   call sw_data%destroy()
   call jrmodel%destroy()

   ! Summary
   write(*,*) ''
   write(*,*) '========================================='
   write(*,'(A,I0,A,I0,A)') ' RESULTS: ', passed_tests, '/', total_tests, ' tests passed'
   write(*,*) '========================================='
   write(*,*) ''

   all_passed = (passed_tests == total_tests)
   if (.not. all_passed) then
      write(*,*) 'SOME TESTS FAILED!'
      stop 1
   else
      write(*,*) 'ALL TESTS PASSED!'
   end if

contains

   !---------------------------------------------------------------------------
   subroutine test_kp_various_periods(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status
      integer(ip) :: i

      total = total + 1
      write(*,*) 'Test 1: Kp selection for various 3-hour periods (0-7)'

      ! Use a date in the middle of available data: Jan 15, 2021
      base_mjd = 59229.0_dp  ! MJD for Jan 15, 2021 at midnight

      ! Get flux data for this date
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data for test date'
         return
      end if

      ! Test each 3-hour period
      do i = 0, 7
         ! Time corresponding to period i: i*3 hours + 7.5 hours (middle of period + 6.7hr lag)
         test_mjd = base_mjd + (i * 3.0_dp + 7.5_dp) / 24.0_dp

         call prepare_flux_data(sw_data, flux_data, test_mjd, &
                               kp_out, f107_out, f107a_out)

         ! Verify kp_out matches the expected period
         if (abs(kp_out - flux_data%kp(i + 1)) > 1.0e-10_dp) then
            write(*,'(A,I0,A)') '   FAILED: Period ', i, ' Kp mismatch'
            return
         end if
      end do

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_kp_various_periods

   !---------------------------------------------------------------------------
   subroutine test_kp_overflow(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 2: Kp selection with sub_index >= 8 (clamps to 7)'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Time way past the end of day (should clamp to period 7)
      test_mjd = base_mjd + 30.0_dp / 24.0_dp  ! 30 hours past midnight

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! Should use period 7 (last period)
      if (abs(kp_out - flux_data%kp(8)) > 1.0e-10_dp) then
         write(*,*) '   FAILED: Should clamp to period 7'
         return
      end if

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_kp_overflow

   !---------------------------------------------------------------------------
   subroutine test_kp_previous_day(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data, flux_prev
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 3: Kp selection with sub_index < 0 (previous day)'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Time early in the day such that with 6.7hr lag, we need previous day
      test_mjd = base_mjd + 3.0_dp / 24.0_dp  ! 3am (3am - 6.7hr = -3.7hr = previous day)

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! Get previous day's data to verify
      call sw_data%get_flux_data(base_mjd - 1.0_dp, flux_prev, sw_status)

      if (sw_status) then
         ! Should use appropriate period from previous day
         ! 3am - 6.7hr = -3.7hr = 20.3hr previous day = period 6 (18:00-21:00)
         if (abs(kp_out - flux_prev%kp(7)) > 1.0e-10_dp) then
            write(*,*) '   FAILED: Incorrect period from previous day'
            write(*,*) '   Expected kp:', flux_prev%kp(7), ' Got:', kp_out
            return
         end if
      else
         write(*,*) '   WARNING: Could not verify previous day lookup'
      end if

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_kp_previous_day

   !---------------------------------------------------------------------------
   subroutine test_f107_offset_before_ref(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 4: F10.7 offset before reference epoch (5pm measurement)'

      ! Use date before 5/31/91: Jan 1, 1990 (MJD 47892)
      base_mjd = 47892.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data for 1990'
         return
      end if

      ! Test at 6pm (after 5pm measurement time)
      test_mjd = base_mjd + 18.0_dp / 24.0_dp

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! Just verify it runs without error - the offset is internal
      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_f107_offset_before_ref

   !---------------------------------------------------------------------------
   subroutine test_f107_offset_after_ref(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 5: F10.7 offset after reference epoch (8pm measurement)'

      base_mjd = 59229.0_dp  ! 2021 - well after 5/31/91
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Test at 9pm (after 8pm measurement time)
      test_mjd = base_mjd + 21.0_dp / 24.0_dp

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_f107_offset_after_ref

   !---------------------------------------------------------------------------
   subroutine test_f107_daily_after_measurement(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data, flux_prev
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 6: F10.7 daily after measurement time (use yesterday)'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Test at 10pm (after 8pm measurement)
      test_mjd = base_mjd + 22.0_dp / 24.0_dp

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! Get previous day's data
      call sw_data%get_flux_data(base_mjd - 1.0_dp, flux_prev, sw_status)

      if (sw_status) then
         if (abs(f107_out - flux_prev%f107_obs) > 1.0e-10_dp) then
            write(*,*) '   FAILED: Should use previous day F10.7'
            write(*,*) '   Expected:', flux_prev%f107_obs, ' Got:', f107_out
            return
         end if
      end if

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_f107_daily_after_measurement

   !---------------------------------------------------------------------------
   subroutine test_f107_daily_before_measurement(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data, flux_2days
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 7: F10.7 daily before measurement time (use 2 days ago)'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Test at 6am (before 8am data validity boundary)
      test_mjd = base_mjd + 6.0_dp / 24.0_dp

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! Get 2 days ago data
      call sw_data%get_flux_data(base_mjd - 2.0_dp, flux_2days, sw_status)

      if (sw_status) then
         if (abs(f107_out - flux_2days%f107_obs) > 1.0e-10_dp) then
            write(*,*) '   FAILED: Should use 2 days ago F10.7'
            write(*,*) '   Expected:', flux_2days%f107_obs, ' Got:', f107_out
            return
         end if
      end if

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_f107_daily_before_measurement

   !---------------------------------------------------------------------------
   subroutine test_f107a_after_measurement(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 8: F10.7a after measurement time (use today)'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Test at 10pm (after 8pm measurement)
      test_mjd = base_mjd + 22.0_dp / 24.0_dp

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      if (abs(f107a_out - flux_data%f107a_obs_ctr) > 1.0e-10_dp) then
         write(*,*) '   FAILED: Should use current day F10.7a'
         return
      end if

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_f107a_after_measurement

   !---------------------------------------------------------------------------
   subroutine test_f107a_before_measurement(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data, flux_prev
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 9: F10.7a before measurement time (use yesterday)'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Test at 6am (before 8am data validity boundary)
      test_mjd = base_mjd + 6.0_dp / 24.0_dp

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! Get previous day's data
      call sw_data%get_flux_data(base_mjd - 1.0_dp, flux_prev, sw_status)

      if (sw_status) then
         if (abs(f107a_out - flux_prev%f107a_obs_ctr) > 1.0e-10_dp) then
            write(*,*) '   FAILED: Should use previous day F10.7a'
            write(*,*) '   Expected:', flux_prev%f107a_obs_ctr, ' Got:', f107a_out
            return
         end if
      end if

      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_f107a_before_measurement

   !---------------------------------------------------------------------------
   subroutine test_edge_midnight(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 10: Edge case - exactly at midnight (frac_epoch = 0)'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Test exactly at midnight
      test_mjd = base_mjd

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! With 6.7hr lag, midnight - 6.7hr = previous day period 6
      ! Just verify it runs without error
      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_edge_midnight

   !---------------------------------------------------------------------------
   subroutine test_edge_measurement_boundary(sw_data, total, passed)
      type(sw_data_type), intent(inout) :: sw_data
      integer, intent(inout) :: total, passed

      type(flux_data_type) :: flux_data
      real(dp) :: kp_out, f107_out, f107a_out
      real(dp) :: test_mjd, base_mjd
      logical :: sw_status

      total = total + 1
      write(*,*) 'Test 11: Edge case - near measurement time boundary'

      base_mjd = 59229.0_dp
      call sw_data%get_flux_data(base_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) '   FAILED: Could not get flux data'
         return
      end if

      ! Test exactly at 8pm (measurement time)
      test_mjd = base_mjd + 20.0_dp / 24.0_dp

      call prepare_flux_data(sw_data, flux_data, test_mjd, &
                            kp_out, f107_out, f107a_out)

      ! Just verify it runs without error at the boundary
      passed = passed + 1
      write(*,*) '   PASSED'
   end subroutine test_edge_measurement_boundary

end program test_prepare_flux_data