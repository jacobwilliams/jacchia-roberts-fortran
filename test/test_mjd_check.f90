program test_mjd_check
   use jacchia_roberts_kinds, only: ip, dp
   use jacchia_roberts_utilities, only: date_to_mjd
   implicit none

   real(dp) :: mjd_result

   ! Test Jan 1, 2000 - should be 51544.0
   call date_to_mjd(2000, 1, 1, mjd_result)
   print '(A,F20.10)', 'Jan 1, 2000: ', mjd_result
   if (mjd_result /= 51544.0_dp) then
      error stop 'Error: Expected MJD 51544.0 for Jan 1, 2000'
   end if

   ! Test with the C++ formula manually
   call test_cpp_formula(2000, 1, 1)

contains

   subroutine test_cpp_formula(year, month, day)
      integer :: year, month, day
      integer :: computeYearMon, computeMonth
      real(dp) :: ModJulianDay

      computeYearMon = (7*(year + (month + 9)/12))/4
      computeMonth = (275 * month)/9
      ModJulianDay = 367.0_dp*year - real(computeYearMon,dp) + real(computeMonth,dp) + &
                     real(day,dp) + 1721013.5_dp - 2400000.5_dp

      print '(A,F20.10)', 'C++ formula:  ', ModJulianDay
      print '(A,I10)', 'computeYearMon: ', computeYearMon
      print '(A,I10)', 'computeMonth:   ', computeMonth

      if (ModJulianDay /= 51544.0_dp) then
         error stop 'Error: Expected MJD 51544.0 for Jan 1, 2000'
      end if

   end subroutine test_cpp_formula

end program test_mjd_check
