!------------------------------------------------------------------------------
!>
!   Unit tests for jacchia_roberts_utilities module
!
!   Tests the mathematical utility functions used in the Jacchia-Roberts
!   atmospheric density model.

program test_jacchia_roberts_utilities
   use jacchia_roberts_kinds, only: ip, dp
   use jacchia_roberts_utilities, only: roots, deflate_polynomial, date_to_mjd
   implicit none

   integer :: n_tests, n_passed, n_failed
   real(dp), parameter :: tol = 1.0e-10_dp

   n_tests = 0
   n_passed = 0
   n_failed = 0

   print *, '=========================================='
   print *, 'Testing jacchia_roberts_utilities module'
   print *, '=========================================='
   print *, ''

   ! Test date_to_mjd
   call test_date_to_mjd_known_dates(n_tests, n_passed, n_failed)
   call test_date_to_mjd_leap_years(n_tests, n_passed, n_failed)

   ! Test roots function
   call test_roots_quadratic(n_tests, n_passed, n_failed)
   call test_roots_cubic(n_tests, n_passed, n_failed)

   ! Test deflate_polynomial
   call test_deflate_polynomial_simple(n_tests, n_passed, n_failed)
   call test_deflate_polynomial_quadratic(n_tests, n_passed, n_failed)

   ! Print summary
   print *, ''
   print *, '=========================================='
   print *, 'Test Summary:'
   print *, '=========================================='
   print '(A,I0)', '  Total tests:  ', n_tests
   print '(A,I0)', '  Passed:       ', n_passed
   print '(A,I0)', '  Failed:       ', n_failed
   print *, ''

   if (n_failed == 0) then
      print *, 'All tests PASSED!'
   else
      error stop 'Some tests FAILED!'
   end if

contains

   !---------------------------------------------------------------------------
   !>
   !   Test date_to_mjd with known reference dates
   subroutine test_date_to_mjd_known_dates(n_tests, n_passed, n_failed)
      integer, intent(inout) :: n_tests, n_passed, n_failed
      real(dp) :: mjd_result, mjd_expected
      logical :: test_pass

      print *, 'Testing date_to_mjd with known dates...'

      ! Test 1: January 1, 2000 (J2000 epoch)
      n_tests = n_tests + 1
      call date_to_mjd(2000, 1, 1, mjd_result)
      mjd_expected = 51544.0_dp
      test_pass = abs(mjd_result - mjd_expected) < tol
      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Jan 1, 2000 (J2000 epoch)'
      else
         n_failed = n_failed + 1
         print '(A,F15.6,A,F15.6)', '  ✗ Jan 1, 2000: Expected ', mjd_expected, ', Got ', mjd_result
      end if

      ! Test 2: January 1, 1957 (Sputnik era)
      n_tests = n_tests + 1
      call date_to_mjd(1957, 10, 1, mjd_result)
      mjd_expected = 36112.0_dp
      test_pass = abs(mjd_result - mjd_expected) < tol
      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Oct 1, 1957 (Sputnik era)'
      else
         n_failed = n_failed + 1
         print '(A,F15.6,A,F15.6)', '  ✗ Oct 1, 1957: Expected ', mjd_expected, ', Got ', mjd_result
      end if

      ! Test 3: January 1, 2021 (test case from space weather)
      n_tests = n_tests + 1
      call date_to_mjd(2021, 1, 1, mjd_result)
      mjd_expected = 59215.0_dp
      test_pass = abs(mjd_result - mjd_expected) < tol
      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Jan 1, 2021'
      else
         n_failed = n_failed + 1
         print '(A,F15.6,A,F15.6)', '  ✗ Jan 1, 2021: Expected ', mjd_expected, ', Got ', mjd_result
      end if

      ! Test 4: May 24, 1968 (Modified Julian Date epoch)
      n_tests = n_tests + 1
      call date_to_mjd(1858, 11, 17, mjd_result)
      mjd_expected = 0.0_dp
      test_pass = abs(mjd_result - mjd_expected) < tol
      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Nov 17, 1858 (MJD epoch)'
      else
         n_failed = n_failed + 1
         print '(A,F15.6,A,F15.6)', '  ✗ Nov 17, 1858: Expected ', mjd_expected, ', Got ', mjd_result
      end if

      print *, ''
   end subroutine test_date_to_mjd_known_dates

   !---------------------------------------------------------------------------
   !>
   !   Test date_to_mjd with leap years
   subroutine test_date_to_mjd_leap_years(n_tests, n_passed, n_failed)
      integer, intent(inout) :: n_tests, n_passed, n_failed
      real(dp) :: mjd1, mjd2
      logical :: test_pass

      print *, 'Testing date_to_mjd with leap years...'

      ! Test: Feb 29, 2000 should be valid (leap year)
      n_tests = n_tests + 1
      call date_to_mjd(2000, 2, 29, mjd1)
      call date_to_mjd(2000, 3, 1, mjd2)
      test_pass = abs(mjd2 - mjd1 - 1.0_dp) < tol
      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Feb 29, 2000 (leap year)'
      else
         n_failed = n_failed + 1
         print '(A)', '  ✗ Feb 29, 2000 leap year test failed'
      end if

      ! Test: 2000 is divisible by 400 (leap year)
      n_tests = n_tests + 1
      call date_to_mjd(2000, 1, 1, mjd1)
      call date_to_mjd(2001, 1, 1, mjd2)
      test_pass = abs(mjd2 - mjd1 - 366.0_dp) < tol
      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Year 2000 has 366 days'
      else
         n_failed = n_failed + 1
         print '(A,F15.6)', '  ✗ Year 2000 day count: ', mjd2 - mjd1
      end if

      ! Test: 1900 is not a leap year (divisible by 100 but not 400)
      n_tests = n_tests + 1
      call date_to_mjd(1900, 1, 1, mjd1)
      call date_to_mjd(1901, 1, 1, mjd2)
      test_pass = abs(mjd2 - mjd1 - 365.0_dp) < tol
      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Year 1900 has 365 days (not leap)'
      else
         n_failed = n_failed + 1
         print '(A,F15.6)', '  ✗ Year 1900 day count: ', mjd2 - mjd1
      end if

      print *, ''
   end subroutine test_date_to_mjd_leap_years

   !---------------------------------------------------------------------------
   !>
   !   Test roots function with quadratic equation
   subroutine test_roots_quadratic(n_tests, n_passed, n_failed)
      integer, intent(inout) :: n_tests, n_passed, n_failed
      real(dp) :: a(3), a_deflated(3), croots(1,2), root1, root2
      logical :: test_pass

      print *, 'Testing roots with quadratic equation...'

      ! Test: x^2 - 5x + 6 = 0, roots are x=2 and x=3
      ! Coefficients: a[0]=6, a[1]=-5, a[2]=1 (lowest power first)
      ! Note: roots() finds ONE root at a time, as used in jacchia_roberts_module
      n_tests = n_tests + 1
      a(1) = 6.0_dp
      a(2) = -5.0_dp
      a(3) = 1.0_dp

      ! Find first root with initial guess
      croots(1,1) = 2.5_dp
      croots(1,2) = 0.0_dp
      call roots(a, 3, croots, 1)
      root1 = croots(1,1)

      ! Deflate polynomial
      call deflate_polynomial(a, 3, root1, a_deflated)

      ! Find second root from deflated polynomial
      croots(1,1) = 3.5_dp
      croots(1,2) = 0.0_dp
      call roots(a_deflated, 2, croots, 1)
      root2 = croots(1,1)

      ! Check if roots are close to 2 and 3 (in any order)
      test_pass = (abs(root1 - 2.0_dp) < 1.0e-6_dp .and. abs(root2 - 3.0_dp) < 1.0e-6_dp) .or. &
                  (abs(root1 - 3.0_dp) < 1.0e-6_dp .and. abs(root2 - 2.0_dp) < 1.0e-6_dp)

      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Quadratic x² - 5x + 6 = 0'
         print '(A,F12.6)', '    Root 1: ', root1
         print '(A,F12.6)', '    Root 2: ', root2
      else
         n_failed = n_failed + 1
         print '(A)', '  ✗ Quadratic roots incorrect'
         print '(A,F12.6,A,F12.6)', '    Got: ', root1, ', ', root2
      end if

      print *, ''
   end subroutine test_roots_quadratic

   !---------------------------------------------------------------------------
   !>
   !   Test roots function with cubic equation
   subroutine test_roots_cubic(n_tests, n_passed, n_failed)
      integer, intent(inout) :: n_tests, n_passed, n_failed
      real(dp) :: a(4), a_def1(4), a_def2(4), croots(1,2), root1, root2, root3
      logical :: test_pass

      print *, 'Testing roots with cubic equation...'

      ! Test: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6 = 0
      ! Coefficients: a[0]=-6, a[1]=11, a[2]=-6, a[3]=1
      ! Find roots one at a time with deflation between (as in real usage)
      n_tests = n_tests + 1
      a(1) = -6.0_dp
      a(2) = 11.0_dp
      a(3) = -6.0_dp
      a(4) = 1.0_dp

      ! Find first root
      croots(1,1) = 0.5_dp
      croots(1,2) = 0.0_dp
      call roots(a, 4, croots, 1)
      root1 = croots(1,1)

      ! Deflate and find second root
      call deflate_polynomial(a, 4, root1, a_def1)
      croots(1,1) = 2.0_dp
      croots(1,2) = 0.0_dp
      call roots(a_def1, 3, croots, 1)
      root2 = croots(1,1)

      ! Deflate and find third root
      call deflate_polynomial(a_def1, 3, root2, a_def2)
      croots(1,1) = 3.5_dp
      croots(1,2) = 0.0_dp
      call roots(a_def2, 2, croots, 1)
      root3 = croots(1,1)

      ! Verify roots by substituting into original polynomial: x^3 - 6x^2 + 11x - 6
      ! P(x) = -6 + 11x - 6x^2 + x^3
      test_pass = abs(root1**3 - 6.0_dp*root1**2 + 11.0_dp*root1 - 6.0_dp) < 1.0e-6_dp .and. &
                  abs(root2**3 - 6.0_dp*root2**2 + 11.0_dp*root2 - 6.0_dp) < 1.0e-6_dp .and. &
                  abs(root3**3 - 6.0_dp*root3**2 + 11.0_dp*root3 - 6.0_dp) < 1.0e-6_dp

      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Cubic (x-1)(x-2)(x-3) = 0'
         print '(A,F12.6)', '    Root 1: ', root1
         print '(A,F12.6)', '    Root 2: ', root2
         print '(A,F12.6)', '    Root 3: ', root3
      else
         n_failed = n_failed + 1
         print '(A)', '  ✗ Cubic roots incorrect'
         print '(A,F12.6,A,F12.6,A,F12.6)', '    Got: ', root1, ', ', root2, ', ', root3
      end if

      print *, ''
   end subroutine test_roots_cubic

   !---------------------------------------------------------------------------
   !>
   !   Test deflate_polynomial with simple case
   subroutine test_deflate_polynomial_simple(n_tests, n_passed, n_failed)
      integer, intent(inout) :: n_tests, n_passed, n_failed
      real(dp) :: c(3), c_new(3), root
      logical :: test_pass

      print *, 'Testing deflate_polynomial with simple case...'

      ! Test: (x - 2)(x - 3) = x^2 - 5x + 6
      ! Divide by (x - 2), should get x - 3
      ! c = [6, -5, 1], root = 2, result should be [-3, 1]
      n_tests = n_tests + 1
      c(1) = 6.0_dp
      c(2) = -5.0_dp
      c(3) = 1.0_dp
      root = 2.0_dp

      call deflate_polynomial(c, 3, root, c_new)

      ! After deflation, c_new should be [-3, 1, ...]
      test_pass = abs(c_new(1) - (-3.0_dp)) < tol .and. &
                  abs(c_new(2) - 1.0_dp) < tol

      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Deflate (x² - 5x + 6) by (x - 2)'
         print '(A,F12.6,A,F12.6)', '    Result: ', c_new(2), 'x + ', c_new(1)
      else
         n_failed = n_failed + 1
         print '(A)', '  ✗ Deflation incorrect'
         print '(A,F12.6,A,F12.6)', '    Got: ', c_new(2), 'x + ', c_new(1)
      end if

      print *, ''
   end subroutine test_deflate_polynomial_simple

   !---------------------------------------------------------------------------
   !>
   !   Test deflate_polynomial with quadratic
   subroutine test_deflate_polynomial_quadratic(n_tests, n_passed, n_failed)
      integer, intent(inout) :: n_tests, n_passed, n_failed
      real(dp) :: c(4), c_new(4), root
      logical :: test_pass

      print *, 'Testing deflate_polynomial with cubic...'

      ! Test: (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6
      ! Divide by (x - 1), should get x^2 - 5x + 6
      ! c = [-6, 11, -6, 1], root = 1, result should be [6, -5, 1]
      n_tests = n_tests + 1
      c(1) = -6.0_dp
      c(2) = 11.0_dp
      c(3) = -6.0_dp
      c(4) = 1.0_dp
      root = 1.0_dp

      call deflate_polynomial(c, 4, root, c_new)

      ! After deflation, c_new should be [6, -5, 1, ...]
      test_pass = abs(c_new(1) - 6.0_dp) < tol .and. &
                  abs(c_new(2) - (-5.0_dp)) < tol .and. &
                  abs(c_new(3) - 1.0_dp) < tol

      if (test_pass) then
         n_passed = n_passed + 1
         print '(A)', '  ✓ Deflate cubic by (x - 1)'
         print '(A,F12.6,A,F12.6,A,F12.6)', '    Result: ', c_new(3), 'x² + ', c_new(2), 'x + ', c_new(1)
      else
         n_failed = n_failed + 1
         print '(A)', '  ✗ Deflation incorrect'
         print '(A,F12.6,A,F12.6,A,F12.6)', '    Got: ', c_new(3), 'x² + ', c_new(2), 'x + ', c_new(1)
      end if

      print *, ''
   end subroutine test_deflate_polynomial_quadratic

end program test_jacchia_roberts_utilities
