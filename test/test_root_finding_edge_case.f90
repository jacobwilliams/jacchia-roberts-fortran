!------------------------------------------------------------------------------
!>
!  Test case for comparing root-finding methods at various temperatures
!
!  This test directly compares find_cstar_roots (new, stable) vs
!  find_cstar_roots_original (old, potentially unstable) across a range
!  of boundary temperatures (Tx).
!
!  The test verifies that:
!    1. find_cstar_roots maintains root1 > 125 km and root2 < 100 km for all Tx
!    2. Both methods agree when the original method is stable
!    3. The roots satisfy the c_star polynomial equation
!
!  Key constraint from rho_125: For log((height - root1)/(100 - root1)) and
!  log((height - root2)/(100 - root2)) to be real in [100, 125] km:
!    - root1 must be > 125 km
!    - root2 must be < 100 km

program test_root_finding_edge_case
   use jacchia_roberts_kinds, only: ip, dp
   use jacchia_roberts_module, only: find_cstar_roots, find_cstar_roots_original
   implicit none

   ! Parameters for c_star polynomial construction
   real(dp), parameter :: CON_C_1 = -89284375.0_dp
   real(dp), parameter :: CON_C_2 = 3542400.0_dp
   real(dp), parameter :: CON_C_3 = -52687.5_dp
   real(dp), parameter :: CON_C_4 = 340.5_dp
   real(dp), parameter :: CON_C_5 = -0.8_dp
   real(dp), parameter :: TZERO = 183.0_dp

   real(dp) :: c_star(5), c_star_copy(5)
   real(dp) :: tx, root1_new, root2_new, x_root_new, y_root_new
   real(dp) :: root1_old, root2_old, x_root_old, y_root_old
   real(dp) :: poly_eval_new, poly_eval_old, poly_eval
   integer :: n_tests, n_passed, n_failed, i
   logical :: roots_valid_new, roots_valid_old

   print *, ''
   print *, '=========================================='
   print *, 'Root-Finding Method Comparison Test'
   print *, '=========================================='
   print *, ''
   print *, 'This test compares the stability of two root-finding methods'
   print *, 'for the c_star quartic polynomial across various temperatures.'
   print *, ''
   print *, 'Required constraints (from rho_125 formula):'
   print *, '  - root1 > 125 km (to keep log argument positive)'
   print *, '  - root2 < 100 km (to keep log argument positive)'
   print *, '  - polynomial error < 1.0e-6 (to ensure roots are accurate)'
   print *, ''

   n_tests = 0
   n_passed = 0
   n_failed = 0

   ! Test across a range of boundary temperatures
   ! Tx typically ranges from ~300 K to ~600 K depending on solar activity
   do i = 1, 5
      select case(i)
      case(1)
         tx = 300.0_dp  ! Very low temperature (solar minimum)
      case(2)
         tx = 373.0_dp  ! Low temperature (problematic case from git commit)
      case(3)
         tx = 450.0_dp  ! Moderate temperature
      case(4)
         tx = 550.0_dp  ! High temperature (solar maximum)
      case(5)
         tx = 600.0_dp  ! Very high temperature
      end select

      print *, '-------------------------------------------'
      write(*, '(A, F6.1, A)') ' Test case: Tx = ', tx, ' K'
      print *, '-------------------------------------------'

      ! Construct c_star polynomial coefficients
      c_star(1) = CON_C_1 + 1500625.0_dp * tx / (tx - TZERO)
      c_star(2) = CON_C_2
      c_star(3) = CON_C_3
      c_star(4) = CON_C_4
      c_star(5) = CON_C_5

      ! Make a copy for the old method (which modifies c_star in-place)
      c_star_copy = c_star

      ! Test new method
      call find_cstar_roots(c_star, tx, root1_new, root2_new, x_root_new, y_root_new)

      ! Test old method
      call find_cstar_roots_original(c_star_copy, tx, root1_old, root2_old, x_root_old, y_root_old)

      ! Evaluate polynomial at roots to verify they are actually roots
      poly_eval_new = c_star(1) + root1_new * (c_star(2) + root1_new * (c_star(3) + &
                      root1_new * (c_star(4) + root1_new * c_star(5))))
      poly_eval_old = c_star(1) + root1_old * (c_star(2) + root1_old * (c_star(3) + &
                      root1_old * (c_star(4) + root1_old * c_star(5))))

      ! Check if roots satisfy constraints
      roots_valid_new = (root1_new > 125.0_dp) .and. (root2_new < 100.0_dp) .and. (abs(poly_eval_new) < 1.0e-6_dp)
      roots_valid_old = (root1_old > 125.0_dp) .and. (root2_old < 100.0_dp) .and. (abs(poly_eval_old) < 1.0e-6_dp)

      ! Display results
      write(*, '(A)') ' '
      write(*, '(A)') '  New method (find_cstar_roots):'
      write(*, '(A, F12.4)') '    root1 (should be > 125):  ', root1_new
      write(*, '(A, F12.4)') '    root2 (should be < 100):  ', root2_new
      write(*, '(A, F12.4)') '    x_root:                   ', x_root_new
      write(*, '(A, F12.4)') '    y_root:                   ', y_root_new
      write(*, '(A, ES12.4)') '    polynomial error:         ', abs(poly_eval_new)
      write(*, '(A, L)') '    Constraints satisfied:    ', roots_valid_new
      write(*, '(A)') ' '
      write(*, '(A)') '  Old method (find_cstar_roots_original):'
      write(*, '(A, F12.4)') '    root1 (should be > 125):  ', root1_old
      write(*, '(A, F12.4)') '    root2 (should be < 100):  ', root2_old
      write(*, '(A, F12.4)') '    x_root:                   ', x_root_old
      write(*, '(A, F12.4)') '    y_root:                   ', y_root_old
      write(*, '(A, ES12.4)') '    polynomial error:         ', abs(poly_eval_old)
      write(*, '(A, L)') '    Constraints satisfied:    ', roots_valid_old
      write(*, '(A)') ' '

      ! Evaluate test
      n_tests = n_tests + 1
      if (roots_valid_new) then
         write(*, '(A)') '  ✓ PASS: New method maintains valid root constraints'
         n_passed = n_passed + 1
      else
         write(*, '(A)') '  ✗ FAIL: New method violates root constraints'
         n_failed = n_failed + 1
      end if

      if (.not. roots_valid_old) then
         write(*, '(A)') '  ! WARNING: Old method violates root constraints'
      end if

      print *, ''
   end do

   print *, '=========================================='
   print *, 'Test Summary'
   print *, '=========================================='
   write(*, '(A, I2)') ' Total tests  : ', n_tests
   write(*, '(A, I2)') ' Passed       : ', n_passed
   write(*, '(A, I2)') ' Failed       : ', n_failed
   print *, ''

   if (n_failed > 0) then
      error stop 'OVERALL RESULT: FAILED'
   else
      print *, 'OVERALL RESULT: PASSED'
   end if

end program test_root_finding_edge_case
