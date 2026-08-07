!------------------------------------------------------------------------------
!>
!   Math utilities for Jacchia-Roberts atmospheric density model.

module jacchia_roberts_utilities
   use jacchia_roberts_kinds, only: ip, dp
   implicit none
   private

   public :: roots
   public :: deflate_polynomial
   public :: date_to_mjd

   contains

   !---------------------------------------------------------------------------
   !>
   !   Finds the roots of a polynomial using Newton's method
   !
   !@note This routine is no longer used.

   pure subroutine roots(a, na, croots, irl)
      real(dp), intent(in) :: a(:) !! Array of polynomial coefficients (lowest power first)
      integer(ip), intent(in) :: na !! Number of coefficients
      real(dp), intent(inout) :: croots(:,:) !! Initial guesses and output roots (real, imaginary)
      integer(ip), intent(in) :: irl !! Number of roots to solve for

      integer(ip) :: i, ir, n1, n2, j, iter
      real(dp) :: z(2), zs(2), cb(2), cc(2), cb_old(2), dif, denom, temp
      real(dp), parameter :: conv_tol = 1.0e-14_dp
      real(dp), parameter :: abs_tol = 1.0e-15_dp
      integer(ip), parameter :: max_iter = 100

      ir = 0
      n1 = na - 1
      n2 = n1 - 1

      do while (ir < irl)
         z(1) = croots(ir+1,1)
         z(2) = croots(ir+1,2)

         iter = 0
         do
            iter = iter + 1
            if (iter > max_iter) then
               ! Exceeded maximum iterations - exit with current approximation
               exit
            end if

            cb(1) = a(n1+1)
            cb(2) = 0.0_dp
            cc(1) = a(n1+1)
            cc(2) = 0.0_dp

            do i = 0, n2
               j = n2 - i
               ! Save old cb values before updating
               cb_old(1) = cb(1)
               cb_old(2) = cb(2)
               ! Update cb (polynomial value)
               temp = (z(1) * cb(1) - z(2) * cb(2)) + a(j+1)
               cb(2) = z(1) * cb(2) + z(2) * cb(1)
               cb(1) = temp
               ! Update cc (derivative) using OLD cb values
               if (j /= 0) then
                  temp = (z(1) * cc(1) - z(2) * cc(2)) + cb_old(1)
                  cc(2) = (z(1) * cc(2) + z(2) * cc(1)) + cb_old(2)
                  cc(1) = temp
               end if
            end do

            zs(1) = z(1)
            zs(2) = z(2)

            ! Newton's method
            denom = cc(1) * cc(1) + cc(2) * cc(2)
            if (denom < 1.0e-30_dp) then
               ! Derivative is zero - cannot continue Newton's method
               exit
            end if
            z(1) = z(1) - ((cb(1) * cc(1) + cb(2) * cc(2)) / denom)
            z(2) = z(2) + ((cb(1) * cc(2) - cb(2) * cc(1)) / denom)

            ! Convergence criterion with protection against division by zero
            if (abs(zs(1)) > abs_tol) then
               dif = abs((zs(1) - z(1)) / zs(1))
            else
               ! Use absolute difference when zs(1) is near zero
               dif = abs(zs(1) - z(1))
            end if

            if (abs(zs(2)) > abs_tol) then
               dif = dif + abs((zs(2) - z(2)) / zs(2))
            else
               ! Use absolute difference when zs(2) is near zero
               dif = dif + abs(zs(2) - z(2))
            end if

            if (dif <= conv_tol) exit
         end do

         croots(ir+1,1) = z(1)
         croots(ir+1,2) = z(2)
         ir = ir + 1
      end do

   end subroutine roots

   !---------------------------------------------------------------------------
   !>
   !   Reduces the order of a polynomial by division
   !
   !@note This routine is no longer used.

   pure subroutine deflate_polynomial(c, n, root, c_new)
      real(dp), intent(in) :: c(:) !! Polynomial coefficients
      integer(ip), intent(in) :: n !! Order + 1 of polynomial
      real(dp), intent(in) :: root !! A single real root of the polynomial
      real(dp), intent(inout) :: c_new(:) !! Output array with new coefficients

      integer(ip) :: i
      real(dp) :: sum, s

      sum = c(n)

      do i = n - 1, 1, -1
         s = c(i)
         c_new(i) = sum
         sum = s + sum * root
      end do

   end subroutine deflate_polynomial

   !---------------------------------------------------------------------------
   !>
   !   Convert calendar date to Modified Julian Date (MJD)

   pure subroutine date_to_mjd(year, month, day, mjd)
      integer(ip), intent(in) :: year, month, day
      real(dp), intent(out) :: mjd

      ! Convert JDN to MJD at midnight (00:00 UTC)
      ! JDN is at noon, so JDN - 0.5 gives midnight
      ! MJD = JD - 2400000.5 = (JDN - 0.5) - 2400000.5 = JDN - 2400001.0
      mjd = real(julian_day(year, month, day) - 2400001, dp)

   end subroutine date_to_mjd

!*****************************************************************************************
!> author: Jacob Williams
!
!  Returns the Julian day number (i.e., the Julian date at Greenwich noon)
!  on the specified YEAR, MONTH, and DAY.
!
!  Valid for any Gregorian calendar date producing a
!  Julian date greater than zero.
!
!### Reference
!   * [USNO](https://aa.usno.navy.mil/faq/JD_formula)

   pure integer function julian_day(y,m,d)

    implicit none

    integer(ip),intent(in) :: y   !! year (YYYY)
    integer(ip),intent(in) :: m   !! month (MM)
    integer(ip),intent(in) :: d   !! day (DD)

    julian_day = d-32075+1461*(y+4800+(m-14)/12)/4+367*&
                 (m-2-(m-14)/12*12)/12-3*((y+4900+(m-14)/12)/100)/4

   end function julian_day
!*****************************************************************************************

end module jacchia_roberts_utilities
