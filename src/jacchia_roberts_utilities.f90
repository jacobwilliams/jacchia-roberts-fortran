!------------------------------------------------------------------------------
!>
!   Math utilities for Jacchia-Roberts atmospheric density model.

module jacchia_roberts_utilities
   use, intrinsic :: iso_fortran_env, only: dp => real64, ip => int32
   implicit none
   private

   public :: roots
   public :: deflate_polynomial
   public :: date_to_mjd

   contains

   !---------------------------------------------------------------------------
   !>
   !   Finds the roots of a polynomial using Newton's method

   subroutine roots(a, na, croots, irl)
      real(dp), intent(in) :: a(:) !! Array of polynomial coefficients (lowest power first)
      integer(ip), intent(in) :: na !! Number of coefficients
      real(dp), intent(inout) :: croots(:,:) !! Initial guesses and output roots (real, imaginary)
      integer(ip), intent(in) :: irl !! Number of roots to solve for

      integer(ip) :: i, ir, n1, n2, j
      real(dp) :: z(2), zs(2), cb(2), cc(2), dif, denom, temp
      real(dp), parameter :: conv_tol = 1.0e-14_dp

      ir = 0
      n1 = na - 1
      n2 = n1 - 1

      do while (ir < irl)
         z(1) = croots(ir+1,1)
         z(2) = croots(ir+1,2)

         do
            cb(1) = a(n1+1)
            cb(2) = 0.0_dp
            cc(1) = a(n1+1)
            cc(2) = 0.0_dp

            do i = 0, n2
               j = n2 - i
               temp = (z(1) * cb(1) - z(2) * cb(2)) + a(j+1)
               cb(2) = z(1) * cb(2) + z(2) * cb(1)
               cb(1) = temp
               if (j /= 0) then
                  temp = (z(1) * cc(1) - z(2) * cc(2)) + cb(1)
                  cc(2) = (z(1) * cc(2) + z(2) * cc(1)) + cb(2)
                  cc(1) = temp
               end if
            end do

            zs(1) = z(1)
            zs(2) = z(2)

            ! Newton's method
            denom = cc(1) * cc(1) + cc(2) * cc(2)
            z(1) = z(1) - ((cb(1) * cc(1) + cb(2) * cc(2)) / denom)
            z(2) = z(2) + ((cb(1) * cc(2) - cb(2) * cc(1)) / denom)

            ! Convergence criterion
            dif = abs((zs(1) - z(1)) / zs(1))
            if (abs(zs(2)) > 1.0e-15_dp) then
               dif = dif + abs((zs(2) - z(2)) / zs(2))
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

   subroutine deflate_polynomial(c, n, root, c_new)
      real(dp), intent(in) :: c(:) !! Polynomial coefficients
      integer(ip), intent(in) :: n !! Order + 1 of polynomial
      real(dp), intent(in) :: root !! A single real root of the polynomial
      real(dp), intent(inout) :: c_new(:) !! Output array with new coefficients

      integer(ip) :: i
      real(dp) :: sum, save

      sum = c(n)

      do i = n - 1, 1, -1
         save = c(i)
         c_new(i) = sum
         sum = save + sum * root
      end do

   end subroutine deflate_polynomial

   !---------------------------------------------------------------------------
   !>
   !   Convert calendar date to Modified Julian Date (MJD)
   pure subroutine date_to_mjd(year, month, day, mjd)
      integer(ip), intent(in) :: year, month, day
      real(dp), intent(out) :: mjd
      integer(ip) :: a, y, m, jdn

      ! Convert to Julian Day Number using standard algorithm
      a = (14 - month) / 12
      y = year + 4800 - a
      m = month + 12 * a - 3

      jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045

      ! Convert JDN to MJD
      mjd = real(jdn, dp) - 2400000.5_dp

   end subroutine date_to_mjd

end module jacchia_roberts_utilities