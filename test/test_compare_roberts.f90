!>
!  Comparison of old INPE/Kuga/Carrara Jacchia-Roberts implementation
!  vs. the new Fortran port, at a range of altitudes.
!
!  Both models use identical space weather data (CSSI file) via
!  `prepare_flux_data`, so the output difference reflects only the
!  algorithmic differences between the two implementations:
!
!    90-100 km : old uses 5-point Newton-Cotes quadrature;
!                new uses Roberts analytic closed-form.
!   100-125 km : same Roberts partial-fraction formula in both.
!    >125 km   : same Roberts analytic formula in both.
!
!  Geomagnetic temperature correction:
!    Old (Kuga 1971): smooth tanh blend centred at 350 km.
!    New (port)     : step function at 200 km.
!  Differences above 200 km are dominated by this choice.
!
!  Semi-annual, seasonal-latitudinal, and helium corrections are
!  numerically identical in both implementations.
!
!  Build with fpm:  fpm test test_compare_roberts

program test_compare_roberts

   use inpe_roberts_module,    only: soflud_init, rdymos_cssi, rsdamo
   use jacchia_roberts_module, only: jacchia_roberts_type
   use jacchia_roberts_kinds,  only: ip, dp

   implicit none

   ! Physical constants
   real(dp), parameter :: PI           = acos(-1.0_dp)
   real(dp), parameter :: EARTH_RADIUS = 6356.766_dp    !! Earth polar radius (km)
   real(dp), parameter :: MJD_1950     = 33282.0_dp     !! Standard MJD - MJD-1950 offset
   real(dp), parameter :: km2m         = 1000.0_dp      !! km to m conversion factor

   ! Space weather file (relative to fpm run directory)
   character(len=*), parameter :: SW_FILE = 'data/SpaceWeather-All-v1.2.txt'

   ! Test epoch: same as rober_test default
   !   Rjud = 22476  (MJD-1950)  =>  std MJD = 55758  (~2011-08-08)
   real(dp),  parameter :: Rjud = 22476.0_dp !! MJD-1950
   real(dp),  parameter :: Dafr = 43200.0_dp !! Noon UT (seconds) — avoids midnight singularity

   ! Altitudes (km) to evaluate
   integer, parameter :: N_ALT = 15
   real(dp), parameter :: ALTITUDES(N_ALT) = [ &
      90.0_dp, 100.0_dp, 110.0_dp, 150.0_dp, 200.0_dp, 250.0_dp, &
      300.0_dp, 400.0_dp, 500.0_dp, 600.0_dp, &
      800.0_dp, 1000.0_dp, 1200.0_dp, 1500.0_dp, 2000.0_dp ]

   ! ---- Old-model I/O ----
   real(dp)  :: Sa(3), Su(2), Gsti
   real(dp)  :: Te_old(2), Ad_old(6), Wmol_old, Rhod_old
   integer  :: old_status

   ! ---- New-model ----
   type(jacchia_roberts_type) :: jr
   real(dp) :: position(3), sun_vector(3), geo_lat, utc_mjd
   real(dp) :: density_new
   integer(ip) :: sw_init_status

   ! ---- Loop ----
   integer  :: i
   real(dp) :: alt_km, rhod_old_dp, log_old, log_new, pct_diff
   character(len=6) :: flag

   ! ---- Also collect results for the static comparison using rsdamo directly ----
   !      (hardcoded SF to isolate algorithm from data source)
   real(dp), parameter :: SF_FIXED(3) = [150.0_dp, 150.0_dp, 3.0_dp]  !! F10.7, F10.7a, Kp
   real(dp) :: Te_stat(2), Ad_stat(6), Wmol_stat, Rhod_stat

   !--------------------------------------------------------------------------
   ! Initialize both models
   !--------------------------------------------------------------------------

   utc_mjd = real(Rjud, dp) + MJD_1950 + real(Dafr, dp) / 86400.0_dp

   call soflud_init(SW_FILE, old_status)
   if (old_status /= 0) then
      error stop 'ERROR: could not initialize old-model space weather data.'
   end if

   call jr%initialize(EARTH_RADIUS, filename=SW_FILE, status=sw_init_status)
   if (sw_init_status /= 0) then
      write(*,'(A)') 'WARNING: new model space weather load failed; using nominal values.'
   end if

   !--------------------------------------------------------------------------
   ! Fixed geometry
   !   Equatorial position, sun at RA = 0, noon UT.
   !   Both models receive the same pointing: satellite on the x-axis (RA=0,
   !   dec=0); sun also at RA=0, dec=0 (overhead at local noon, equator).
   !--------------------------------------------------------------------------

   Su(1) = 0.0_dp  ! sun RA  = 0
   Su(2) = 0.0_dp  ! sun dec = 0
   Sa(1) = 0.0_dp  ! satellite RA  = 0
   Sa(2) = 0.0_dp  ! satellite geocentric lat = 0 (equator)
   Gsti  = 0.0_dp  ! Greenwich sidereal time (not used by rsdamo)

   sun_vector = [1.0_dp, 0.0_dp, 0.0_dp]   ! unit vector to sun (RA=0)
   geo_lat    = 0.0_dp                       ! geodetic latitude (degrees)

   !--------------------------------------------------------------------------
   ! Table 1: Live space weather — rdymos_cssi (old algo) vs jr%density (new)
   !--------------------------------------------------------------------------

   write(*,'(A)') ''
   write(*,'(A)') '======================================================================='
   write(*,'(A)') '  Table 1: Dynamic comparison — same CSSI space weather data'
   write(*,'(A)') '           Old: INPE/Kuga/Carrara dyjrmo    New: Fortran port'
   write(*,'(A)') '======================================================================='
   write(*,'(2A,F10.1)') '  Epoch (MJD-1950): ', '   Rjud = ', Rjud
   write(*,'(2A,F10.1)') '  Epoch (std MJD) : ', ' utc_mjd = ', utc_mjd
   write(*,'(2A,F8.1)') '  Time of day (s)  : ', '   Dafr = ', Dafr
   write(*,'(A)') '  Geometry: equatorial, sat and sun both at RA=0'
   write(*,'(A)') ''
   write(*,'(A)') '  Alt(km)  Tinf_old(K)  log10(rho_old)  log10(rho_new)  Delta_rho(%)'
   write(*,'(A)') '  -------  -----------  --------------  --------------  ------------'

   do i = 1, N_ALT
      alt_km   = ALTITUDES(i)
      Sa(3)    = alt_km * km2m      ! metres for old model
      position = [(EARTH_RADIUS + alt_km), 0.0_dp, 0.0_dp]   ! km for new model

      ! Old model via rdymos_cssi (uses same prepare_flux_data timing)
      call rdymos_cssi(Sa, Su, Rjud, Dafr, Gsti, Te_old, Ad_old, Wmol_old, Rhod_old, old_status)
      if (old_status /= 0) then
         write(*,'(2X,F7.1,A)') alt_km, '  [old model lookup failed - skipping]'
         cycle
      end if

      ! New model
      density_new = jr%density(alt_km, position, sun_vector, geo_lat, utc_mjd)

      rhod_old_dp = real(Rhod_old, dp)
      log_old  = log10(rhod_old_dp)
      log_new  = log10(density_new)
      pct_diff = 100.0_dp * (density_new - rhod_old_dp) / rhod_old_dp

      ! Flag rows where the two models disagree by more than 2 log-units
      ! (indicates numerical failure, not an algorithmic difference)
      if (abs(log_new - log_old) > 2.0_dp) then
         flag = '  <BUG'
      else
         flag = ''
      end if

      write(*,'(2X,F7.1,2X,F11.2,2X,F14.6,2X,F14.6,2X,F12.4,A)') &
         alt_km, real(Te_old(1), dp), log_old, log_new, pct_diff, trim(flag)
   end do

   write(*,'(A)') ''
   write(*,'(A)') '  Delta_rho = 100 * (rho_new - rho_old) / rho_old'
   write(*,'(A)') '  <BUG = |log10(new) - log10(old)| > 2 (not an algorithm difference)'
   write(*,'(A)') ''
   write(*,'(A)') '  Algorithmic sources of difference:'
   write(*,'(A)') '    90-200 km : Kp correction step at 200 km (new) vs tanh@350 km (old)'
   write(*,'(A)') '    90-100 km : Roberts analytic 100-km formula (new) vs'
   write(*,'(A)') '                5-point Newton-Cotes quadrature (old)'
   write(*,'(A)') '   100-125 km : identical Roberts partial-fraction formula'
   write(*,'(A)') '    >125 km   : identical Roberts analytic formula; Kp-correction method differs'

   !--------------------------------------------------------------------------
   ! Table 2: Algorithm isolation — rsdamo with hardcoded SF vs new density
   !          Removes any space-weather timing difference; isolates pure
   !          algorithmic difference (tanh@350 vs step@200, quadrature vs analytic)
   !--------------------------------------------------------------------------

   write(*,'(A)') ''
   write(*,'(A)') '======================================================================='
   write(*,'(A)') '  Table 2: Algorithm isolation — hardcoded SF = [150, 150, 3]'
   write(*,'(A)') '           Old: rsdamo direct    New: jr%density (reads from file)'
   write(*,'(A)') '  Note: new model still reads file; this table shows rsdamo baseline'
   write(*,'(A)') '======================================================================='
   write(*,'(A)') ''
   write(*,'(A)') '  Alt(km)  Tinf_old(K)  log10(rho_old)   Wmol_old  rho_old(kg/m3)'
   write(*,'(A)') '  -------  -----------  --------------   --------  ---------------'

   do i = 1, N_ALT
      alt_km   = ALTITUDES(i)
      Sa(3)    = alt_km * km2m

      ! rsdamo requires Rjud in MJD-1950; rsdamo internally does amjd = Rjud+33282+Dafr/86400
      call rsdamo(Sa, Su, SF_FIXED, Rjud, Dafr, Gsti, Te_stat, Ad_stat, Wmol_stat, Rhod_stat)

      write(*,'(2X,F7.1,2X,F11.2,2X,F14.6,2X,F10.4,2X,ES15.6)') &
         alt_km, real(Te_stat(1), dp), log10(real(Rhod_stat, dp)), &
         real(Wmol_stat, dp), real(Rhod_stat, dp)
   end do

   write(*,'(A)') ''
   write(*,'(A)') '  Species log10 number densities at 300 km (F10.7=150, F10.7a=150, Kp=3):'
   write(*,'(A)') '  (He, O2, N2, Ar, O, H) from old rsdamo'
   write(*,'(A)') ''

   ! Print species densities at 300 km for inspection
   alt_km = 300.0_dp
   Sa(3)  = alt_km * km2m
   call rsdamo(Sa, Su, SF_FIXED, Rjud, Dafr, Gsti, Te_stat, Ad_stat, Wmol_stat, Rhod_stat)
   write(*,'(2X,A,F7.2,A)') 'At alt = ', alt_km, ' km:'
   write(*,'(4X,A,6F9.4)') 'log10 n: ', real(Ad_stat(1:6), dp)
   write(*,'(4X,A,6(3X,A))') 'Species: ', '    He', '    O2', '    N2', '    Ar', '     O', '     H'

   ! Clean up
   call jr%destroy()

end program test_compare_roberts
