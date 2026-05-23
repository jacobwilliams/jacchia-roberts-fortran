!>
!  Comparison of old INPE/Kuga/Carrara Jacchia-Roberts implementation
!  vs. the new Fortran port, over a range of inputs.
!
!  Table 1 : Full altitude sweep, reference epoch and geometry.
!  Table 2 : rsdamo baseline with fixed space weather (algorithm isolation).
!  Table 3 : Epoch sensitivity   — solar min (~1996), moderate (~2011), solar max (~2014).
!  Table 4 : Sun-declination sensitivity — winter solstice, equinox, summer solstice.
!  Table 5 : Geographic-latitude sensitivity — equatorial, mid-lat, sub-polar.
!  Table 6 : Local-time sensitivity — noon, evening, midnight, dawn.
!
!  Δ% = 100 × (ρ_new − ρ_old) / ρ_old   (new = Fortran port, old = INPE)
!
!  Known algorithmic differences:
!    90-100 km : new uses original GTDS/Roberts analytical formula; old uses
!                5-pt Newton-Cotes numerical integration. Both methods produce
!                proper exponential decay and agree to within 0.4%.
!                If using Vallado's formula, it gave 143-244%
!                errors. Roberts (1971) claimed "identical" results.
!   100-125 km : identical Roberts partial-fraction formula.
!    >125 km   : same formula; Kp correction differs:
!                new: step at 200 km    old: tanh blend centred at 350 km.
!    ~300 km   : largest difference (~1.3%) from the Kp-correction boundary.
!
!  Geometry mapping old -> new:
!    Sa(2) [geocentric lat, rad]  ≈ geo_lat [geodetic deg]  (small difference)
!    Su(2) [sun dec, rad]         -> sun_vector via [cos(d),0,sin(d)]
!    h = Sa(1)−Su(1) [hour angle] -> satellite xy-position angle
!
!  Build:  fpm test test_compare_roberts

program test_compare_roberts_new

   use inpe_roberts_module,    only: soflud_init, rdymos_cssi, rsdamo
   use jacchia_roberts_module, only: jacchia_roberts_type
   use jacchia_roberts_kinds,  only: ip, dp

   implicit none

   real(dp), parameter :: PI           = acos(-1.0_dp)
   real(dp), parameter :: EARTH_RADIUS = 6356.766_dp    !! Earth polar radius (km)
   real(dp), parameter :: MJD_1950     = 33282.0_dp     !! Standard MJD − MJD-1950 offset
   character(len=*), parameter :: SW_FILE = 'data/SpaceWeather-All-v1.2.txt'
   real(dp), parameter :: km2m         = 1000.0_dp      !! km to m conversion factor

   !--- Reference epoch for Tables 1 & 2 ---
   real(dp),  parameter :: Rjud_ref = 22476.0_dp  !! MJD-1950 (~2011-08-08)
   real(dp),  parameter :: Dafr     = 43200.0_dp  !! Noon UT (seconds)

   !--- Full altitude array for Table 1 ---
   integer,  parameter :: N_ALT = 17
   real(dp), parameter :: ALTITUDES(N_ALT) = [ &
      90.0_dp, 95.0_dp, 97.0_dp, 100.0_dp, 110.0_dp, 150.0_dp, 200.0_dp, 250.0_dp, &
      300.0_dp, 400.0_dp, 500.0_dp, 600.0_dp, &
      800.0_dp, 1000.0_dp, 1200.0_dp, 1500.0_dp, 2000.0_dp ]

   !--- Key altitudes for compact sensitivity tables (Tables 3-6) ---
   integer,  parameter :: N_KEY = 5
   real(dp), parameter :: KEY_ALTS(N_KEY) = &
      [90.0_dp, 110.0_dp, 200.0_dp, 300.0_dp, 1000.0_dp]

   !--- Sensitivity parameters ---

   ! Table 3: Epochs (MJD-1950) — solar min, moderate, solar max
   integer,  parameter :: N_EPOCHS = 3
   real(dp),  parameter :: TEST_RJUD(N_EPOCHS) = &
      [16800.0_dp, 22476.0_dp, 23618.0_dp]
   character(len=22), parameter :: EPOCH_DESC(N_EPOCHS) = [ character(len=22) :: &
      '~1996 solar min       ', &
      '~2011 moderate        ', &
      '~2014 solar max       ']

   ! Table 4: Sun declinations (degrees)
   integer,  parameter :: N_SUNS = 3
   real(dp), parameter :: SUN_DEC_DEG(N_SUNS) = [-23.45_dp, 0.0_dp, 23.45_dp]
   character(len=22), parameter :: SUN_DESC(N_SUNS) = [ character(len=22) :: &
      'winter solst. -23.45  ', &
      'equinox         0.00  ', &
      'summer solst. +23.45  ']

   ! Table 5: Geographic latitudes (degrees)
   integer,  parameter :: N_LATS = 3
   real(dp), parameter :: GEO_LAT_DEG(N_LATS) = [0.0_dp, 45.0_dp, 75.0_dp]
   character(len=22), parameter :: LAT_DESC(N_LATS) = [ character(len=22) :: &
      'equatorial   lat=  0N ', &
      'mid-latitude lat= 45N ', &
      'sub-polar    lat= 75N ']

   ! Table 6: Local times — satellite RA relative to sun at RA=0 (hour angle)
   integer,  parameter :: N_TIMES = 4
   real(dp),  parameter :: SAT_RA_RAD(N_TIMES) = &
      [0.0_dp, 1.5707963268_dp, 3.1415926536_dp, -1.5707963268_dp]
   character(len=22), parameter :: TIME_DESC(N_TIMES) = [ character(len=22) :: &
      'noon      h=   0 deg  ', &
      'evening   h=  +90 deg ', &
      'midnight  h= 180 deg  ', &
      'dawn      h=  -90 deg ']

   !--- Old model I/O ---
   real(dp)  :: Sa(3), Su(2), Gsti
   real(dp)  :: Te_old(2), Ad_old(6), Wmol_old, Rhod_old
   integer  :: old_status

   !--- New model ---
   type(jacchia_roberts_type) :: jr
   real(dp) :: position(3), sun_vector(3), geo_lat, utc_mjd
   real(dp) :: density_new
   integer(ip) :: sw_init_status

   !--- Loop and working variables ---
   integer  :: i, j, k
   real(dp) :: alt_km, rhod_old_dp, log_old, log_new, pct_diff, abs_error
   character(len=6) :: flag
   real(dp) :: pct_arr(N_KEY)
   real(dp)  :: Rjud_local
   real(dp) :: utc_mjd_local, dec_rad, lat_rad, sat_ra

   !--- Static rsdamo baseline (Table 2) ---
   real(dp),  parameter :: SF_FIXED(3) = [150.0_dp, 150.0_dp, 3.0_dp]
   real(dp)  :: Te_stat(2), Ad_stat(6), Wmol_stat, Rhod_stat

   !--------------------------------------------------------------------------
   ! Initialize both models
   !--------------------------------------------------------------------------

   utc_mjd = real(Rjud_ref, dp) + MJD_1950 + real(Dafr, dp) / 86400.0_dp

   call soflud_init(SW_FILE, old_status)
   if (old_status /= 0) then
      error stop 'ERROR: could not initialize old-model space weather data.'
   end if

   call jr%initialize(EARTH_RADIUS, filename=SW_FILE, status=sw_init_status)
   if (sw_init_status /= 0) then
      write(*,'(A)') 'WARNING: new model space weather load failed; using nominal values.'
   end if

   Gsti = 0.0_dp  ! not used by rsdamo/rdymos_cssi but required by signature

   !--- Reference geometry: equatorial, noon, sun at equinox ---
   Su(1) = 0.0_dp;  Su(2) = 0.0_dp   ! sun RA=0, dec=0
   Sa(1) = 0.0_dp;  Sa(2) = 0.0_dp   ! satellite RA=0, geocentric lat=0
   sun_vector = [1.0_dp, 0.0_dp, 0.0_dp]
   geo_lat    = 0.0_dp

   !==========================================================================
   ! Table 1 : Full altitude sweep — reference epoch & geometry
   !==========================================================================

   write(*,'(A)') ''
   write(*,'(A)') '======================================================================='
   write(*,'(A)') '  Table 1: Full altitude sweep — CSSI space weather, reference geometry'
   write(*,'(A)') '           Old: INPE/Kuga/Carrara dyjrmo    New: Fortran port'
   write(*,'(A)') '======================================================================='
   write(*,'(2A,F10.1)') '  Epoch (MJD-1950):', '   Rjud = ', Rjud_ref
   write(*,'(2A,F10.1)') '  Epoch (std MJD) :', ' utc_mjd = ', utc_mjd
   write(*,'(2A,F8.1)')  '  Time of day (s) :', '    Dafr = ', Dafr
   write(*,'(A)') '  Geometry: equatorial, noon (h=0), sun at equinox (dec=0)'
   write(*,'(A)') ''
   write(*,'(A)') '  Alt(km)  Tinf_old(K)  rho_old         rho_new       Delta_rho(%)  abs_error   '
   write(*,'(A)') '  -------  -----------  --------------  ------------  ------------  ------------'

   do i = 1, N_ALT
      alt_km   = ALTITUDES(i)
      Sa(3)    = alt_km * km2m
      position = [(EARTH_RADIUS + alt_km), 0.0_dp, 0.0_dp]

      call rdymos_cssi(Sa, Su, Rjud_ref, Dafr, Gsti, Te_old, Ad_old, Wmol_old, Rhod_old, old_status)
      if (old_status /= 0) then
         write(*,'(2X,F7.1,A)') alt_km, '  [data lookup failed]'
         cycle
      end if

      density_new = jr%density(alt_km, position, sun_vector, geo_lat, utc_mjd)

      rhod_old_dp = real(Rhod_old, dp)
      log_old  = log10(rhod_old_dp)
      log_new  = log10(density_new)
      pct_diff = 100.0_dp * (density_new - rhod_old_dp) / rhod_old_dp
      abs_error = abs(density_new - rhod_old_dp)

      if (abs(log_new - log_old) > 2.0_dp) then
         flag = '  <BUG'
      else
         flag = ''
      end if

      write(*,'(2X,F7.1,2X,F11.2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6,A)') &
         alt_km, real(Te_old(1), dp), rhod_old_dp, density_new, pct_diff, abs_error, trim(flag)
   end do

   write(*,'(A)') ''
   write(*,'(A)') '  Delta_rho = 100 * (rho_new - rho_old) / rho_old'
   write(*,'(A)') '  <BUG = |log10(new) - log10(old)| > 2 (numerical failure, not algorithm)'

   !========================================
   ! 90 to 100 tests
   !========================================
   write(*,'(A)') ''
   write(*,'(A)') '  Alt(km)  Tinf_old(K)  rho_old         rho_new       Delta_rho(%)  abs_error   '
   write(*,'(A)') '  -------  -----------  --------------  ------------  ------------  ------------'
   do i = 90, 100
      alt_km   = real(i, dp)
      Sa(3)    = alt_km * km2m
      position = [(EARTH_RADIUS + alt_km), 0.0_dp, 0.0_dp]

      call rdymos_cssi(Sa, Su, Rjud_ref, Dafr, Gsti, Te_old, Ad_old, Wmol_old, Rhod_old, old_status)
      if (old_status /= 0) then
         write(*,'(2X,F7.1,A)') alt_km, '  [data lookup failed]'
         cycle
      end if

      density_new = jr%density(alt_km, position, sun_vector, geo_lat, utc_mjd)

      rhod_old_dp = real(Rhod_old, dp)
      log_old  = log10(rhod_old_dp)
      log_new  = log10(density_new)
      pct_diff = 100.0_dp * (density_new - rhod_old_dp) / rhod_old_dp
      abs_error = abs(density_new - rhod_old_dp)

      if (abs(log_new - log_old) > 2.0_dp) then
         flag = '  <BUG'
      else
         flag = ''
      end if

      write(*,'(2X,F7.1,2X,F11.2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6,A)') &
         alt_km, real(Te_old(1), dp), rhod_old_dp, density_new, pct_diff, abs_error, trim(flag)
   end do

   !==========================================================================
   ! Table 2 : rsdamo baseline — fixed SF, algorithm isolation
   !==========================================================================

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
      alt_km = ALTITUDES(i)
      Sa(3)  = alt_km * km2m
      call rsdamo(Sa, Su, SF_FIXED, Rjud_ref, Dafr, Gsti, Te_stat, Ad_stat, Wmol_stat, Rhod_stat)
      write(*,'(2X,F7.1,2X,F11.2,2X,F14.6,2X,F10.4,2X,ES15.6)') &
         alt_km, real(Te_stat(1), dp), log10(real(Rhod_stat, dp)), &
         real(Wmol_stat, dp), real(Rhod_stat, dp)
   end do

   write(*,'(A)') ''
   write(*,'(A)') '  Species log10 number densities at 300 km (F10.7=150, F10.7a=150, Kp=3):'
   write(*,'(A)') '  (He, O2, N2, Ar, O, H) from old rsdamo'
   write(*,'(A)') ''

   alt_km = 300.0_dp
   Sa(3)  = alt_km * km2m
   call rsdamo(Sa, Su, SF_FIXED, Rjud_ref, Dafr, Gsti, Te_stat, Ad_stat, Wmol_stat, Rhod_stat)
   write(*,'(2X,A,F7.2,A)') 'At alt = ', alt_km, ' km:'
   write(*,'(4X,A,6F9.4)') 'log10 n: ', real(Ad_stat(1:6), dp)
   write(*,'(4X,A,6(3X,A))') 'Species: ', '    He', '    O2', '    N2', '    Ar', '     O', '     H'

   !==========================================================================
   ! Sensitivity tables (Tables 3-6)
   !
   ! One parameter is varied per table; the others are held at the reference:
   !   epoch  = Rjud_ref (~2011)
   !   sun    = equinox (dec=0), sun RA=0
   !   lat    = equatorial (0 deg)
   !   h      = noon (satellite RA=0, h=0)
   !==========================================================================

   !==========================================================================
   ! Table 3 : Epoch sensitivity
   !==========================================================================

   write(*,'(A)') ''
   write(*,'(A)') '======================================================================='
   write(*,'(A)') '  Table 3: Epoch sensitivity   Delta% = 100*(rho_new - rho_old)/rho_old'
   write(*,'(A)') '  Fixed geometry: equatorial, noon, sun at equinox'
   write(*,'(A)') '  Rjud 16800 ~ 1996-01 (solar min);  23618 ~ 2014-10 (solar max)'
   write(*,'(A)') '======================================================================='
   write(*,'(A)') ''
   write(*,'(A)') '  Epoch                      90km   110km   200km   300km  1000km'
   write(*,'(A)') '  ----------------------   ------  ------  ------  ------  ------'

   Su(1) = 0.0_dp;  Su(2) = 0.0_dp
   Sa(1) = 0.0_dp;  Sa(2) = 0.0_dp
   sun_vector = [1.0_dp, 0.0_dp, 0.0_dp]
   geo_lat    = 0.0_dp

   do k = 1, N_EPOCHS
      Rjud_local    = TEST_RJUD(k)
      utc_mjd_local = real(Rjud_local, dp) + MJD_1950 + real(Dafr, dp) / 86400.0_dp
      do j = 1, N_KEY
         alt_km   = KEY_ALTS(j)
         Sa(3)    = alt_km * km2m
         position = [(EARTH_RADIUS + alt_km), 0.0_dp, 0.0_dp]
         call rdymos_cssi(Sa, Su, Rjud_local, Dafr, Gsti, &
                          Te_old, Ad_old, Wmol_old, Rhod_old, old_status)
         if (old_status /= 0) then
            pct_arr(j) = -999.0_dp
         else
            density_new = jr%density(alt_km, position, sun_vector, geo_lat, utc_mjd_local)
            pct_arr(j)  = 100.0_dp * (density_new - real(Rhod_old, dp)) / real(Rhod_old, dp)
         end if
      end do
      write(*,'(2X,A22,5F8.2)') EPOCH_DESC(k), (pct_arr(j), j=1,N_KEY)
   end do

   !==========================================================================
   ! Table 4 : Sun-declination sensitivity
   !==========================================================================

   write(*,'(A)') ''
   write(*,'(A)') '======================================================================='
   write(*,'(A)') '  Table 4: Sun-declination sensitivity   Delta% = 100*(rho_new - rho_old)/rho_old'
   write(*,'(A)') '  Fixed: reference epoch (~2011), equatorial, noon'
   write(*,'(A)') '======================================================================='
   write(*,'(A)') ''
   write(*,'(A)') '  Sun declination (deg)      90km   110km   200km   300km  1000km'
   write(*,'(A)') '  ----------------------   ------  ------  ------  ------  ------'

   Sa(1) = 0.0_dp;  Sa(2) = 0.0_dp
   Rjud_local    = Rjud_ref
   utc_mjd_local = utc_mjd
   geo_lat       = 0.0_dp

   do k = 1, N_SUNS
      dec_rad    = SUN_DEC_DEG(k) * (PI / 180.0_dp)
      Su(2)      = real(dec_rad, 8)
      sun_vector = [cos(dec_rad), 0.0_dp, sin(dec_rad)]
      do j = 1, N_KEY
         alt_km   = KEY_ALTS(j)
         Sa(3)    = alt_km * km2m
         position = [(EARTH_RADIUS + alt_km), 0.0_dp, 0.0_dp]
         call rdymos_cssi(Sa, Su, Rjud_local, Dafr, Gsti, &
                          Te_old, Ad_old, Wmol_old, Rhod_old, old_status)
         if (old_status /= 0) then
            pct_arr(j) = -999.0_dp
         else
            density_new = jr%density(alt_km, position, sun_vector, geo_lat, utc_mjd_local)
            pct_arr(j)  = 100.0_dp * (density_new - real(Rhod_old, dp)) / real(Rhod_old, dp)
         end if
      end do
      write(*,'(2X,A22,5F8.2)') SUN_DESC(k), (pct_arr(j), j=1,N_KEY)
   end do
   Su(2) = 0.0_dp;  sun_vector = [1.0_dp, 0.0_dp, 0.0_dp]   ! restore equinox

   !==========================================================================
   ! Table 5 : Geographic-latitude sensitivity
   !==========================================================================

   write(*,'(A)') ''
   write(*,'(A)') '======================================================================='
   write(*,'(A)') '  Table 5: Geographic-latitude sensitivity   Delta% = 100*(rho_new - rho_old)/rho_old'
   write(*,'(A)') '  Fixed: reference epoch (~2011), sun at equinox, noon'
   write(*,'(A)') '  Sa(2) = geocentric lat (rad) ~ geo_lat = geodetic lat (deg)'
   write(*,'(A)') '======================================================================='
   write(*,'(A)') ''
   write(*,'(A)') '  Geographic latitude        90km   110km   200km   300km  1000km'
   write(*,'(A)') '  ----------------------   ------  ------  ------  ------  ------'

   Su(1) = 0.0_dp;  Su(2) = 0.0_dp
   Sa(1) = 0.0_dp
   Rjud_local    = Rjud_ref
   utc_mjd_local = utc_mjd

   do k = 1, N_LATS
      lat_rad = GEO_LAT_DEG(k) * (PI / 180.0_dp)
      Sa(2)   = real(lat_rad, 8)
      geo_lat = GEO_LAT_DEG(k)
      do j = 1, N_KEY
         alt_km   = KEY_ALTS(j)
         Sa(3)    = alt_km * km2m
         ! satellite at RA=0, geocentric lat = lat_rad
         position = [(EARTH_RADIUS + alt_km)*cos(lat_rad), 0.0_dp, &
                     (EARTH_RADIUS + alt_km)*sin(lat_rad)]
         call rdymos_cssi(Sa, Su, Rjud_local, Dafr, Gsti, &
                          Te_old, Ad_old, Wmol_old, Rhod_old, old_status)
         if (old_status /= 0) then
            pct_arr(j) = -999.0_dp
         else
            density_new = jr%density(alt_km, position, sun_vector, geo_lat, utc_mjd_local)
            pct_arr(j)  = 100.0_dp * (density_new - real(Rhod_old, dp)) / real(Rhod_old, dp)
         end if
      end do
      write(*,'(2X,A22,5F8.2)') LAT_DESC(k), (pct_arr(j), j=1,N_KEY)
   end do
   Sa(2) = 0.0_dp;  geo_lat = 0.0_dp   ! restore equatorial

   !==========================================================================
   ! Table 6 : Local-time (hour-angle) sensitivity
   !==========================================================================

   write(*,'(A)') ''
   write(*,'(A)') '======================================================================='
   write(*,'(A)') '  Table 6: Local-time sensitivity   Delta% = 100*(rho_new - rho_old)/rho_old'
   write(*,'(A)') '  Fixed: reference epoch (~2011), equatorial, sun at equinox'
   write(*,'(A)') '  Hour angle h = sat_RA - sun_RA  (sun fixed at RA=0)'
   write(*,'(A)') '======================================================================='
   write(*,'(A)') ''
   write(*,'(A)') '  Local time                 90km   110km   200km   300km  1000km'
   write(*,'(A)') '  ----------------------   ------  ------  ------  ------  ------'

   Su(1) = 0.0_dp;  Su(2) = 0.0_dp
   Sa(2) = 0.0_dp
   Rjud_local    = Rjud_ref
   utc_mjd_local = utc_mjd
   geo_lat       = 0.0_dp

   do k = 1, N_TIMES
      Sa(1)  = SAT_RA_RAD(k)
      sat_ra = real(SAT_RA_RAD(k), dp)
      do j = 1, N_KEY
         alt_km   = KEY_ALTS(j)
         Sa(3)    = alt_km * km2m
         ! equatorial position at satellite RA
         position = [(EARTH_RADIUS + alt_km)*cos(sat_ra), &
                     (EARTH_RADIUS + alt_km)*sin(sat_ra), 0.0_dp]
         call rdymos_cssi(Sa, Su, Rjud_local, Dafr, Gsti, &
                          Te_old, Ad_old, Wmol_old, Rhod_old, old_status)
         if (old_status /= 0) then
            pct_arr(j) = -999.0_dp
         else
            density_new = jr%density(alt_km, position, sun_vector, geo_lat, utc_mjd_local)
            pct_arr(j)  = 100.0_dp * (density_new - real(Rhod_old, dp)) / real(Rhod_old, dp)
         end if
      end do
      write(*,'(2X,A22,5F8.2)') TIME_DESC(k), (pct_arr(j), j=1,N_KEY)
   end do

   ! Clean up
   call jr%destroy()

end program test_compare_roberts_new
