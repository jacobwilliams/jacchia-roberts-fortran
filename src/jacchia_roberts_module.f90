!------------------------------------------------------------------------------
!>
!   Jacchia-Roberts atmospheric density model.
!
!   The Jacchia-Roberts atmosphere model is valid for altitudes from
!   100 km to 2500 km.
!
!### Reference:
!   * Roberts, C. E., Jr., "An Analytic Model for Upper Atmosphere Densities
!     Based Upon Jacchia's 1970 Models", Celestial Mechanics, Vol. 4,
!     pp. 368-377, 1971.
!   * See GMAT: General Mission Analysis Tool files:
!     `JacchiaRobertsAtmosphere.cpp` and `SolarFluxReader.cpp`.
!   * See also: https://github.com/jacobwilliams/INPE-atmosphere-models

module jacchia_roberts_module
   use jacchia_roberts_kinds, only: ip, dp
   use space_weather_module, only: sw_data_type, flux_data_type
   implicit none
   private

   ! Mathematical constants
   real(dp), parameter :: PI = acos(-1.0_dp) !! Pi constant
   real(dp), parameter :: RAD_PER_DEG = PI / 180.0_dp !! degree to radian conversion factor

   !---------------------------------------------------------------------------
   ! Physical constants
   !---------------------------------------------------------------------------

   real(dp), parameter :: RHO_ZERO = 3.46e-9_dp     !! Low altitude density in g/cm^3
   real(dp), parameter :: TZERO    = 183.0_dp       !! Temperature in degrees Kelvin at height of 90 km
   real(dp), parameter :: G_ZERO   = 9.80665_dp     !! Earth gravitational constant m/s^2
   real(dp), parameter :: GAS_CON  = 8.31432_dp     !! Gas constant (joules/(degK-mole))
   real(dp), parameter :: AVOGADRO = 6.022045e23_dp !! Avogadro's number

   !---------------------------------------------------------------------------
   ! Values to use if the space weather file cannot be loaded
   ! or data cannot be retrieved for the given date
   !---------------------------------------------------------------------------

   real(dp),parameter :: nominalKp    = 3.0_dp   !! Nominal Kp index
   real(dp),parameter :: nominalF107  = 150.0_dp !! Nominal F10.7 value (solar flux units)
   real(dp),parameter :: nominalF107a = 150.0_dp !! Nominal F10.7a value (solar flux units)

   !---------------------------------------------------------------------------
   ! Model constants
   !---------------------------------------------------------------------------

   real(dp), parameter :: CON_C(5) = [ &
      -89284375.0_dp, &
         3542400.0_dp, &
         -52687.5_dp, &
            340.5_dp, &
             -0.8_dp ] !! Constants for series expansion

   real(dp), parameter :: CON_L(5) = [ &
       0.1031445e5_dp, &
       0.2341230e1_dp, &
       0.1579202e-2_dp, &
      -0.1252487e-5_dp, &
       0.2462708e-9_dp ] !! Constants for series expansion

   real(dp), parameter :: MZERO = 28.82678_dp !! Constants required between 90 km and 100 km

   real(dp), parameter :: M_CON(7) = [ &
      -435093.363387_dp, &
        28275.5646391_dp, &
         -765.33466108_dp, &
           11.043387545_dp, &
           -0.08958790995_dp, &
            0.00038737586_dp, &
           -0.000000697444_dp ] !! Constants for series expansion of M(z) function

   real(dp), parameter :: S_CON(6) = [ &
       3144902516.672729_dp, &
       -123774885.4832917_dp, &
          1816141.096520398_dp, &
          -11403.31079489267_dp, &
              24.36498612105595_dp, &
               0.008957502869707995_dp ] !! Constants for series expansion of S(z) function

   real(dp), parameter :: S_BETA(6) = [ &
      -52864482.17910969_dp, &
         -16632.50847336828_dp, &
             -1.308252378125_dp, &
              0.0_dp, &
              0.0_dp, &
              0.0_dp ] !! Constants for series expansion of S(z) function - temperature part

   real(dp), parameter :: OMEGA = -0.94585589_dp !! Constants required for altitudes between 100 km and 125 km

   real(dp), parameter :: ZETA_CON(7) = [ &
       0.1985549e-10_dp, &
      -0.1833490e-14_dp, &
       0.1711735e-17_dp, &
      -0.1021474e-20_dp, &
       0.3727894e-24_dp, &
      -0.7734110e-28_dp, &
       0.7026942e-32_dp ] !! Constants for series expansion

   !> Molecular masses of atmospheric constituents in grams/mole
   real(dp), parameter :: MOL_MASS(6) = [ &
      28.0134_dp, &   ! Nitrogen
      39.948_dp, &    ! Argon
       4.0026_dp, &   ! Helium
      31.9988_dp, &   ! Molecular Oxygen
      15.9994_dp, &   ! Atomic Oxygen
       1.00797_dp ]   ! Hydrogen

   !> Number density divided by Avogadro's number of atmospheric constituents
   real(dp), parameter :: NUM_DENS(5) = [ &
      0.78110_dp, &       ! Nitrogen
      0.93432e-2_dp, &    ! Argon
      0.61471e-5_dp, &    ! Helium
      0.161778_dp, &      ! Molecular Oxygen
      0.95544e-1_dp ]     ! Atomic Oxygen

   !> Polynomial coefficients for constituent densities of each atmospheric gas
   real(dp), parameter :: CON_DEN(5,7) = reshape([ &
      ! Nitrogen
       0.1093155e2_dp, &
       0.1186783e-2_dp, &
      -0.1677341e-5_dp, &
       0.1420228e-8_dp, &
      -0.7139785e-12_dp, &
       0.1969715e-15_dp, &
      -0.2296182e-19_dp, &
      ! Argon
       0.8049405e1_dp, &
       0.2382822e-2_dp, &
      -0.3391366e-5_dp, &
       0.2909714e-8_dp, &
      -0.1481702e-11_dp, &
       0.4127600e-15_dp, &
      -0.4837461e-19_dp, &
      ! Helium
       0.7646886e1_dp, &
      -0.4383486e-3_dp, &
       0.4694319e-6_dp, &
      -0.2894886e-9_dp, &
       0.9451989e-13_dp, &
      -0.1270838e-16_dp, &
       0.0_dp, &
      ! Molecular Oxygen
       0.9924237e1_dp, &
       0.1600311e-2_dp, &
      -0.2274761e-5_dp, &
       0.1938454e-8_dp, &
      -0.9782183e-12_dp, &
       0.2698450e-15_dp, &
      -0.3131808e-19_dp, &
      ! Atomic Oxygen
       0.1097083e2_dp, &
       0.6118742e-4_dp, &
      -0.1165003e-6_dp, &
       0.9239354e-10_dp, &
      -0.3490739e-13_dp, &
       0.5116298e-17_dp, &
       0.0_dp ], [5, 7], order=[2, 1])

   !---------------------------------------------------------------------------
   ! Type definitions
   !---------------------------------------------------------------------------

   type,public :: geoparms_type
      !! Geomagnetic parameters
      real(dp) :: tkp    = 0.0_dp  !! Geomagnetic index Kp
      real(dp) :: xtemp  = 0.0_dp  !! Exospheric temperature (K)
   end type geoparms_type

   type,public :: jacchia_roberts_type
      !! Jacchia-Roberts atmosphere model type
      private
      ! Module state variables
      real(dp) :: cb_polar_radius  = 0.0_dp   !! Central body polar radius (km)
      real(dp) :: cb_polar_squared = 0.0_dp   !! Polar radius squared (km^2)
      real(dp) :: root1            = 0.0_dp   !! Auxiliary temperature root
      real(dp) :: root2            = 0.0_dp   !! Auxiliary temperature root
      real(dp) :: x_root           = 0.0_dp   !! Complex root real part
      real(dp) :: y_root           = 0.0_dp   !! Complex root imaginary part (absolute value)
      real(dp) :: t_infinity       = 0.0_dp   !! Exospheric temperature at infinity
      real(dp) :: tx               = 0.0_dp   !! Temperature at boundary
      real(dp) :: sum              = 0.0_dp   !! Intermediate temperature sum
      type(sw_data_type) :: sw_data  !! Space weather data type (from space_weather_module)
      contains
      private
      procedure, public :: density     => jacchia_roberts_density
      procedure, public :: initialize  => jr_init
      procedure, public :: destroy     => jr_cleanup
      procedure :: exotherm
      procedure :: rho_100
      procedure :: rho_125
      procedure :: rho_high
   end type jacchia_roberts_type

   ! Public procedures for testing
   public :: prepare_flux_data

contains

   !---------------------------------------------------------------------------
   !>
   !   Initialize the Jacchia-Roberts module with central body parameters

   subroutine jr_init(me, cb_polar_radius, filename, status)
      class(jacchia_roberts_type), intent(inout) :: me
      real(dp), intent(in) :: cb_polar_radius !!  Polar radius of central body (km)
      character(len=*), intent(in) :: filename !! Path to CSSI space weather file
      integer(ip), intent(out) :: status  !! Output status (0=success, non-zero=error)

      me%cb_polar_radius = cb_polar_radius
      me%cb_polar_squared = cb_polar_radius * cb_polar_radius

      ! Initialize other state variables
      me%root1 = 0.0_dp
      me%root2 = 0.0_dp
      me%x_root = 0.0_dp
      me%y_root = 0.0_dp
      me%t_infinity = 0.0_dp
      me%tx = 0.0_dp
      me%sum = 0.0_dp

      call me%sw_data%initialize(filename, status)

   end subroutine jr_init

   !---------------------------------------------------------------------------
   !>
   !   Clean up module resources

   subroutine jr_cleanup(me)
      class(jacchia_roberts_type), intent(inout) :: me
      call me%sw_data%destroy()
   end subroutine jr_cleanup

   !---------------------------------------------------------------------------
   !>
   !   Compute atmospheric density using the Jacchia-Roberts model
   !
   !   This version automatically retrieves space weather data (F10.7, Kp)
   !   from the loaded space weather file based on the provided MJD.

   function jacchia_roberts_density(me, height, position, sun_vector, geo_lat, &
                                    utc_mjd) result(density)
      class(jacchia_roberts_type), intent(inout) :: me
      real(dp), intent(in) :: height !! Spacecraft height above reference ellipsoid (km)
      real(dp), intent(in) :: position(3) !! Spacecraft position vector (km, TOD/GCI)
      real(dp), intent(in) :: sun_vector(3) !! Unit vector to Sun (TOD/GCI)
      real(dp), intent(in) :: geo_lat !! Geodetic latitude (degrees)
      real(dp), intent(in) :: utc_mjd !! UTC Modified Julian Date
      real(dp) :: density !! Atmospheric density (kg/m^3)

      real(dp) :: sun_dec !! Sun declination (radians)
      real(dp) :: temperature !! Local exospheric temperature (K)
      real(dp) :: t_500 !! Temperature at 500 km (K)
      real(dp) :: geo_lat_rad !! Geodetic latitude in radians
      real(dp) :: raw_density !! Density before unit conversion (g/cm^3)
      real(dp) :: f107_daily !! Daily F10.7 value (adjusted for measurement time)
      real(dp) :: f107a_avg !! 81-day average F10.7a (adjusted for measurement time)
      real(dp) :: tkp !! Kp index
      type(geoparms_type) :: geo !! Geomagnetic parameters
      type(flux_data_type) :: flux_data !! Space weather flux data
      logical :: sw_status !! Status for space weather data retrieval

      ! sun declination: atan2(Z, sqrt(X^2 + Y^2))
      sun_dec = atan2(sun_vector(3), sqrt(sun_vector(1)**2 + sun_vector(2)**2))

      ! Get space weather data for this date
      call me%sw_data%get_flux_data(utc_mjd, flux_data, sw_status)

      if (.not. sw_status) then
         write(*,*) 'WARNING: Using nominal space weather values'
         ! Use nominal values as fallback
         geo%xtemp = 379.0_dp + 3.24_dp * nominalF107a + 1.3_dp * (nominalF107 - nominalF107a)
         geo%tkp = nominalKp
      else
         ! Prepare flux data with proper timing adjustments.
         ! Kp: 6.7hr lag (matches PrepareKpData). F10.7: previous day (matches PrepareKpData).
         ! F10.7a: detected day (matches PrepareApData intent; see prepare_flux_data docstring).
         call prepare_flux_data(me%sw_data, flux_data, utc_mjd, tkp, f107_daily, f107a_avg)
         ! Calculate exospheric temperature from adjusted F10.7 data
         geo%xtemp = 379.0_dp + 3.24_dp * f107a_avg + 1.3_dp * (f107_daily - f107a_avg)
         ! Set Kp index in geomagnetic parameters:
         geo%tkp = tkp
      end if

      ! Convert geodetic latitude to radians
      geo_lat_rad = geo_lat * RAD_PER_DEG

      ! Compute height-dependent density (matches C++ algorithm exactly)
      if (height <= 90.0_dp) then
         ! At or below 90 km: use constant density
         raw_density = RHO_ZERO
      else if (height < 100.0_dp) then
         ! 90-100 km: use rho_100
         temperature = me%exotherm(position, sun_vector, geo, height, &
                                   sun_dec, geo_lat_rad)
         raw_density = me%rho_100(height, temperature)
      else if (height <= 125.0_dp) then
         ! 100-125 km: use rho_125
         temperature = me%exotherm(position, sun_vector, geo, height, &
                                   sun_dec, geo_lat_rad)
         raw_density = me%rho_125(height, temperature)
      else if (height <= 2500.0_dp) then
         t_500 = me%exotherm(position, sun_vector, geo, 500.0_dp, &
                             sun_dec, geo_lat_rad)
         temperature = me%exotherm(position, sun_vector, geo, height, &
                                   sun_dec, geo_lat_rad)
         raw_density = me%rho_high(height, temperature, t_500, sun_dec, geo_lat_rad)
      else
         raw_density = 0.0_dp
      end if

      ! Apply density correction
      raw_density = raw_density * rho_cor(height, utc_mjd, geo_lat_rad, geo)

      ! Convert from g/cm^3 to kg/m^3
      density = raw_density * 1.0e3_dp

   end function jacchia_roberts_density

   !---------------------------------------------------------------------------
   !>
   !   Compute the temperature of Earth's atmosphere and auxiliary
   !   temperature-related quantities at a given altitude

   function exotherm(me, position, sun_vector, geo, height, sun_dec, geo_lat) &
                     result(exotemp)
      class(jacchia_roberts_type), intent(inout) :: me
      real(dp), intent(in) :: position(3) !! Spacecraft position (km)
      real(dp), intent(in) :: sun_vector(3) !! Sun unit vector
      type(geoparms_type), intent(in) :: geo !! Geomagnetic parameters
      real(dp), intent(in) :: height !! Spacecraft height (km)
      real(dp), intent(in) :: sun_dec !! Sun declination (radians)
      real(dp), intent(in) :: geo_lat !! Geodetic latitude (radians)
      real(dp) :: exotemp !! Local exospheric temperature (K)

      real(dp) :: expkp, hour_angle, cross_denom
      real(dp) :: theta, eta, tau, th22, t1, sun_denom, cos_denom
      real(dp) :: c_star(5)
      real(dp) :: cosAlpha
      integer(ip) :: i

      real(dp), parameter :: error_tolerance = 1.0e-14_dp
      real(dp), parameter :: real_tol = 1.0e-15_dp
      real(dp), parameter :: d37 = 37.0_dp * RAD_PER_DEG  !! 37 degrees in radians
      real(dp), parameter :: d6  = 6.0_dp  * RAD_PER_DEG  !! 6 degrees in radians
      real(dp), parameter :: d43 = 43.0_dp * RAD_PER_DEG  !! 43 degrees in radians
      real(dp), parameter :: a    = 371.6678_dp !! Constant for temperature calculation
      real(dp), parameter :: b    = 0.0518806_dp !! Constant for temperature calculation
      real(dp), parameter :: c    = -294.3505_dp !! Constant for temperature calculation
      real(dp), parameter :: kbar = -0.00216222_dp !! Constant for temperature calculation

      ! Compute hour angle of the sun
      sun_denom = sqrt(sun_vector(1)**2 + sun_vector(2)**2)
      cos_denom = sqrt(position(1)**2 + position(2)**2)

      ! Check for invalid position (zero magnitude in xy-plane)
      if (cos_denom < real_tol) then
         write(*,*) 'ERROR: Position vector has zero magnitude in xy-plane'
         exotemp = TZERO
         return
      end if

      ! Check for degenerate sun vector (sun at coordinate pole, XY magnitude = 0).
      ! For Earth-Sun geometry this cannot happen (sun declination <= 23.4 deg gives
      ! sun_denom >= 0.92), but guard against it to avoid NaN from 0/0.
      if (sun_denom < real_tol) then
         hour_angle = 0.0_dp
      else

      ! Calculate cosine of angle between vectors
      cosAlpha = (sun_vector(1) * position(1) + sun_vector(2) * position(2)) / &
                 (sun_denom * cos_denom)

      ! Compute cross product for sign determination
      cross_denom = abs(sun_vector(1) * position(2) - sun_vector(2) * position(1))

      ! Calculate hour angle with boundary checking
      ! When vectors are nearly collinear, cosAlpha determines the angle directly
      if (cosAlpha >= 1.0_dp - error_tolerance) then
         ! Vectors aligned: hour angle = 0
         hour_angle = 0.0_dp
      else if (cosAlpha <= -1.0_dp + error_tolerance) then
         ! Vectors opposite: hour angle = +/- PI
         ! If cross product is too small (collinear case), default to +PI
         if (cross_denom < real_tol) then
            hour_angle = PI
         else
            hour_angle = sign(1.0_dp, sun_vector(1) * position(2) - &
                             sun_vector(2) * position(1)) * PI
         end if
      else
         ! General case: use cross product to determine sign
         ! If vectors are collinear (shouldn't happen here due to cosAlpha checks),
         ! default to 0
         if (cross_denom < real_tol) then
            hour_angle = 0.0_dp
         else
            hour_angle = sign(1.0_dp, sun_vector(1) * position(2) - &
                         sun_vector(2) * position(1)) * acos(cosAlpha)
         end if
      end if

      end if  ! sun_denom guard

      ! Compute sun and spacecraft position dependent part of temperature
      theta = 0.5_dp * abs(geo_lat + sun_dec)
      eta = 0.5_dp * abs(geo_lat - sun_dec)
      tau = hour_angle - d37 + d6 * sin(hour_angle + d43)
      if (tau < -PI) then
         tau = tau + 2.0_dp * PI
      else if (tau > PI) then
         tau = tau - 2.0_dp * PI
      end if

      th22 = sin(theta)**2.2_dp
      t1 = geo%xtemp * (1.0_dp + 0.3_dp * (th22 + cos(0.5_dp * tau)**3.0_dp * &
                        (cos(eta)**2.2_dp - th22)))
      expkp = exp(geo%tkp)

      ! Compute t_infinity based on altitude
      if (height < 200.0_dp) then
         me%t_infinity = t1 + 14.0_dp * geo%tkp + 0.02_dp * expkp
      else
         me%t_infinity = t1 + 28.0_dp * geo%tkp + 0.03_dp * expkp
      end if

      me%tx = a + b * me%t_infinity + c * exp(kbar * me%t_infinity)

      ! Compute temperature based on altitude regime
      if (height < 125.0_dp) then
         ! Compute height dependent polynomial
         me%sum = CON_C(5)
         do i = 4, 1, -1
            me%sum = CON_C(i) + me%sum * height
         end do

         ! Compute temperature
         exotemp = me%tx + (me%tx - TZERO) * me%sum / 1.500625e6_dp

      else if (height > 125.0_dp) then
         ! Compute temperature dependent polynomial
         me%sum = CON_L(5)
         do i = 4, 1, -1
            me%sum = CON_L(i) + me%sum * me%t_infinity
         end do

         ! Compute temperature
         exotemp = me%t_infinity - (me%t_infinity - me%tx) * &
                   exp(-(me%tx - TZERO) / (me%t_infinity - me%tx) * &
                       (height - 125.0_dp) / 35.0_dp * me%sum / &
                       (me%cb_polar_radius + height))
      else
         exotemp = me%tx
      end if

      ! Compute auxiliary quantities for heights <= 125 km
      if (height <= 125.0_dp) then
         ! Obtain coefficients of polynomial for auxiliary quantities
         c_star(1) = CON_C(1) + 1500625.0_dp * me%tx / (me%tx - TZERO)
         do i = 2, 5
            c_star(i) = CON_C(i)
         end do

         call find_cstar_roots(c_star, me%tx, me%root1, me%root2, me%x_root, me%y_root)
      end if

   end function exotherm

   !---------------------------------------------------------------------------
   !>
   !   Find the two real roots and the complex conjugate root pair of the
   !   degree-4 temperature-profile polynomial c_star(z).
   !
   !### Algorithm:
   !  Simultaneous Newton-Raphson (INPE/Kuga 1985).  Temperature-
   !  dependent initial guesses keep root1 > 125 km and root2 < 100 km for
   !  all physical temperatures, avoiding log(negative) in [[rho_125]].  The
   !  complex pair is recovered analytically via Vieta's formulas.
   !
   !### Reference
   !  * Kuga, h.k. "Reformulacao computacional do modelo de jacchia-roberts para a densidade atmosferica".
   !    inpe, sao jose dos campos, 1985. a ser publicado
   !
   !### See also
   !  * [[roots]] - original routine that was used.
   !
   !@note This replaces the original C++ code's, which was found to be unstable
   !      for some input temperatures (e.g. 300 K).

   subroutine find_cstar_roots(c_star, tx, root1, root2, x_root, y_root)
      real(dp), intent(in)  :: c_star(5) !! Quartic polynomial coefficients `c_star(1)..c_star(5)`
      real(dp), intent(in)  :: tx        !! Temperature at 125 km (K), used for initial guesses
      real(dp), intent(out) :: root1     !! Larger real root (> 125 km)
      real(dp), intent(out) :: root2     !! Smaller real root (< 100 km)
      real(dp), intent(out) :: x_root    !! Real part of complex conjugate root pair
      real(dp), intent(out) :: y_root    !! Imaginary part of complex conjugate root pair

      real(dp) :: temp_norm, pz1, pz2, dpz1, dpz2, r1n, r2n, x2y2
      integer(ip) :: i !! counter

      integer,parameter ::  MAX_ITER = 100 !! 7 in the original code
      real(dp),parameter :: CONVERGENCE_TOL = 1.0e-14_dp !! 1e-7 in the original code

      ! Temperature-dependent initial guesses (valid for all physical T_inf):
      !   root1 ~ 125-170 km (large real root)
      !   root2 ~  55- 65 km (small real root)
      temp_norm = (tx - 300.0_dp) / 200.0_dp
      root1 = 167.77_dp - 3.35_dp * temp_norm
      root2 =  57.34_dp + 7.95_dp * temp_norm

      do i = 1, MAX_ITER
         pz1  = c_star(1) + root1*(c_star(2) + root1*(c_star(3) + root1*(c_star(4) + root1*c_star(5))))
         dpz1 = c_star(2) + root1*(2.0_dp*c_star(3) + root1*(3.0_dp*c_star(4) + root1*4.0_dp*c_star(5)))
         pz2  = c_star(1) + root2*(c_star(2) + root2*(c_star(3) + root2*(c_star(4) + root2*c_star(5))))
         dpz2 = c_star(2) + root2*(2.0_dp*c_star(3) + root2*(3.0_dp*c_star(4) + root2*4.0_dp*c_star(5)))
         r1n = root1 - pz1 / dpz1
         r2n = root2 - pz2 / dpz2
         if (abs(r1n - root1) < CONVERGENCE_TOL .and. abs(r2n - root2) < CONVERGENCE_TOL) exit
         root1 = r1n
         root2 = r2n
      end do
      root1 = r1n
      root2 = r2n

      ! Complex conjugate pair via Vieta's formulas:
      !   sum  of all 4 roots = -c_star(4)/c_star(5) = 340.5/0.8 = 425.625
      !   prod of all 4 roots =  c_star(1)/c_star(5)
      x_root = 0.5_dp * (425.625_dp - root1 - root2)
      x2y2   = c_star(1) / (c_star(5) * root1 * root2)
      y_root = sqrt(abs(x2y2 - x_root**2))

   end subroutine find_cstar_roots

   !---------------------------------------------------------------------------
   !>
   !   Compute density of the atmosphere between 90 and 100 km

   function rho_100(me, height, temperature) result(density)
      class(jacchia_roberts_type), intent(inout) :: me
      real(dp), intent(in) :: height !! Spacecraft altitude (km)
      real(dp), intent(in) :: temperature !! Exospheric temperature (K)
      real(dp) :: density !! Atmospheric density (g/cm^3)

      real(dp) :: m_poly, s_poly, f2, p1, p2, p3, p4, p5, p6, b(6)
      real(dp) :: x_star, u(2), w(2), v, factor_k, roots_2, log_f1
      integer(ip) :: i

      ! Compute M(z) polynomial
      m_poly = M_CON(7)
      do i = 6, 1, -1
         m_poly = m_poly * height + M_CON(i)
      end do

      ! Compute temperature dependent coefficients
      do i = 1, 6
         b(i) = S_CON(i) + S_BETA(i) * me%tx / (me%tx - TZERO)
      end do

      ! Compute functions of auxiliary temperature values
      roots_2 = me%x_root**2 + me%y_root**2
      x_star = -2.0_dp * me%root1 * me%root2 * me%cb_polar_radius * &
               (me%cb_polar_squared + 2.0_dp * me%cb_polar_radius * &
                me%x_root + roots_2)
      v = (me%cb_polar_radius + me%root1) * &
          (me%cb_polar_radius + me%root2) * &
          (me%cb_polar_squared + 2.0_dp * me%cb_polar_radius * &
           me%x_root + roots_2)
      u(1) = (me%root1 - me%root2) * &
             (me%root1 + me%cb_polar_radius)**2 * &
             (me%root1**2 - 2.0_dp * me%root1 * me%x_root + roots_2)
      u(2) = (me%root1 - me%root2) * &
             (me%root2 + me%cb_polar_radius)**2 * &
             (me%root2**2 - 2.0_dp * me%root2 * me%x_root + roots_2)
      w(1) = me%root1 * me%root2 * me%cb_polar_radius * &
             (me%cb_polar_radius + me%root1) * &
             (me%cb_polar_radius + roots_2 / me%root1)
      w(2) = me%root1 * me%root2 * me%cb_polar_radius * &
             (me%cb_polar_radius + me%root2) * &
             (me%cb_polar_radius + roots_2 / me%root2)

      ! Compute S(z) polynomial for z = root1
      s_poly = b(6)
      do i = 5, 1, -1
         s_poly = s_poly * me%root1 + b(i)
      end do
      p2 = s_poly / u(1)

      ! Compute S(z) polynomial for z = root2
      s_poly = b(6)
      do i = 5, 1, -1
         s_poly = s_poly * me%root2 + b(i)
      end do
      p3 = -s_poly / u(2)

      ! Compute S(z) polynomial for z = negative earth radius
      s_poly = b(6)
      do i = 5, 1, -1
         s_poly = -s_poly * me%cb_polar_radius + b(i)
      end do
      p5 = s_poly / v

      ! Compute power of fourth quantity in f1 function
      p4 = (b(1) - me%root1 * me%root2 * me%cb_polar_squared * &
           (b(5) + b(6) * (2.0_dp * me%x_root + me%root1 + &
            me%root2 - me%cb_polar_radius)) + w(1) * p2 + w(2) * p3 - &
            me%root1 * me%root2 * b(6) * me%cb_polar_radius * roots_2 + &
            me%root1 * me%root2 * (me%cb_polar_squared - roots_2) * p5) / x_star

      ! Compute power of first quantity in f1 function
      p1 = b(6) - 2.0_dp * p4 - p3 - p2

      ! Compute p6 factor in f2 function
      p6 = b(5) + b(6) * (2.0_dp * me%x_root + me%root1 + &
           me%root2 - me%cb_polar_radius) - p5 - &
           2.0_dp * (me%x_root + me%cb_polar_radius) * p4 - &
           (me%root2 + me%cb_polar_radius) * p3 - &
           (me%root1 + me%cb_polar_radius) * p2

      ! Compute natural log of f1 function
      log_f1 = p1 * log((height + me%cb_polar_radius) / (90.0_dp + me%cb_polar_radius)) + &
               p2 * log((height - me%root1) / (90.0_dp - me%root1)) + &
               p3 * log((height - me%root2) / (90.0_dp - me%root2)) + &
               p4 * log((height**2 - 2.0_dp * me%x_root * height + roots_2) / &
                       (8100.0_dp - 180.0_dp * me%x_root + roots_2))

      ! Compute f2 function
      f2 = (height - 90.0_dp) * (M_CON(7) + p5 / ((height + me%cb_polar_radius) * &
           (90.0_dp + me%cb_polar_radius))) + &
           p6 * atan(me%y_root * (height - 90.0_dp) / &
           (me%y_root**2 + (height - me%x_root) * (90.0_dp - me%x_root))) / &
           me%y_root

      ! Compute factor_k
      factor_k = -G_ZERO / (GAS_CON * (me%tx - TZERO))

      ! Compute final density
      density = RHO_ZERO * TZERO * m_poly * exp(factor_k * (log_f1 + f2)) / &
                (MZERO * temperature)

   end function rho_100

   !---------------------------------------------------------------------------
   !>
   !   Compute density of the atmosphere between 100 and 125 km

   function rho_125(me, height, temperature) result(density)
      class(jacchia_roberts_type), intent(inout) :: me
      real(dp), intent(in) :: height !! Spacecraft altitude (km)
      real(dp), intent(in) :: temperature !! Exospheric temperature (K)
      real(dp) :: density !! Atmospheric density (g/cm^3)

      real(dp) :: f4, q1, q2, q3, q4, q5, q6
      real(dp) :: rho_prime, x_star, u(2), w(2), v, factor_k, t_100
      real(dp) :: rho_sum, rhoi, roots_2, log_f3
      integer(ip) :: i

      ! Compute base density polynomial
      rho_prime = ZETA_CON(7)
      do i = 6, 1, -1
         rho_prime = rho_prime * me%t_infinity + ZETA_CON(i)
      end do

      ! Compute base temperature
      t_100 = me%tx + OMEGA * (me%tx - TZERO)

      ! Compute functions of auxiliary temperature values
      roots_2 = me%x_root**2 + me%y_root**2
      x_star = -2.0_dp * me%root1 * me%root2 * me%cb_polar_radius * &
               (me%cb_polar_squared + 2.0_dp * me%cb_polar_radius * &
                me%x_root + roots_2)
      v = (me%cb_polar_radius + me%root1) * &
          (me%cb_polar_radius + me%root2) * &
          (me%cb_polar_squared + 2.0_dp * me%cb_polar_radius * &
           me%x_root + roots_2)
      u(1) = (me%root1 - me%root2) * &
             (me%root1 + me%cb_polar_radius)**2 * &
             (me%root1**2 - 2.0_dp * me%root1 * me%x_root + roots_2)
      u(2) = (me%root1 - me%root2) * &
             (me%root2 + me%cb_polar_radius)**2 * &
             (me%root2**2 - 2.0_dp * me%root2 * me%x_root + roots_2)
      w(1) = me%root1 * me%root2 * me%cb_polar_radius * &
             (me%cb_polar_radius + me%root1) * &
             (me%cb_polar_radius + roots_2 / me%root1)
      w(2) = me%root1 * me%root2 * me%cb_polar_radius * &
             (me%cb_polar_radius + me%root2) * &
             (me%cb_polar_radius + roots_2 / me%root2)

      ! Compute powers
      q2 = 1.0_dp / u(1)
      q3 = -1.0_dp / u(2)
      q5 = 1.0_dp / v
      q4 = (1.0_dp + w(1) * q2 + w(2) * q3 + me%root1 * me%root2 * &
           (me%cb_polar_squared - roots_2) * q5) / x_star
      q1 = -2.0_dp * q4 - q3 - q2
      q6 = -q5 - 2.0_dp * (me%x_root + me%cb_polar_radius) * q4 - &
           (me%root2 + me%cb_polar_radius) * q3 - &
           (me%root1 + me%cb_polar_radius) * q2

      ! Compute log of f3 function
      log_f3 = q1 * log((height + me%cb_polar_radius) / (100.0_dp + me%cb_polar_radius)) + &
               q2 * log((height - me%root1) / (100.0_dp - me%root1)) + &
               q3 * log((height - me%root2) / (100.0_dp - me%root2)) + &
               q4 * log((height**2 - 2.0_dp * me%x_root * height + roots_2) / &
                       (1.0e4_dp - 200.0_dp * me%x_root + roots_2))

      ! Compute f4 function
      f4 = (height - 100.0_dp) * q5 / ((height + me%cb_polar_radius) * &
           (100.0_dp + me%cb_polar_radius)) + &
           q6 * atan(me%y_root * (height - 100.0_dp) / &
           (me%y_root**2 + (height - me%x_root) * (100.0_dp - me%x_root))) / &
           me%y_root

      ! Compute factor_k (diffusion scale height coefficient)
      factor_k = -1500625.0_dp * G_ZERO * me%cb_polar_squared / &
                 (GAS_CON * CON_C(5) * (me%tx - TZERO))

      ! Compute mass-dependent sum
      rho_sum = 0.0_dp
      do i = 1, 5
         rhoi = MOL_MASS(i) * NUM_DENS(i) * exp(MOL_MASS(i) * factor_k * (f4 + log_f3))

         if (i == 3) then  ! Helium
            rhoi = rhoi * (t_100 / temperature)**(-0.38_dp)
         end if

         rho_sum = rho_sum + rhoi
      end do

      density = rho_sum * rho_prime * t_100 / temperature

   end function rho_125

   !---------------------------------------------------------------------------
   !>
   !   Compute density of the atmosphere between 125 and 2500 km

   function rho_high(me, height, temperature, t_500, sun_dec, geo_lat) result(density)
      class(jacchia_roberts_type), intent(inout) :: me
      real(dp), intent(in) :: height !! Spacecraft altitude (km)
      real(dp), intent(in) :: temperature !! Exospheric temperature at spacecraft altitude (K)
      real(dp), intent(in) :: t_500 !! Exospheric temperature at altitude of 500 km (K)
      real(dp), intent(in) :: sun_dec !! Declination of the sun in TOD coordinates (radians)
      real(dp), intent(in) :: geo_lat !! Geodetic latitude of spacecraft (radians)
      real(dp) :: density

      real(dp) :: f, log_di, gamma, exp1, r
      real(dp) :: di, polar125
      integer(ip) :: i, j

      real(dp), parameter :: SIN_PI_OVER_4_CUBED = sin(0.25_dp * PI)**3

      density = 0.0_dp
      polar125 = me%cb_polar_radius + 125.0_dp

      do i = 1, 6
         ! Compute constituent density sum for this atmospheric component
         if (i <= 5) then  ! Skip hydrogen (i=6) for initial calculation
            log_di = CON_DEN(i,7)
            do j = 6, 1, -1
               log_di = log_di * me%t_infinity + CON_DEN(i,j)
            end do
            di = 10.0_dp**log_di / AVOGADRO
         end if

         ! Compute second exponent in density expression
         gamma = 35.0_dp * MOL_MASS(i) * G_ZERO * me%cb_polar_squared * &
                 (me%t_infinity - me%tx) / &
                 (GAS_CON * me%sum * me%t_infinity * &
                  (me%tx - TZERO) * polar125)

         ! Compute first exponent in density expression
         exp1 = 1.0_dp + gamma

         ! A factor which is non-unity only for helium
         f = 1.0_dp

         ! Compute corrections for helium
         if (i == 3) then  ! Helium
            exp1 = exp1 - 0.38_dp
            ! Handle sun_dec = 0 case to avoid NaN from 0/0
            if (abs(sun_dec) < 1.0e-15_dp) then
               ! When sun_dec = 0, the formula simplifies to f = 10^0 = 1
               f = 1.0_dp
            else
               f = 4.9914_dp * abs(sun_dec) * &
                  (sin(0.25_dp * PI - 0.5_dp * geo_lat * sun_dec / abs(sun_dec))**3 - &
                  SIN_PI_OVER_4_CUBED) / PI
               f = 10.0_dp**f
            end if
         end if

         ! Add corrections to density, skip hydrogen unless above 500 km
         if (height > 500.0_dp .and. i == 6) then
            r = MOL_MASS(6) * 10.0_dp**(73.13_dp - (39.4_dp - 5.5_dp * log10(t_500)) * &
                log10(t_500)) * (t_500 / temperature)**exp1 * &
                ((me%t_infinity - temperature) / (me%t_infinity - t_500))**gamma / &
                AVOGADRO
            density = density + r
         else if (i <= 5) then
            r = f * MOL_MASS(i) * di * (me%tx / temperature)**exp1 * &
                ((me%t_infinity - temperature) / (me%t_infinity - me%tx))**gamma
            density = density + r
         end if
      end do

   end function rho_high

   !---------------------------------------------------------------------------
   !>
   !   Calculates the density correction factor

   pure function rho_cor(height, utc_mjd, geo_lat, geo) result(correction)
      real(dp), intent(in) :: height !! Spacecraft altitude (km)
      real(dp), intent(in) :: utc_mjd !! UTC Modified Julian Date
      real(dp), intent(in) :: geo_lat !! Geodetic latitude of spacecraft (radians)
      type(geoparms_type), intent(in) :: geo !! Geomagnetic parameters
      real(dp) :: correction !! Correction factor (multiplicative)

      real(dp) :: geo_cor, semian_cor, slat_cor, f, g, day_58, tausa, alpha
      real(dp) :: sin_lat, eta_lat

      ! Compute geomagnetic activity correction
      if (height < 200.0_dp) then
         geo_cor = 0.012_dp * geo%tkp + 0.000012_dp * exp(geo%tkp)
      else
         geo_cor = 0.0_dp
      end if

      ! Compute semiannual variation correction
      f = (5.876e-7_dp * height**2.331_dp + 0.06328_dp) * exp(-0.002868_dp * height)

      ! Compute years since Jan 1, 1958 midnight (MJD 36204.0)
      ! Note: The C++ code uses 6204.5 because GMAT uses GSFC MJD (reference epoch
      ! Jan 5, 1941 noon = MJD 29999.5). In GSFC MJD, Jan 1, 1958 midnight = 6204.5.
      ! This Fortran implementation uses standard MJD where Jan 1, 1958 midnight = 36204.0.
      ! Verified: 36204.0 - 29999.5 = 6204.5
      day_58 = (utc_mjd - 36204.0_dp) / 365.2422_dp

      tausa = day_58 + 0.09544_dp * &
              ((0.5_dp * (1.0_dp + sin(2.0_dp * PI * day_58 + 6.035_dp)))**1.65_dp - 0.5_dp)
      alpha = sin(4.0_dp * PI * tausa + 4.259_dp)
      g = 0.02835_dp + (0.3817_dp + 0.17829_dp * sin(2.0_dp * PI * tausa + 4.137_dp)) * alpha
      semian_cor = f * g

      ! Compute seasonal latitudinal variation
      sin_lat = sin(geo_lat)
      eta_lat = sin(2.0_dp * PI * day_58 + 1.72_dp) * sin_lat * abs(sin_lat)
      slat_cor = 0.014_dp * (height - 90.0_dp) * eta_lat * &
                 exp(-0.0013_dp * (height - 90.0_dp)**2)

      correction = 10.0_dp**(geo_cor + semian_cor + slat_cor)

   end function rho_cor

   !---------------------------------------------------------------------------
   !>
   !   Prepare flux data for the current epoch, accounting for measurement timing.
   !   Kp selection matches the GMAT `PrepareKpData` historic-data branch:
   !
   !     1. Kp: 6.7 hour lag (Vallado & Finkleman)
   !     2. F10.7 daily: previous day (matches PrepareKpData)
   !     3. F10.7a average: detected day (matches PrepareApData / PrepareKpData comment
   !        intent; PrepareKpData's code uses f107index-1 for F10.7a by copy-paste error)

   subroutine prepare_flux_data(sw_data, flux_data, utc_mjd, kp_out, f107_out, f107a_out)
      type(sw_data_type), intent(inout) :: sw_data !! Space weather data manager
      type(flux_data_type), intent(in) :: flux_data !! Current day's flux data
      real(dp), intent(in) :: utc_mjd !! Current epoch (MJD)
      real(dp), intent(out) :: kp_out !! Selected Kp index for current epoch
      real(dp), intent(out) :: f107_out !! Selected F10.7 daily value for current epoch
      real(dp), intent(out) :: f107a_out !! Selected F10.7a average for current epoch

      real(dp) :: frac_epoch !! Fractional day from midnight of flux_data%mjd
      real(dp) :: frac_epoch_kp !! Fractional day for Kp selection (accounts for 6.7 hour lag)
      real(dp) :: f107_offset !! Time of day offset for F10.7 validity boundary
      integer(ip) :: sub_index !! Index for 3-hour Kp period (0-7 for 8 periods per day)
      integer(ip) :: f107_index !! Index for F10.7 selection (0 for current day, -1 for previous day)
      logical :: sw_status !! Status flag for SW data retrieval
      type(flux_data_type) :: flux_data_prev !! Previous day's flux data (for Kp and F10.7 fallback)
      type(flux_data_type) :: flux_data_2days_ago !! Flux data from 2 days ago (for F10.7 fallback)

      real(dp), parameter :: F107_REF_EPOCH = 48407.5_dp  !! MJD for 5/31/91 noon (matches C++ GSFC MJD 18408.0)

      ! Calculate fractional day from midnight
      frac_epoch = utc_mjd - flux_data%mjd

      ! Apply 6.7 hour lag for Kp per Vallado & Finkleman
      frac_epoch_kp = frac_epoch - 6.7_dp / 24.0_dp

      ! Determine which 3-hour period (0-7 for 8 periods per day)
      sub_index = min(floor(frac_epoch_kp * 8.0_dp), 7)

      ! F10.7 is measured at 8pm (5pm before 5/31/91), but the data row is
      ! valid from 8am on the current day to 8am the next day (5am before ref epoch)
      ! So f107_offset represents the DATA VALIDITY BOUNDARY, not measurement time
      if (utc_mjd < F107_REF_EPOCH) then
         f107_offset = 5.0_dp / 24.0_dp  ! 5am validity boundary
      else
         f107_offset = 8.0_dp / 24.0_dp  ! 8am validity boundary
      end if

      ! Determine f107_index based on time of day
      if (frac_epoch < f107_offset) then
         f107_index = -1  ! Before validity boundary - use previous day's index
      else
         f107_index = 0   ! After validity boundary - use current day's index
      end if

      ! ===== KP SELECTION =====
      if (sub_index >= 0) then
         ! Use Kp from current day (Fortran uses 1-based indexing)
         kp_out = flux_data%kp(sub_index + 1)
      else
         ! Need previous day's data
         call sw_data%get_flux_data(flux_data%mjd - 1.0_dp, flux_data_prev, sw_status)

         if (sw_status) then
            ! Use appropriate Kp from previous day
            ! sub_index is negative, so 8 + sub_index gives the right period
            kp_out = flux_data_prev%kp(8 + sub_index + 1)
         else
            ! If we can't get previous day, use first period of current day
            kp_out = flux_data%kp(1)
         end if
      end if

      ! ===== F10.7 DAILY VALUE (always uses previous day) =====
      if (f107_index == 0) then
         ! Current time is after measurement time - use yesterday's F10.7
         call sw_data%get_flux_data(flux_data%mjd - 1.0_dp, flux_data_prev, sw_status)
         if (sw_status) then
            f107_out = flux_data_prev%f107_obs
         else
            f107_out = flux_data%f107_obs  ! Fallback to current
         end if
      else
         ! Current time is before measurement time - use 2 days ago F10.7
         call sw_data%get_flux_data(flux_data%mjd - 2.0_dp, flux_data_2days_ago, sw_status)
         if (sw_status) then
            f107_out = flux_data_2days_ago%f107_obs
         else
            ! Try 1 day ago as fallback
            call sw_data%get_flux_data(flux_data%mjd - 1.0_dp, flux_data_prev, sw_status)
            if (sw_status) then
               f107_out = flux_data_prev%f107_obs
            else
               f107_out = flux_data%f107_obs  ! Last resort
            end if
         end if
      end if

      ! ===== F10.7a AVERAGE (centered 81-day average) =====
      ! if (f107_index == 0) then
      !    ! GMAT bug?: uses previous day's F10.7a even after 8am
      !    call sw_data%get_flux_data(flux_data%mjd - 1.0_dp, flux_data_prev, sw_status)
      !    if (sw_status) then
      !       f107a_out = flux_data_prev%f107a_obs_ctr
      !    else
      !       f107a_out = flux_data%f107a_obs_ctr  ! Fallback
      !    end if
      if (f107_index == 0) then
         ! Current time is after measurement time - use today's F10.7a
         f107a_out = flux_data%f107a_obs_ctr
      else
         ! Current time is before measurement time - use yesterday's F10.7a
         call sw_data%get_flux_data(flux_data%mjd - 1.0_dp, flux_data_prev, sw_status)
         if (sw_status) then
            f107a_out = flux_data_prev%f107a_obs_ctr
         else
            f107a_out = flux_data%f107a_obs_ctr  ! Fallback to current
         end if
      end if

   end subroutine prepare_flux_data

end module jacchia_roberts_module
