!------------------------------------------------------------------------------
! Module: jacchia_roberts_module
!------------------------------------------------------------------------------
! Description:
!   Jacchia-Roberts atmospheric density model
!   Ported from GMAT C++ implementation
!
!   The Jacchia-Roberts atmosphere model is valid for altitudes from
!   100 km to 2500 km.
!
! Reference:
!   Roberts, C. E., Jr., "An Analytic Model for Upper Atmosphere Densities
!   Based Upon Jacchia's 1970 Models", Celestial Mechanics, Vol. 4,
!   pp. 368-377, 1971.
!
! Original C++ code:
!   GMAT: General Mission Analysis Tool
!   Copyright (c) 2002-2026 United States Government
!   Licensed under the Apache License, Version 2.0
!
! Author: Converted to Fortran 2008
! Date: 2026
!------------------------------------------------------------------------------

module jacchia_roberts_module
   use, intrinsic :: iso_fortran_env, only: real64, int32
   use space_weather_module, only: sw_init, sw_get_flux_data, sw_cleanup, flux_data_type
   implicit none
   private

   ! Public interfaces
   public :: jacchia_roberts_density
   public :: jr_init
   public :: jr_load_space_weather
   public :: jr_cleanup
   public :: geoparms_type

   ! Double precision kind
   integer, parameter :: dp = real64
   integer, parameter :: ip = int32

   ! Mathematical constants
   real(dp), parameter :: PI = acos(-1.0_dp)
   real(dp), parameter :: RAD_PER_DEG = PI / 180.0_dp

   !---------------------------------------------------------------------------
   ! Physical constants
   !---------------------------------------------------------------------------

   ! Low altitude density in g/cm^3
   real(dp), parameter :: RHO_ZERO = 3.46e-9_dp

   ! Temperature in degrees Kelvin at height of 90 km
   real(dp), parameter :: TZERO = 183.0_dp

   ! Earth gravitational constant m/s^2
   real(dp), parameter :: G_ZERO = 9.80665_dp

   ! Gas constant (joules/(degK-mole))
   real(dp), parameter :: GAS_CON = 8.31432_dp

   ! Avogadro's number
   real(dp), parameter :: AVOGADRO = 6.022045e23_dp

   !---------------------------------------------------------------------------
   ! Model constants
   !---------------------------------------------------------------------------

   ! Constants for series expansion
   real(dp), parameter :: CON_C(5) = [ &
      -89284375.0_dp, &
         3542400.0_dp, &
         -52687.5_dp, &
            340.5_dp, &
             -0.8_dp ]

   ! Constants for series expansion
   real(dp), parameter :: CON_L(5) = [ &
       0.1031445e5_dp, &
       0.2341230e1_dp, &
       0.1579202e-2_dp, &
      -0.1252487e-5_dp, &
       0.2462708e-9_dp ]

   ! Constants required between 90 km and 100 km
   real(dp), parameter :: MZERO = 28.82678_dp

   ! Constants for series expansion of M(z) function
   real(dp), parameter :: M_CON(7) = [ &
      -435093.363387_dp, &
        28275.5646391_dp, &
         -765.33466108_dp, &
           11.043387545_dp, &
           -0.08958790995_dp, &
            0.00038737586_dp, &
           -0.000000697444_dp ]

   ! Constants for series expansion of S(z) function
   real(dp), parameter :: S_CON(6) = [ &
       3144902516.672729_dp, &
       -123774885.4832917_dp, &
          1816141.096520398_dp, &
          -11403.31079489267_dp, &
              24.36498612105595_dp, &
               0.008957502869707995_dp ]

   ! Constants for series expansion of S(z) function - temperature part
   real(dp), parameter :: S_BETA(6) = [ &
      -52864482.17910969_dp, &
         -16632.50847336828_dp, &
             -1.308252378125_dp, &
              0.0_dp, &
              0.0_dp, &
              0.0_dp ]

   ! Constants required for altitudes between 100 km and 125 km
   real(dp), parameter :: OMEGA = -0.94585589_dp

   ! Constants for series expansion
   real(dp), parameter :: ZETA_CON(7) = [ &
       0.1985549e-10_dp, &
      -0.1833490e-14_dp, &
       0.1711735e-17_dp, &
      -0.1021474e-20_dp, &
       0.3727894e-24_dp, &
      -0.7734110e-28_dp, &
       0.7026942e-32_dp ]

   ! Molecular masses of atmospheric constituents in grams/mole
   real(dp), parameter :: MOL_MASS(6) = [ &
      28.0134_dp, &   ! Nitrogen
      39.948_dp, &    ! Argon
       4.0026_dp, &   ! Helium
      31.9988_dp, &   ! Molecular Oxygen
      15.9994_dp, &   ! Atomic Oxygen
       1.00797_dp ]   ! Hydrogen

   ! Number density divided by Avogadro's number of atmospheric constituents
   real(dp), parameter :: NUM_DENS(5) = [ &
      0.78110_dp, &       ! Nitrogen
      0.93432e-2_dp, &    ! Argon
      0.61471e-5_dp, &    ! Helium
      0.161778_dp, &      ! Molecular Oxygen
      0.95544e-1_dp ]     ! Atomic Oxygen

   ! Polynomial coefficients for constituent densities of each atmospheric gas
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

   ! Geomagnetic parameters
   type :: geoparms_type
      real(dp) :: tkp        ! Geomagnetic index Kp
      real(dp) :: xtemp      ! Exospheric temperature
   end type geoparms_type

   ! Module state variables
   type :: jr_state_type
      real(dp) :: cb_polar_radius      ! Central body polar radius (km)
      real(dp) :: cb_polar_squared     ! Polar radius squared
      real(dp) :: root1                ! Auxiliary temperature root
      real(dp) :: root2                ! Auxiliary temperature root
      real(dp) :: x_root               ! Complex root real part
      real(dp) :: y_root               ! Complex root imaginary part (absolute value)
      real(dp) :: t_infinity           ! Exospheric temperature at infinity
      real(dp) :: tx                   ! Temperature at boundary
      real(dp) :: sum                  ! Intermediate temperature sum
   end type jr_state_type

   ! Module-level state
   type(jr_state_type), save :: jr_state

contains

   !---------------------------------------------------------------------------
   ! Subroutine: jr_init
   !---------------------------------------------------------------------------
   ! Description:
   !   Initialize the Jacchia-Roberts module with central body parameters
   !
   ! Arguments:
   !   cb_polar_radius - Polar radius of central body (km)
   !---------------------------------------------------------------------------
   subroutine jr_init(cb_polar_radius)
      real(dp), intent(in) :: cb_polar_radius

      jr_state%cb_polar_radius = cb_polar_radius
      jr_state%cb_polar_squared = cb_polar_radius * cb_polar_radius

      ! Initialize other state variables
      jr_state%root1 = 0.0_dp
      jr_state%root2 = 0.0_dp
      jr_state%x_root = 0.0_dp
      jr_state%y_root = 0.0_dp
      jr_state%t_infinity = 0.0_dp
      jr_state%tx = 0.0_dp
      jr_state%sum = 0.0_dp

   end subroutine jr_init

   !---------------------------------------------------------------------------
   ! Subroutine: jr_load_space_weather
   !---------------------------------------------------------------------------
   ! Description:
   !   Load space weather data from a CSSI CSV file
   !
   ! Arguments:
   !   filename - Path to CSSI space weather CSV file
   !   status   - Output status (0=success, non-zero=error)
   !---------------------------------------------------------------------------
   subroutine jr_load_space_weather(filename, status)
      character(len=*), intent(in) :: filename
      integer(ip), intent(out) :: status

      call sw_init(filename, status)

   end subroutine jr_load_space_weather

   !---------------------------------------------------------------------------
   ! Subroutine: jr_cleanup
   !---------------------------------------------------------------------------
   ! Description:
   !   Clean up module resources
   !---------------------------------------------------------------------------
   subroutine jr_cleanup()
      call sw_cleanup()
   end subroutine jr_cleanup

   !---------------------------------------------------------------------------
   ! Function: jacchia_roberts_density
   !---------------------------------------------------------------------------
   ! Description:
   !   Compute atmospheric density using the Jacchia-Roberts model
   !
   !   This version automatically retrieves space weather data (F10.7, Kp)
   !   from the loaded space weather file based on the provided MJD.
   !
   ! Arguments:
   !   height       - Spacecraft height above reference ellipsoid (km)
   !   position     - Spacecraft position vector (km, TOD/GCI)
   !   sun_vector   - Unit vector to Sun (TOD/GCI)
   !   geo_lat      - Geodetic latitude (degrees)
   !   sun_dec      - Sun declination (radians)
   !   utc_mjd      - UTC Modified Julian Date
   !
   ! Returns:
   !   density - Atmospheric density (kg/m^3)
   !---------------------------------------------------------------------------
   function jacchia_roberts_density(height, position, sun_vector, geo_lat, &
                                    sun_dec, utc_mjd) result(density)
      real(dp), intent(in) :: height
      real(dp), intent(in) :: position(3)
      real(dp), intent(in) :: sun_vector(3)
      real(dp), intent(in) :: geo_lat
      real(dp), intent(in) :: sun_dec
      real(dp), intent(in) :: utc_mjd
      real(dp) :: density

      real(dp) :: temperature, t_500
      real(dp) :: geo_lat_rad
      real(dp) :: raw_density
      type(geoparms_type) :: geo
      type(flux_data_type) :: flux_data
      integer(ip) :: sw_status

      ! Check altitude validity
      if (height < 100.0_dp) then
         write(*,*) 'ERROR: Jacchia-Roberts model not valid below 100 km'
         density = 0.0_dp
         return
      end if

      ! Get space weather data for this date
      call sw_get_flux_data(utc_mjd, flux_data, sw_status)

      if (sw_status /= 0) then
         write(*,*) 'WARNING: Using nominal space weather values'
         ! Use nominal values as fallback
         geo%xtemp = 379.0_dp + 3.24_dp * 150.0_dp + 1.3_dp * (150.0_dp - 150.0_dp)
         geo%tkp = 3.0_dp
      else
         ! Calculate exospheric temperature from F10.7 data
         ! Formula from GMAT: geo.xtemp = 379.0 + 3.24 * F107a + 1.3 * (F107 - F107a)
         geo%xtemp = 379.0_dp + 3.24_dp * flux_data%f107a_obs_ctr + &
                     1.3_dp * (flux_data%f107_obs - flux_data%f107a_obs_ctr)

         ! Use the first Kp value (3-hour period starting at 00:00 UTC)
         ! For more accuracy, could select Kp based on time of day
         geo%tkp = flux_data%kp(1)
      end if

      ! Convert geodetic latitude to radians
      geo_lat_rad = geo_lat * RAD_PER_DEG

      ! Compute height-dependent density
      if (height < 100.0_dp) then
         raw_density = RHO_ZERO
      else if (height < 125.0_dp) then
         temperature = exotherm(position, sun_vector, geo, height, &
                               sun_dec, geo_lat_rad)
         raw_density = rho_100(height, temperature)
      else if (height <= 2500.0_dp) then
         t_500 = exotherm(position, sun_vector, geo, 500.0_dp, &
                         sun_dec, geo_lat_rad)
         temperature = exotherm(position, sun_vector, geo, height, &
                               sun_dec, geo_lat_rad)
         raw_density = rho_high(height, temperature, t_500, sun_dec, geo_lat_rad)
      else
         raw_density = 0.0_dp
      end if

      ! Apply density correction
      raw_density = raw_density * rho_cor(height, utc_mjd, geo_lat_rad, geo)

      ! Convert from g/cm^3 to kg/m^3
      density = raw_density * 1.0e3_dp

   end function jacchia_roberts_density

   !---------------------------------------------------------------------------
   ! Function: exotherm
   !---------------------------------------------------------------------------
   ! Description:
   !   Compute the temperature of Earth's atmosphere and auxiliary
   !   temperature-related quantities at a given altitude
   !
   ! Arguments:
   !   position     - Spacecraft position (km)
   !   sun_vector   - Sun unit vector
   !   geo          - Geomagnetic parameters
   !   height       - Spacecraft height (km)
   !   sun_dec      - Sun declination (radians)
   !   geo_lat      - Geodetic latitude (radians)
   !
   ! Returns:
   !   exotemp - Local exospheric temperature (K)
   !---------------------------------------------------------------------------
   function exotherm(position, sun_vector, geo, height, sun_dec, geo_lat) &
                     result(exotemp)
      real(dp), intent(in) :: position(3)
      real(dp), intent(in) :: sun_vector(3)
      type(geoparms_type), intent(in) :: geo
      real(dp), intent(in) :: height
      real(dp), intent(in) :: sun_dec
      real(dp), intent(in) :: geo_lat
      real(dp) :: exotemp

      real(dp) :: expkp, hour_angle, cross_denom
      real(dp) :: theta, eta, tau, th22, t1, sun_denom, cos_denom
      real(dp) :: c_star(5), aux(4,2)
      real(dp) :: cosAlpha
      real(dp), parameter :: error_tolerance = 1.0e-14_dp
      real(dp), parameter :: real_tol = 1.0e-15_dp
      integer(ip) :: i, na

      na = 5

      ! Compute hour angle of the sun
      sun_denom = sqrt(sun_vector(1)**2 + sun_vector(2)**2)
      cos_denom = sqrt(position(1)**2 + position(2)**2)

      ! Check for invalid position (zero magnitude in xy-plane)
      if (cos_denom < real_tol) then
         write(*,*) 'ERROR: Position vector has zero magnitude in xy-plane'
         exotemp = TZERO
         return
      end if

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
         ! Vectors opposite: hour angle = ±π
         ! If cross product is too small (collinear case), default to +π
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

      ! Compute sun and spacecraft position dependent part of temperature
      theta = 0.5_dp * abs(geo_lat + sun_dec)
      eta = 0.5_dp * abs(geo_lat - sun_dec)
      tau = hour_angle - 0.64577182325_dp + 0.10471975512_dp * &
            sin(hour_angle + 0.75049157836_dp)

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
         jr_state%t_infinity = t1 + 14.0_dp * geo%tkp + 0.02_dp * expkp
      else
         jr_state%t_infinity = t1 + 28.0_dp * geo%tkp + 0.03_dp * expkp
      end if

      jr_state%tx = 371.6678_dp + 0.0518806_dp * jr_state%t_infinity - &
                    294.3505_dp * exp(-0.00216222_dp * jr_state%t_infinity)

      ! Compute temperature based on altitude regime
      if (height < 125.0_dp) then
         ! Compute height dependent polynomial
         jr_state%sum = CON_C(5)
         do i = 4, 1, -1
            jr_state%sum = CON_C(i) + jr_state%sum * height
         end do

         ! Compute temperature
         exotemp = jr_state%tx + (jr_state%tx - TZERO) * jr_state%sum / 1.500625e6_dp

      else if (height > 125.0_dp) then
         ! Compute temperature dependent polynomial
         jr_state%sum = CON_L(5)
         do i = 4, 1, -1
            jr_state%sum = CON_L(i) + jr_state%sum * jr_state%t_infinity
         end do

         ! Compute temperature
         exotemp = jr_state%t_infinity - (jr_state%t_infinity - jr_state%tx) * &
                   exp(-(jr_state%tx - TZERO) / (jr_state%t_infinity - jr_state%tx) * &
                       (height - 125.0_dp) / 35.0_dp * jr_state%sum / &
                       (jr_state%cb_polar_radius + height))
      else
         exotemp = jr_state%tx
      end if

      ! Compute auxiliary quantities for heights <= 125 km
      if (height <= 125.0_dp) then
         ! Obtain coefficients of polynomial for auxiliary quantities
         c_star(1) = CON_C(1) + 1500625.0_dp * jr_state%tx / (jr_state%tx - TZERO)
         do i = 2, 5
            c_star(i) = CON_C(i)
         end do

         ! Get 1st real root
         aux(1,1) = 125.0_dp
         aux(1,2) = 0.0_dp
         call roots(c_star, na, aux, 1)
         jr_state%root1 = aux(1,1)

         call deflate_polynomial(c_star, na, jr_state%root1, c_star)

         ! Get 2nd real root
         aux(1,1) = 200.0_dp
         aux(1,2) = 0.0_dp
         call roots(c_star, na-1, aux, 1)
         jr_state%root2 = aux(1,1)

         call deflate_polynomial(c_star, na-1, jr_state%root2, c_star)

         ! Get remaining complex roots
         aux(1,1) = 10.0_dp
         aux(1,2) = 125.0_dp
         call roots(c_star, na-2, aux, 1)
         jr_state%x_root = aux(1,1)
         jr_state%y_root = abs(aux(1,2))
      end if

   end function exotherm

   !---------------------------------------------------------------------------
   ! Function: rho_100
   !---------------------------------------------------------------------------
   ! Description:
   !   Compute density of the atmosphere between 90 and 100 km
   !
   ! Arguments:
   !   height      - Spacecraft altitude (km)
   !   temperature - Exospheric temperature (K)
   !
   ! Returns:
   !   Density in g/cm^3
   !---------------------------------------------------------------------------
   function rho_100(height, temperature) result(density)
      real(dp), intent(in) :: height
      real(dp), intent(in) :: temperature
      real(dp) :: density

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
         b(i) = S_CON(i) + S_BETA(i) * jr_state%tx / (jr_state%tx - TZERO)
      end do

      ! Compute functions of auxiliary temperature values
      roots_2 = jr_state%x_root**2 + jr_state%y_root**2
      x_star = -2.0_dp * jr_state%root1 * jr_state%root2 * jr_state%cb_polar_radius * &
               (jr_state%cb_polar_squared + 2.0_dp * jr_state%cb_polar_radius * &
                jr_state%x_root + roots_2)
      v = (jr_state%cb_polar_radius + jr_state%root1) * &
          (jr_state%cb_polar_radius + jr_state%root2) * &
          (jr_state%cb_polar_squared + 2.0_dp * jr_state%cb_polar_radius * &
           jr_state%x_root + roots_2)
      u(1) = (jr_state%root1 - jr_state%root2) * &
             (jr_state%root1 + jr_state%cb_polar_radius)**2 * &
             (jr_state%root1**2 - 2.0_dp * jr_state%root1 * jr_state%x_root + roots_2)
      u(2) = (jr_state%root1 - jr_state%root2) * &
             (jr_state%root2 + jr_state%cb_polar_radius)**2 * &
             (jr_state%root2**2 - 2.0_dp * jr_state%root2 * jr_state%x_root + roots_2)
      w(1) = jr_state%root1 * jr_state%root2 * jr_state%cb_polar_radius * &
             (jr_state%cb_polar_radius + jr_state%root1) * &
             (jr_state%cb_polar_radius + roots_2 / jr_state%root1)
      w(2) = jr_state%root1 * jr_state%root2 * jr_state%cb_polar_radius * &
             (jr_state%cb_polar_radius + jr_state%root2) * &
             (jr_state%cb_polar_radius + roots_2 / jr_state%root2)

      ! Compute S(z) polynomial for z = root1
      s_poly = b(6)
      do i = 5, 1, -1
         s_poly = s_poly * jr_state%root1 + b(i)
      end do
      p2 = s_poly / u(1)

      ! Compute S(z) polynomial for z = root2
      s_poly = b(6)
      do i = 5, 1, -1
         s_poly = s_poly * jr_state%root2 + b(i)
      end do
      p3 = -s_poly / u(2)

      ! Compute S(z) polynomial for z = negative earth radius
      s_poly = b(6)
      do i = 5, 1, -1
         s_poly = -s_poly * jr_state%cb_polar_radius + b(i)
      end do
      p5 = s_poly / v

      ! Compute power of fourth quantity in f1 function
      p4 = (b(1) - jr_state%root1 * jr_state%root2 * jr_state%cb_polar_squared * &
           (b(5) + b(6) * (2.0_dp * jr_state%x_root + jr_state%root1 + &
            jr_state%root2 - jr_state%cb_polar_radius)) + w(1) * p2 + w(2) * p3 - &
            jr_state%root1 * jr_state%root2 * b(6) * jr_state%cb_polar_radius * roots_2 + &
            jr_state%root1 * jr_state%root2 * (jr_state%cb_polar_squared - roots_2) * p5) / x_star

      ! Compute power of first quantity in f1 function
      p1 = b(6) - 2.0_dp * p4 - p3 - p2

      ! Compute p6 factor in f2 function
      p6 = b(5) + b(6) * (2.0_dp * jr_state%x_root + jr_state%root1 + &
           jr_state%root2 - jr_state%cb_polar_radius) - p5 - &
           2.0_dp * (jr_state%x_root + jr_state%cb_polar_radius) * p4 - &
           (jr_state%root2 + jr_state%cb_polar_radius) * p3 - &
           (jr_state%root1 + jr_state%cb_polar_radius) * p2

      ! Compute natural log of f1 function
      log_f1 = p1 * log((height + jr_state%cb_polar_radius) / (90.0_dp + jr_state%cb_polar_radius)) + &
               p2 * log((height - jr_state%root1) / (90.0_dp - jr_state%root1)) + &
               p3 * log((height - jr_state%root2) / (90.0_dp - jr_state%root2)) + &
               p4 * log((height**2 - 2.0_dp * jr_state%x_root * height + roots_2) / &
                       (8100.0_dp - 180.0_dp * jr_state%x_root + roots_2))

      ! Compute f2 function
      f2 = (height - 90.0_dp) * (M_CON(7) + p5 / ((height + jr_state%cb_polar_radius) * &
           (90.0_dp + jr_state%cb_polar_radius))) + &
           p6 * atan(jr_state%y_root * (height - 90.0_dp) / &
           (jr_state%y_root**2 + (height - jr_state%x_root) * (90.0_dp - jr_state%x_root))) / &
           jr_state%y_root

      ! Compute factor_k
      factor_k = -G_ZERO / (GAS_CON * (jr_state%tx - TZERO))

      ! Compute final density
      density = RHO_ZERO * TZERO * m_poly * exp(factor_k * (log_f1 + f2)) / &
                (MZERO * temperature)

   end function rho_100

   !---------------------------------------------------------------------------
   ! Function: rho_125
   !---------------------------------------------------------------------------
   ! Description:
   !   Compute density of the atmosphere between 100 and 125 km
   !
   ! Arguments:
   !   height      - Spacecraft altitude (km)
   !   temperature - Exospheric temperature (K)
   !
   ! Returns:
   !   Density in g/cm^3
   !---------------------------------------------------------------------------
   function rho_125(height, temperature) result(density)
      real(dp), intent(in) :: height
      real(dp), intent(in) :: temperature
      real(dp) :: density

      real(dp) :: f4, q1, q2, q3, q4, q5, q6
      real(dp) :: rho_prime, x_star, u(2), w(2), v, factor_k, t_100
      real(dp) :: rho_sum, rhoi, roots_2, log_f3
      integer(ip) :: i

      ! Compute base density polynomial
      rho_prime = ZETA_CON(7)
      do i = 6, 1, -1
         rho_prime = rho_prime * jr_state%t_infinity + ZETA_CON(i)
      end do

      ! Compute base temperature
      t_100 = jr_state%tx + OMEGA * (jr_state%tx - TZERO)

      ! Compute functions of auxiliary temperature values
      roots_2 = jr_state%x_root**2 + jr_state%y_root**2
      x_star = -2.0_dp * jr_state%root1 * jr_state%root2 * jr_state%cb_polar_radius * &
               (jr_state%cb_polar_squared + 2.0_dp * jr_state%cb_polar_radius * &
                jr_state%x_root + roots_2)
      v = (jr_state%cb_polar_radius + jr_state%root1) * &
          (jr_state%cb_polar_radius + jr_state%root2) * &
          (jr_state%cb_polar_squared + 2.0_dp * jr_state%cb_polar_radius * &
           jr_state%x_root + roots_2)
      u(1) = (jr_state%root1 - jr_state%root2) * &
             (jr_state%root1 + jr_state%cb_polar_radius)**2 * &
             (jr_state%root1**2 - 2.0_dp * jr_state%root1 * jr_state%x_root + roots_2)
      u(2) = (jr_state%root1 - jr_state%root2) * &
             (jr_state%root2 + jr_state%cb_polar_radius)**2 * &
             (jr_state%root2**2 - 2.0_dp * jr_state%root2 * jr_state%x_root + roots_2)
      w(1) = jr_state%root1 * jr_state%root2 * jr_state%cb_polar_radius * &
             (jr_state%cb_polar_radius + jr_state%root1) * &
             (jr_state%cb_polar_radius + roots_2 / jr_state%root1)
      w(2) = jr_state%root1 * jr_state%root2 * jr_state%cb_polar_radius * &
             (jr_state%cb_polar_radius + jr_state%root2) * &
             (jr_state%cb_polar_radius + roots_2 / jr_state%root2)

      ! Compute powers
      q2 = 1.0_dp / u(1)
      q3 = -1.0_dp / u(2)
      q5 = 1.0_dp / v
      q4 = (1.0_dp + w(1) * q2 + w(2) * q3 + jr_state%root1 * jr_state%root2 * &
           (jr_state%cb_polar_squared - roots_2) * q5) / x_star
      q1 = -2.0_dp * q4 - q3 - q2
      q6 = -q5 - 2.0_dp * (jr_state%x_root + jr_state%cb_polar_radius) * q4 - &
           (jr_state%root2 + jr_state%cb_polar_radius) * q3 - &
           (jr_state%root1 + jr_state%cb_polar_radius) * q2

      ! Compute log of f3 function
      log_f3 = q1 * log((height + jr_state%cb_polar_radius) / (100.0_dp + jr_state%cb_polar_radius)) + &
               q2 * log((height - jr_state%root1) / (100.0_dp - jr_state%root1)) + &
               q3 * log((height - jr_state%root2) / (100.0_dp - jr_state%root2)) + &
               q4 * log((height**2 - 2.0_dp * jr_state%x_root * height + roots_2) / &
                       (1.0e4_dp - 200.0_dp * jr_state%x_root + roots_2))

      ! Compute f4 function
      f4 = (height - 100.0_dp) * q5 / ((height + jr_state%cb_polar_radius) * &
           (100.0_dp + jr_state%cb_polar_radius)) + &
           q6 * atan(jr_state%y_root * (height - 100.0_dp) / &
           (jr_state%y_root**2 + (height - jr_state%x_root) * (100.0_dp - jr_state%x_root))) / &
           jr_state%y_root

      ! Compute f3 power
      factor_k = -1500625.0_dp * G_ZERO * jr_state%cb_polar_squared / &
                 (GAS_CON * CON_C(5) * (jr_state%tx - TZERO))

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
   ! Function: rho_high
   !---------------------------------------------------------------------------
   ! Description:
   !   Compute density of the atmosphere between 125 and 2500 km
   !
   ! Arguments:
   !   height      - Spacecraft altitude (km)
   !   temperature - Exospheric temperature (K)
   !   t_500       - Temperature at 500 km (K)
   !   sun_dec     - Sun declination (radians)
   !   geo_lat     - Geodetic latitude (radians)
   !
   ! Returns:
   !   Density in g/cm^3
   !---------------------------------------------------------------------------
   function rho_high(height, temperature, t_500, sun_dec, geo_lat) result(density)
      real(dp), intent(in) :: height
      real(dp), intent(in) :: temperature
      real(dp), intent(in) :: t_500
      real(dp), intent(in) :: sun_dec
      real(dp), intent(in) :: geo_lat
      real(dp) :: density

      real(dp) :: f, log_di, gamma, exp1, r
      real(dp) :: di, polar125
      integer(ip) :: i, j

      density = 0.0_dp
      polar125 = jr_state%cb_polar_radius + 125.0_dp

      do i = 1, 6
         ! Compute constituent density sum for this atmospheric component
         if (i <= 5) then  ! Skip hydrogen (i=6) for initial calculation
            log_di = CON_DEN(i,7)
            do j = 6, 1, -1
               log_di = log_di * jr_state%t_infinity + CON_DEN(i,j)
            end do
            di = 10.0_dp**log_di / AVOGADRO
         end if

         ! Compute second exponent in density expression
         gamma = 35.0_dp * MOL_MASS(i) * G_ZERO * jr_state%cb_polar_squared * &
                 (jr_state%t_infinity - jr_state%tx) / &
                 (GAS_CON * jr_state%sum * jr_state%t_infinity * &
                  (jr_state%tx - TZERO) * polar125)

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
                    0.35355_dp) / PI
               f = 10.0_dp**f
            end if
         end if

         ! Add corrections to density, skip hydrogen unless above 500 km
         if (height > 500.0_dp .and. i == 6) then
            r = MOL_MASS(6) * 10.0_dp**(73.13_dp - (39.4_dp - 5.5_dp * log10(t_500)) * &
                log10(t_500)) * (t_500 / temperature)**exp1 * &
                ((jr_state%t_infinity - temperature) / (jr_state%t_infinity - t_500))**gamma / &
                AVOGADRO
            density = density + r
         else if (i <= 5) then
            r = f * MOL_MASS(i) * di * (jr_state%tx / temperature)**exp1 * &
                ((jr_state%t_infinity - temperature) / (jr_state%t_infinity - jr_state%tx))**gamma
            density = density + r
         end if
      end do

   end function rho_high

   !---------------------------------------------------------------------------
   ! Function: rho_cor
   !---------------------------------------------------------------------------
   ! Description:
   !   Calculates the density correction factor
   !
   ! Arguments:
   !   height  - Spacecraft altitude (km)
   !   utc_mjd - UTC Modified Julian Date
   !   geo_lat - Geodetic latitude (radians)
   !   geo     - Geomagnetic parameters
   !
   ! Returns:
   !   Correction factor (multiplicative)
   !---------------------------------------------------------------------------
   function rho_cor(height, utc_mjd, geo_lat, geo) result(correction)
      real(dp), intent(in) :: height
      real(dp), intent(in) :: utc_mjd
      real(dp), intent(in) :: geo_lat
      type(geoparms_type), intent(in) :: geo
      real(dp) :: correction

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
      day_58 = (utc_mjd - 6204.5_dp) / 365.2422_dp
      tausa = day_58 + 0.09544_dp * &
              (0.5_dp * (1.0_dp + sin(2.0_dp * PI * day_58 + 6.035_dp))**1.65_dp - 0.5_dp)
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
   ! Subroutine: roots
   !---------------------------------------------------------------------------
   ! Description:
   !   Finds the roots of a polynomial using Newton's method
   !
   ! Arguments:
   !   a      - Array of polynomial coefficients (lowest power first)
   !   na     - Number of coefficients
   !   croots - Initial guesses and output roots (real, imaginary)
   !   irl    - Number of roots to solve for
   !---------------------------------------------------------------------------
   subroutine roots(a, na, croots, irl)
      real(dp), intent(in) :: a(:)
      integer(ip), intent(in) :: na
      real(dp), intent(inout) :: croots(:,:)
      integer(ip), intent(in) :: irl

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
   ! Subroutine: deflate_polynomial
   !---------------------------------------------------------------------------
   ! Description:
   !   Reduces the order of a polynomial by division
   !
   ! Arguments:
   !   c     - Polynomial coefficients
   !   n     - Order + 1 of polynomial
   !   root  - A single real root of the polynomial
   !   c_new - Output array with new coefficients
   !---------------------------------------------------------------------------
   subroutine deflate_polynomial(c, n, root, c_new)
      real(dp), intent(in) :: c(:)
      integer(ip), intent(in) :: n
      real(dp), intent(in) :: root
      real(dp), intent(inout) :: c_new(:)

      integer(ip) :: i
      real(dp) :: sum, save

      sum = c(n)

      do i = n - 1, 1, -1
         save = c(i)
         c_new(i) = sum
         sum = save + sum * root
      end do

   end subroutine deflate_polynomial

end module jacchia_roberts_module
