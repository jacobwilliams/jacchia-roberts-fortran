program gmat_comparison

use csv_module
use jacchia_roberts_module
use jacchia_roberts_kinds, only: ip, dp
use space_weather_module, only: sw_data_type, flux_data_type

character(len=*),parameter :: gmat_file = 'test/JR_density_test_data.csv'
real(dp),parameter :: rad_earth = 6356.766_dp ! Earth polar radius (km)
character(len=*),parameter :: sw_file = 'data/SpaceWeather-All-v1.2.txt' ! space weather file to

type(csv_file) :: f
character(len=30),dimension(:),allocatable :: header
real(dp),dimension(:),allocatable :: mjd, a1mjd, ttmjd, x, y, z, alt, lat, density, sun_x, sun_y, sun_z
real(dp),dimension(:),allocatable :: rho_model, percent_errors
real(dp),dimension(:),allocatable :: rho_model_shifted, percent_errors_shifted
real(dp) :: rho, percent_error, max_error, max_error2
real(dp) :: rho_shifted, percent_error_shifted, max_error_shifted
logical :: status_ok
real(dp),dimension(3) :: r, r_sun_hat
integer :: i, max_i, max_i2 !! counter and top-error rows
type(jacchia_roberts_type) :: jr
integer(ip) :: status
type(sw_data_type) :: sw_debug
type(flux_data_type) :: flux_debug
real(dp) :: tkp_debug, f107_debug, f107a_debug
real(dp) :: frac_epoch, frac_epoch_kp
integer(ip) :: sub_index
logical :: sw_status

real(dp),parameter :: gmat_mjd_offset = 29999.5_dp ! GMAT's MJD is offset from the standard MJD by this amount
real(dp),parameter :: F107_REF_EPOCH = 48407.5_dp  ! MJD for 5/31/91 noon
real(dp),parameter :: weather_time_bias_seconds = -10.0_dp ! probe GMAT/Fortran epoch alignment

! read the file
call f%read(gmat_file,header_row=1,status_ok=status_ok)
if (.not. status_ok) error stop 'error reading file'

! get some data
! Sat.UTCModJulian,Sat.X,Sat.Y,Sat.Z,Sun.EarthMJ2000Eq.X,Sun.EarthMJ2000Eq.Y,Sun.EarthMJ2000Eq.Z,Sat.Earth.Altitude,Sat.Earth.Latitude,Sat.FM.AtmosDensity
call f%get(1,mjd,status_ok);      if (.not. status_ok) error stop 'error reading mjd'
call f%get(2,a1mjd,status_ok);    if (.not. status_ok) error stop 'error reading a1mjd'
call f%get(3,ttmjd,status_ok);    if (.not. status_ok) error stop 'error reading ttmjd'
call f%get(4,x,status_ok);        if (.not. status_ok) error stop 'error reading x'
call f%get(5,y,status_ok);        if (.not. status_ok) error stop 'error reading y'
call f%get(6,z,status_ok);        if (.not. status_ok) error stop 'error reading z'
call f%get(7,sun_x,status_ok);    if (.not. status_ok) error stop 'error reading sun_x'
call f%get(8,sun_y,status_ok);    if (.not. status_ok) error stop 'error reading sun_y'
call f%get(9,sun_z,status_ok);    if (.not. status_ok) error stop 'error reading sun_z'
call f%get(10,alt,status_ok);      if (.not. status_ok) error stop 'error reading alt' ! km
call f%get(11,lat,status_ok);     if (.not. status_ok) error stop 'error reading lat' ! deg
call f%get(12,density,status_ok); if (.not. status_ok) error stop 'error reading density' ! kg/km^3

! initialize the model:
call jr%initialize(rad_earth, sw_file, status)
if (status /= 0_ip) then
   error stop 'Error initializing model'
end if

! mjd = mjd + gmat_mjd_offset ! adjust GMAT/GSCF MJD to standard MJD

! Replicate GMAT's faulty A1→UTC conversion (gives 40.034s instead of 34.034s for 2011):
!   a1utc_sec = floor(-9.24696 + 0.001925 * round(a1mjd_gsfc)) + 0.03437805175781120
! Then: gmat_utc_gsfc = a1mjd_gsfc - a1utc_sec / 86400
mjd = (a1mjd - (floor(-9.24696_dp + 0.001925_dp * anint(a1mjd)) + &
       0.03437805175781120_dp) / 86400.0_dp) + gmat_mjd_offset

density = density * 1.0e-9_dp ! convert kg/km^3 -> kg/m^3
allocate(rho_model(size(mjd)))
allocate(percent_errors(size(mjd)))
allocate(rho_model_shifted(size(mjd)))
allocate(percent_errors_shifted(size(mjd)))

! compare each point to the results from the jacchia-roberts model
write(*,*) ''
write(*,'(*(a16,1x))') 'mjd', 'alt', 'lat', 'gmat', 'jr', '% error'
max_error = -1.0_dp
max_error_shifted = -1.0_dp
max_error2 = -1.0_dp
max_i = 1
max_i2 = 1
do i = 1, size(mjd)
   r = [x(i), y(i), z(i)]
   r_sun_hat = unit([sun_x(i), sun_y(i), sun_z(i)]) ! unit vector from earth to sun
   rho = jr%density(alt(i), r, r_sun_hat, lat(i), mjd(i))  ! kg/m^3
   rho_shifted = jr%density(alt(i), r, r_sun_hat, lat(i), mjd(i) + weather_time_bias_seconds / 86400.0_dp)
   percent_error = abs(rho - density(i)) / density(i) * 100.0_dp
   percent_error_shifted = abs(rho_shifted - density(i)) / density(i) * 100.0_dp
   rho_model(i) = rho
   rho_model_shifted(i) = rho_shifted
   percent_errors(i) = percent_error
   percent_errors_shifted(i) = percent_error_shifted
   if (percent_error_shifted > max_error_shifted) max_error_shifted = percent_error_shifted
   write(*,'(*(E16.9,1x))') mjd(i), alt(i), lat(i), density(i), rho, percent_error
    if (percent_error > max_error) then
        max_error2 = max_error
        max_i2 = max_i
        max_error = percent_error
        max_i = i
    else if (percent_error > max_error2) then
        max_error2 = percent_error
        max_i2 = i
    end if
end do

write(*,*) ''
write(*,'(a,E24.16)') 'Max % error (nominal epoch)      : ', max_error
write(*,'(a,E24.16)') 'Max % error (epoch -10s shifted) : ', max_error_shifted

! Targeted diagnostics for the two worst-error rows.
call sw_debug%initialize(sw_file, status)
if (status /= 0_ip) error stop 'Error initializing debug space weather data'
call print_point_diagnostics(max_i, 'Point A (largest error)')
if (max_i2 /= max_i) then
    call print_point_diagnostics(max_i2, 'Point B (second largest error)')
end if

call sw_debug%destroy()

! destroy the file
call f%destroy()

contains

    subroutine print_point_diagnostics(idx, label)
        integer, intent(in) :: idx
        character(len=*), intent(in) :: label
        integer :: j, j0, j1
        real(dp) :: kp_phase, frac_in_bin, sec_from_prev_edge, sec_to_next_edge
        real(dp) :: frac_epoch_probe, frac_epoch_kp_probe, kp_phase_probe
        real(dp) :: dt_day, rho_probe_minus, rho_probe_plus
        real(dp) :: err_probe_minus, err_probe_plus
        real(dp) :: r_probe(3), sun_probe(3)
        real(dp) :: tkp_probe_minus, tkp_probe_plus
        real(dp) :: f107_probe_minus, f107_probe_plus
        real(dp) :: f107a_probe_minus, f107a_probe_plus
        integer(ip) :: sub_index_probe_minus, sub_index_probe_plus
        type(flux_data_type) :: flux_probe

        call sw_debug%get_flux_data(mjd(idx), flux_debug, sw_status)
        if (.not. sw_status) error stop 'Error reading debug flux data'
        call prepare_flux_data(sw_debug, flux_debug, mjd(idx), tkp_debug, f107_debug, f107a_debug)

        frac_epoch = mjd(idx) - flux_debug%mjd
        frac_epoch_kp = frac_epoch - 6.7_dp / 24.0_dp
        sub_index = floor(frac_epoch_kp * 8.0_dp)
        if (sub_index >= 8) sub_index = 7_ip

        kp_phase = frac_epoch_kp * 8.0_dp
        frac_in_bin = kp_phase - floor(kp_phase)
        sec_from_prev_edge = frac_in_bin * 10800.0_dp
        sec_to_next_edge = (1.0_dp - frac_in_bin) * 10800.0_dp

        ! Probe sensitivity to tiny timestamp shifts across Kp bin boundaries.
        dt_day = 10.0_dp / 86400.0_dp
        r_probe = [x(idx), y(idx), z(idx)]
        sun_probe = unit([sun_x(idx), sun_y(idx), sun_z(idx)])
        rho_probe_minus = jr%density(alt(idx), r_probe, sun_probe, lat(idx), mjd(idx) - dt_day)
        rho_probe_plus  = jr%density(alt(idx), r_probe, sun_probe, lat(idx), mjd(idx) + dt_day)
        err_probe_minus = abs(rho_probe_minus - density(idx)) / density(idx) * 100.0_dp
        err_probe_plus  = abs(rho_probe_plus  - density(idx)) / density(idx) * 100.0_dp

        call sw_debug%get_flux_data(mjd(idx) - dt_day, flux_probe, sw_status)
        if (.not. sw_status) error stop 'Error reading debug flux data at -10s probe'
        call prepare_flux_data(sw_debug, flux_probe, mjd(idx) - dt_day, tkp_probe_minus, f107_probe_minus, f107a_probe_minus)
        frac_epoch_probe = (mjd(idx) - dt_day) - flux_probe%mjd
        frac_epoch_kp_probe = frac_epoch_probe - 6.7_dp / 24.0_dp
        kp_phase_probe = frac_epoch_kp_probe * 8.0_dp
        sub_index_probe_minus = floor(kp_phase_probe)
        if (sub_index_probe_minus >= 8) sub_index_probe_minus = 7_ip

        call sw_debug%get_flux_data(mjd(idx) + dt_day, flux_probe, sw_status)
        if (.not. sw_status) error stop 'Error reading debug flux data at +10s probe'
        call prepare_flux_data(sw_debug, flux_probe, mjd(idx) + dt_day, tkp_probe_plus, f107_probe_plus, f107a_probe_plus)
        frac_epoch_probe = (mjd(idx) + dt_day) - flux_probe%mjd
        frac_epoch_kp_probe = frac_epoch_probe - 6.7_dp / 24.0_dp
        kp_phase_probe = frac_epoch_kp_probe * 8.0_dp
        sub_index_probe_plus = floor(kp_phase_probe)
        if (sub_index_probe_plus >= 8) sub_index_probe_plus = 7_ip

        write(*,*) ''
        write(*,'(a)') 'High-error diagnostics: '//trim(label)
        write(*,'(a,E24.16)') '  % error          : ', percent_errors(idx)
        write(*,'(a,I0)')      '  row index        : ', idx
        write(*,'(a,E24.16)') '  mjd              : ', mjd(idx)
        write(*,'(a,E24.16)') '  flux row mjd     : ', flux_debug%mjd
        write(*,'(a,E24.16)') '  frac_epoch       : ', frac_epoch
        write(*,'(a,E24.16)') '  frac_epoch_kp    : ', frac_epoch_kp
        write(*,'(a,E24.16)') '  kp_phase(8/day)  : ', kp_phase
        write(*,'(a,I0)')      '  kp sub_index     : ', sub_index
        write(*,'(a,E24.16)') '  sec from edge    : ', min(sec_from_prev_edge, sec_to_next_edge)
        write(*,'(a,E24.16)') '  selected Kp      : ', tkp_debug
        write(*,'(a,E24.16)') '  selected F10.7   : ', f107_debug
        write(*,'(a,E24.16)') '  selected F10.7a  : ', f107a_debug
        write(*,'(a,E24.16)') '  shifted -10s %err: ', percent_errors_shifted(idx)
        write(*,'(a,E24.16)') '  probe -10s %err  : ', err_probe_minus
        write(*,'(a,I0)')      '  probe -10s bin   : ', sub_index_probe_minus
        write(*,'(a,E24.16)') '  probe -10s Kp    : ', tkp_probe_minus
        write(*,'(a,E24.16)') '  probe +10s %err  : ', err_probe_plus
        write(*,'(a,I0)')      '  probe +10s bin   : ', sub_index_probe_plus
        write(*,'(a,E24.16)') '  probe +10s Kp    : ', tkp_probe_plus
        if (mjd(idx) < F107_REF_EPOCH) then
            write(*,'(a)')      '  F10.7 boundary   : pre-1991 (5am)'
        else
            write(*,'(a)')      '  F10.7 boundary   : post-1991 (8am)'
        end if

        write(*,'(a)') '  local neighborhood (i-1:i+1):'
        write(*,'(*(a16,1x))') 'row', 'dt[s]', 'mjd', 'gmat', 'jr', '% error', 'jr(-10s)', '%err(-10s)'
        j0 = max(1, idx - 1)
        j1 = min(size(mjd), idx + 1)
        do j = j0, j1
            write(*,'(I16,1x,F16.6,1x,6E24.16,1x)') j, (mjd(j) - mjd(idx)) * 86400.0_dp, mjd(j), density(j), &
                                                     rho_model(j), percent_errors(j), rho_model_shifted(j), &
                                                     percent_errors_shifted(j)
        end do
    end subroutine print_point_diagnostics

    pure function unit(vec) result(u)
        real(dp), intent(in) :: vec(3)
        real(dp) :: u(3)
        u = vec / norm2(vec)
    end function unit
end program gmat_comparison