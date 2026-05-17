!------------------------------------------------------------------------------
!>
!   Read and interpolate CSSI Space Weather data files (legacy text format).
!
!   Compatible with CSSI Space Weather files from https://celestrak.org/SpaceData/
!
!   The CSSI Space Weather format contains:
!    * F10.7 solar flux (observed and adjusted)
!    * F10.7a 81-day centered average
!    * Geomagnetic Kp indices (8 values per day, 3-hour periods)
!    * Ap indices (geomagnetic activity)
!
!   Format: `FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,F4.1,I2,I4,F6.1,I2,5F6.1)`
!
!### Reference:
!   * [CSSI Space Weather Data Format Specification](http://celestrak.com/SpaceData/SpaceWx-format.asp)

module space_weather_module
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none
   private

   ! Public interfaces
   public :: sw_init
   public :: sw_get_flux_data
   public :: sw_cleanup
   public :: flux_data_type

   ! Double precision kind
   integer, parameter :: dp = real64
   integer, parameter :: ip = int32

   type :: flux_data_type
      !! Space weather data type
      real(dp) :: mjd                !! Modified Julian Date
      real(dp) :: f107_obs           !! Observed F10.7 solar flux
      real(dp) :: f107_adj           !! Adjusted F10.7 solar flux
      real(dp) :: f107a_obs_ctr      !! Observed F10.7 81-day centered average
      real(dp) :: f107a_adj_ctr      !! Adjusted F10.7 81-day centered average
      real(dp) :: kp(8)              !! Kp indices for 8 3-hour periods
      real(dp) :: ap_avg             !! Average Ap for the day
   end type flux_data_type

   type :: sw_data_type
      !! Module state
      integer(ip) :: n_records = 0
      real(dp), allocatable :: mjd(:)
      real(dp), allocatable :: f107_obs(:)
      real(dp), allocatable :: f107_adj(:)
      real(dp), allocatable :: f107a_obs_ctr(:)
      real(dp), allocatable :: f107a_adj_ctr(:)
      real(dp), allocatable :: kp(:,:)  !! (8, n_records)
      real(dp), allocatable :: ap_avg(:)
      logical :: initialized = .false.
   end type sw_data_type

   type(sw_data_type), save :: sw_data

contains

   !---------------------------------------------------------------------------
   !>
   !   Initialize the space weather module by reading a CSSI file

   subroutine sw_init(filename, status)
      character(len=*), intent(in) :: filename !! Path to CSSI space weather file
      integer(ip), intent(out) :: status !! Output status (0=success, non-zero=error)

      integer(ip) :: unit, io_stat, i
      character(len=500) :: line
      integer(ip) :: year, month, day, bsrn, nd, kp_sum, kp_int(8), ap_int(8), ap_avg_int
      integer(ip) :: c9, isn, qual
      real(dp) :: cp, f107_adj, f107_adj_ctr81, f107_adj_lst81
      real(dp) :: f107_obs, f107_obs_ctr81
      real(dp) :: mjd_val
      integer(ip) :: n_lines
      logical :: in_data_section

      status = 0

      ! Deallocate if already allocated
      if (allocated(sw_data%mjd)) then
         deallocate(sw_data%mjd)
         deallocate(sw_data%f107_obs)
         deallocate(sw_data%f107_adj)
         deallocate(sw_data%f107a_obs_ctr)
         deallocate(sw_data%f107a_adj_ctr)
         deallocate(sw_data%kp)
         deallocate(sw_data%ap_avg)
      end if

      ! Open file
      open(newunit=unit, file=filename, status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) then
         write(*,'(A,A)') 'ERROR: Cannot open space weather file: ', trim(filename)
         status = 1
         return
      end if

      ! First pass: count valid data records
      in_data_section = .false.
      n_lines = 0

      do
         read(unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle

         ! Check for BEGIN sections (OBSERVED, DAILY_PREDICTED, MONTHLY_PREDICTED)
         if (index(line, 'BEGIN OBSERVED') > 0 .or. &
             index(line, 'BEGIN DAILY_PREDICTED') > 0 .or. &
             index(line, 'BEGIN MONTHLY_PREDICTED') > 0) then
            in_data_section = .true.
            cycle
         end if

         ! Check for END sections
         if (index(line, 'END OBSERVED') > 0 .or. &
             index(line, 'END DAILY_PREDICTED') > 0 .or. &
             index(line, 'END MONTHLY_PREDICTED') > 0) then
            in_data_section = .false.
            cycle
         end if

         ! Stop at monthly fit section (like C++ stops at END MONTHLY_FIT)
         if (index(line, 'BEGIN MONTHLY_FIT') > 0) exit

         if (in_data_section .and. len_trim(line) > 80) then
            n_lines = n_lines + 1
         end if
      end do

      if (n_lines == 0) then
         write(*,'(A)') 'ERROR: No data records found in file'
         status = 4
         close(unit)
         return
      end if

      ! Allocate arrays with exact size needed
      allocate(sw_data%mjd(n_lines))
      allocate(sw_data%f107_obs(n_lines))
      allocate(sw_data%f107_adj(n_lines))
      allocate(sw_data%f107a_obs_ctr(n_lines))
      allocate(sw_data%f107a_adj_ctr(n_lines))
      allocate(sw_data%kp(8, n_lines))
      allocate(sw_data%ap_avg(n_lines))

      ! Rewind file for second pass
      rewind(unit)

      ! Second pass: read the data
      in_data_section = .false.
      sw_data%n_records = 0

      do
         read(unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         ! Skip empty lines and comments
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle

         ! Check for BEGIN sections (OBSERVED, DAILY_PREDICTED, MONTHLY_PREDICTED)
         if (index(line, 'BEGIN OBSERVED') > 0 .or. &
             index(line, 'BEGIN DAILY_PREDICTED') > 0 .or. &
             index(line, 'BEGIN MONTHLY_PREDICTED') > 0) then
            in_data_section = .true.
            cycle
         end if

         ! Check for END sections
         if (index(line, 'END OBSERVED') > 0 .or. &
             index(line, 'END DAILY_PREDICTED') > 0 .or. &
             index(line, 'END MONTHLY_PREDICTED') > 0) then
            in_data_section = .false.
            cycle
         end if

         ! Stop at monthly fit section (like C++ stops at END MONTHLY_FIT)
         if (index(line, 'BEGIN MONTHLY_FIT') > 0) exit

         ! Parse data lines using fixed-width format
         if (in_data_section .and. len_trim(line) >= 130) then
            ! FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,F4.1,I2,I4,F6.1,I2,5F6.1)
            ! Column positions per CelesTrak documentation
            read(line, '(I4,1X,I2,1X,I2,1X,I4,1X,I2,1X,8(I2,1X),I3,1X,8(I3,1X),I3,1X,F3.1,1X,I1,1X,I3,1X,F5.1,1X,I1,1X,5(F5.1,1X))', iostat=io_stat) &
                 year, month, day, bsrn, nd, &
                 kp_int(1), kp_int(2), kp_int(3), kp_int(4), &
                 kp_int(5), kp_int(6), kp_int(7), kp_int(8), kp_sum, &
                 ap_int(1), ap_int(2), ap_int(3), ap_int(4), &
                 ap_int(5), ap_int(6), ap_int(7), ap_int(8), ap_avg_int, &
                 cp, c9, isn, f107_adj, qual, f107_adj_ctr81, f107_adj_lst81, &
                 f107_obs, f107_obs_ctr81

            if (io_stat /= 0) cycle

            ! Convert date to MJD
            call date_to_mjd(year, month, day, mjd_val)

            sw_data%n_records = sw_data%n_records + 1

            ! Store data (Kp values are stored as Kp*10, convert back)
            sw_data%mjd(sw_data%n_records) = mjd_val
            do i = 1, 8
               sw_data%kp(i, sw_data%n_records) = real(kp_int(i), dp) / 10.0_dp
            end do
            sw_data%ap_avg(sw_data%n_records) = real(ap_avg_int, dp)
            sw_data%f107_adj(sw_data%n_records) = f107_adj
            sw_data%f107a_adj_ctr(sw_data%n_records) = f107_adj_ctr81
            sw_data%f107_obs(sw_data%n_records) = f107_obs
            sw_data%f107a_obs_ctr(sw_data%n_records) = f107_obs_ctr81
         end if
      end do

      close(unit)

      if (sw_data%n_records == 0) then
         write(*,'(A)') 'ERROR: No valid data records found'
         status = 4
         return
      end if

      sw_data%initialized = .true.

      write(*,'(A,I0,A)') 'Space weather data loaded: ', sw_data%n_records, ' records'
      write(*,'(A,I0,A,I0)') '  (Allocated for ', n_lines, ' lines, successfully parsed ', sw_data%n_records, ')'
      write(*,'(A,F12.2,A,F12.2)') '  Date range (MJD): ', sw_data%mjd(1), ' to ', &
                                    sw_data%mjd(sw_data%n_records)

   end subroutine sw_init

   !---------------------------------------------------------------------------
   !>
   !   Get space weather data for a given Modified Julian Date
   !   Performs linear interpolation if exact date not found

   subroutine sw_get_flux_data(mjd, flux_data, status)
      real(dp), intent(in) :: mjd !! Modified Julian Date
      type(flux_data_type), intent(out) :: flux_data !! Output flux data structure
      integer(ip), intent(out) :: status !! Output status (0=success, 1=not initialized, 2=out of range)

      integer(ip) :: idx
      real(dp) :: frac

      status = 0

      if (.not. sw_data%initialized) then
         write(*,'(A)') 'ERROR: Space weather module not initialized'
         status = 1
         return
      end if

      ! Check if MJD is in range
      if (mjd < sw_data%mjd(1) .or. mjd > sw_data%mjd(sw_data%n_records)) then
         write(*,'(A,F12.2)') 'WARNING: MJD out of range: ', mjd
         write(*,'(A,F12.2,A,F12.2)') '  Valid range: ', sw_data%mjd(1), ' to ', &
                                       sw_data%mjd(sw_data%n_records)
         status = 2
         ! Return data from nearest boundary
         if (mjd < sw_data%mjd(1)) then
            idx = 1
         else
            idx = sw_data%n_records
         end if
         call copy_record(idx, flux_data)
         return
      end if

      ! Find the record (binary search would be better for large datasets)
      do idx = 1, sw_data%n_records
         if (sw_data%mjd(idx) >= mjd) exit
      end do

      ! If exact match or first record
      if (abs(sw_data%mjd(idx) - mjd) < 1.0e-6_dp .or. idx == 1) then
         call copy_record(idx, flux_data)
         return
      end if

      ! Linear interpolation between idx-1 and idx
      frac = (mjd - sw_data%mjd(idx-1)) / (sw_data%mjd(idx) - sw_data%mjd(idx-1))

      flux_data%mjd = mjd
      flux_data%f107_obs = sw_data%f107_obs(idx-1) + frac * &
                           (sw_data%f107_obs(idx) - sw_data%f107_obs(idx-1))
      flux_data%f107_adj = sw_data%f107_adj(idx-1) + frac * &
                           (sw_data%f107_adj(idx) - sw_data%f107_adj(idx-1))
      flux_data%f107a_obs_ctr = sw_data%f107a_obs_ctr(idx-1) + frac * &
                                (sw_data%f107a_obs_ctr(idx) - sw_data%f107a_obs_ctr(idx-1))
      flux_data%f107a_adj_ctr = sw_data%f107a_adj_ctr(idx-1) + frac * &
                                (sw_data%f107a_adj_ctr(idx) - sw_data%f107a_adj_ctr(idx-1))
      flux_data%ap_avg = sw_data%ap_avg(idx-1) + frac * &
                        (sw_data%ap_avg(idx) - sw_data%ap_avg(idx-1))

      ! For Kp, use the values from the closest record (no interpolation)
      if (frac < 0.5_dp) then
         flux_data%kp = sw_data%kp(:, idx-1)
      else
         flux_data%kp = sw_data%kp(:, idx)
      end if

   end subroutine sw_get_flux_data

   !---------------------------------------------------------------------------
   !>
   !   Clean up allocated memory

   subroutine sw_cleanup()
      if (allocated(sw_data%mjd)) then
         deallocate(sw_data%mjd)
         deallocate(sw_data%f107_obs)
         deallocate(sw_data%f107_adj)
         deallocate(sw_data%f107a_obs_ctr)
         deallocate(sw_data%f107a_adj_ctr)
         deallocate(sw_data%kp)
         deallocate(sw_data%ap_avg)
      end if
      sw_data%initialized = .false.
      sw_data%n_records = 0
   end subroutine sw_cleanup

   !---------------------------------------------------------------------------
   !>
   !   Copy a single record from the space weather data

   subroutine copy_record(idx, flux_data)
      integer(ip), intent(in) :: idx !! Index of the record to copy
      type(flux_data_type), intent(out) :: flux_data !! Output flux data structure

      flux_data%mjd = sw_data%mjd(idx)
      flux_data%f107_obs = sw_data%f107_obs(idx)
      flux_data%f107_adj = sw_data%f107_adj(idx)
      flux_data%f107a_obs_ctr = sw_data%f107a_obs_ctr(idx)
      flux_data%f107a_adj_ctr = sw_data%f107a_adj_ctr(idx)
      flux_data%kp = sw_data%kp(:, idx)
      flux_data%ap_avg = sw_data%ap_avg(idx)
   end subroutine copy_record

   !---------------------------------------------------------------------------
   !>
   !   Convert calendar date to Modified Julian Date (MJD)
   subroutine date_to_mjd(year, month, day, mjd)
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

end module space_weather_module
