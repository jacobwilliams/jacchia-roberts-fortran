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
   use, intrinsic :: iso_fortran_env, only: dp => real64, ip => int32
   use jacchia_roberts_utilities, only: date_to_mjd
   implicit none
   private

   type,public :: flux_data_type
      !! Space weather data type
      real(dp) :: mjd            = 0.0_dp   !! Modified Julian Date
      real(dp) :: f107_obs       = 0.0_dp   !! Observed F10.7 solar flux
      real(dp) :: f107_adj       = 0.0_dp   !! Adjusted F10.7 solar flux
      real(dp) :: f107a_obs_ctr  = 0.0_dp   !! Observed F10.7 81-day centered average
      real(dp) :: f107a_adj_ctr  = 0.0_dp   !! Adjusted F10.7 81-day centered average
      real(dp) :: kp(8)          = 0.0_dp   !! Kp indices for 8 3-hour periods
      real(dp) :: ap_avg         = 0.0_dp   !! Average Ap for the day
   end type flux_data_type

   type,public :: sw_data_type
      !! Space weather data management type
      private
      integer(ip) :: n_records = 0
      real(dp), allocatable :: mjd(:)
      real(dp), allocatable :: f107_obs(:)
      real(dp), allocatable :: f107_adj(:)
      real(dp), allocatable :: f107a_obs_ctr(:)
      real(dp), allocatable :: f107a_adj_ctr(:)
      real(dp), allocatable :: kp(:,:)  !! (8, n_records)
      real(dp), allocatable :: ap_avg(:)
      logical :: initialized = .false.
      real(dp) :: historic_start = -1.0_dp     !! First epoch in data
      real(dp) :: historic_end = -1.0_dp       !! Last daily epoch + 1
      real(dp) :: historic_daily_end = -1.0_dp !! Last contiguous daily data epoch
      logical :: warn_epoch_before = .true.
      logical :: warn_epoch_after = .true.
      contains
      private
      procedure,public :: initialize    => sw_init
      procedure,public :: get_flux_data => sw_get_flux_data
      procedure,public :: destroy       => sw_cleanup
      procedure :: copy_record
   end type sw_data_type

contains

   !---------------------------------------------------------------------------
   !>
   !   Initialize the space weather module by reading a CSSI file

   subroutine sw_init(me, filename, status)
      class(sw_data_type), intent(inout) :: me
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
      if (allocated(me%mjd)) then
         deallocate(me%mjd)
         deallocate(me%f107_obs)
         deallocate(me%f107_adj)
         deallocate(me%f107a_obs_ctr)
         deallocate(me%f107a_adj_ctr)
         deallocate(me%kp)
         deallocate(me%ap_avg)
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

         ! Stop at monthly fit section
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
      allocate(me%mjd(n_lines))
      allocate(me%f107_obs(n_lines))
      allocate(me%f107_adj(n_lines))
      allocate(me%f107a_obs_ctr(n_lines))
      allocate(me%f107a_adj_ctr(n_lines))
      allocate(me%kp(8, n_lines))
      allocate(me%ap_avg(n_lines))

      ! Rewind file for second pass
      rewind(unit)

      ! Second pass: read the data
      in_data_section = .false.
      me%n_records = 0

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

         ! Stop at monthly fit section
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
            ! Space weather files use noon (add 0.5 to midnight MJD)
            call date_to_mjd(year, month, day, mjd_val)
            mjd_val = mjd_val + 0.5_dp

            me%n_records = me%n_records + 1

            ! Store data (Kp values are stored as Kp*10, convert back)
            me%mjd(me%n_records) = mjd_val
            do i = 1, 8
               me%kp(i, me%n_records) = real(kp_int(i), dp) / 10.0_dp
            end do
            me%ap_avg(me%n_records) = real(ap_avg_int, dp)
            me%f107_adj(me%n_records) = f107_adj
            me%f107a_adj_ctr(me%n_records) = f107_adj_ctr81
            me%f107_obs(me%n_records) = f107_obs
            me%f107a_obs_ctr(me%n_records) = f107_obs_ctr81
         end if
      end do

      close(unit)

      if (me%n_records == 0) then
         write(*,'(A)') 'ERROR: No valid data records found'
         status = 4
         return
      end if

      ! Set epoch range tracking
      me%historic_start = me%mjd(1)
      ! Note: add 1.0 to last epoch to reach end of day
      me%historic_end = me%mjd(me%n_records) + 1.0_dp

      ! Find where daily data ends and monthly data begins (gaps > 1 day)
      me%historic_daily_end = me%historic_end
      do i = 2, me%n_records
         if (me%mjd(i) - me%mjd(i-1) > 1.01_dp) then
            me%historic_daily_end = me%mjd(i-1)
            exit
         end if
      end do

      me%initialized = .true.

      write(*,'(A,I0,A)') 'Space weather data loaded: ', me%n_records, ' records'
      write(*,'(A,I0,A,I0)') '  (Allocated for ', n_lines, ' lines, successfully parsed ', me%n_records, ')'
      write(*,'(A,F12.2,A,F12.2)') '  Date range (MJD): ', me%mjd(1), ' to ', &
                                    me%mjd(me%n_records)

   end subroutine sw_init

   !---------------------------------------------------------------------------
   !>
   !   Get space weather data for a given Modified Julian Date.
   !   Uses direct indexing for daily data, no interpolation

   subroutine sw_get_flux_data(me, mjd, flux_data, status)
      class(sw_data_type), intent(inout) :: me
      real(dp), intent(in) :: mjd !! Modified Julian Date
      type(flux_data_type), intent(out) :: flux_data !! Output flux data structure
      integer(ip), intent(out) :: status !! Output status (0=success, 1=not initialized, 2=out of range)

      integer(ip) :: idx, i
      integer(ip) :: daily_end_idx

      status = 0

      if (.not. me%initialized) then
         write(*,'(A)') 'ERROR: Space weather module not initialized'
         status = 1
         return
      end if

      ! Handle epoch before data starts
      if (mjd < me%historic_start) then
         if (me%warn_epoch_before) then
            write(*,'(A)') 'WARNING: Requested epoch is earlier than space weather data start'
            write(*,'(A)') '         Using first file entry'
            me%warn_epoch_before = .false.
         end if
         call me%copy_record(1, flux_data)
         return
      end if

      ! Handle epoch after data ends
      if (mjd >= me%historic_end) then
         if (me%warn_epoch_after) then
            write(*,'(A)') 'WARNING: Requested epoch is later than space weather data end'
            write(*,'(A)') '         Using last file entry'
            me%warn_epoch_after = .false.
         end if
         call me%copy_record(me%n_records, flux_data)
         return
      end if

      ! Within data range
      ! For daily data: direct index calculation (O(1))
      ! For monthly data: linear search through monthly section only

      if (mjd <= me%historic_daily_end) then
         ! Daily data: use direct index calculation
         idx = int(mjd - me%historic_start) + 1

         ! Bounds check
         if (idx < 1) idx = 1
         if (idx > me%n_records) idx = me%n_records

         call me%copy_record(idx, flux_data)
      else
         ! Monthly data: search from daily_end to end of array
         ! Find the index where daily data ends
         daily_end_idx = int(me%historic_daily_end - me%historic_start) + 1

         ! Search monthly data for the record at or before mjd
         idx = daily_end_idx
         do i = daily_end_idx, me%n_records
            if (me%mjd(i) > mjd) exit
            idx = i
         end do

         call me%copy_record(idx, flux_data)
      end if

   end subroutine sw_get_flux_data

   !---------------------------------------------------------------------------
   !>
   !   Clean up allocated memory

   subroutine sw_cleanup(me)
      class(sw_data_type), intent(inout) :: me
      if (allocated(me%mjd)) then
         deallocate(me%mjd)
         deallocate(me%f107_obs)
         deallocate(me%f107_adj)
         deallocate(me%f107a_obs_ctr)
         deallocate(me%f107a_adj_ctr)
         deallocate(me%kp)
         deallocate(me%ap_avg)
      end if
      me%initialized = .false.
      me%n_records = 0
   end subroutine sw_cleanup

   !---------------------------------------------------------------------------
   !>
   !   Copy a single record from the space weather data

   subroutine copy_record(me, idx, flux_data)
      class(sw_data_type), intent(inout) :: me
      integer(ip), intent(in) :: idx !! Index of the record to copy
      type(flux_data_type), intent(out) :: flux_data !! Output flux data structure

      flux_data%mjd           = me%mjd(idx)
      flux_data%f107_obs      = me%f107_obs(idx)
      flux_data%f107_adj      = me%f107_adj(idx)
      flux_data%f107a_obs_ctr = me%f107a_obs_ctr(idx)
      flux_data%f107a_adj_ctr = me%f107a_adj_ctr(idx)
      flux_data%kp            = me%kp(:, idx)
      flux_data%ap_avg        = me%ap_avg(idx)
   end subroutine copy_record

end module space_weather_module
