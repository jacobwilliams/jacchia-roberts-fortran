!------------------------------------------------------------------------------
! Program: test_density_profiles
!------------------------------------------------------------------------------
! Description:
!   Generate density vs altitude profiles for Solar Cycle 22 epochs and plot
!   them using pyplot-fortran. Each line represents a different epoch (solar
!   minimum, maximum, and moderate activity), showing how atmospheric density
!   varies with altitude and solar activity.
!
!   Solar Cycle 22: September 1986 - May 1996 (peak July 1989)
!
! Output:
!   - Console output with progress information
!   - Plot file: density_profiles.png
!------------------------------------------------------------------------------

program test_density_profiles
   use jacchia_roberts_module, only: jacchia_roberts_type
   use jacchia_roberts_kinds, only: ip, dp
   use pyplot_module, only: pyplot
   implicit none

   ! Constants
   real(dp), parameter :: PI = acos(-1.0_dp)
   real(dp), parameter :: RAD_PER_DEG = PI / 180.0_dp
   real(dp), parameter :: EARTH_POLAR_RADIUS = 6356.766_dp
   real(dp), parameter :: EARTH_RADIUS = EARTH_POLAR_RADIUS

   ! Variables
   type(jacchia_roberts_type) :: atm
   type(pyplot) :: plt
   integer(ip) :: sw_status, i, j
   real(dp) :: density
   real(dp) :: height, r_mag
   real(dp) :: position(3), sun_vector(3), geo_lat

   ! Altitude range
   integer(ip), parameter :: n_altitudes = 100
   real(dp) :: altitudes(n_altitudes)
   real(dp) :: h_min, h_max, dh

   ! Epoch range
   integer(ip), parameter :: n_epochs = 3
   real(dp) :: epochs(n_epochs)
   character(len=20) :: epoch_labels(n_epochs)
   real(dp) :: densities(n_altitudes, n_epochs)

   ! Plot styling
   character(len=10) :: colors(n_epochs)
   character(len=100) :: label_str

   ! Banner
   write(*,'(A)') ''
   write(*,'(A)') '=========================================================='
   write(*,'(A)') '  Jacchia-Roberts Density Profiles - Solar Cycle 22'
   write(*,'(A)') '=========================================================='
   write(*,'(A)') ''

   ! Initialize the module with space weather file
   call atm%initialize(EARTH_POLAR_RADIUS, filename='data/SpaceWeather-All-v1.2.txt', status=sw_status)
   if (sw_status /= 0) then
      write(*,'(A)') 'ERROR: Could not load space weather file'
      write(*,'(A)') '       Please ensure data/SpaceWeather-All-v1.2.txt exists'
      stop 1
   end if
   write(*,'(A)') 'Module initialized successfully'
   write(*,'(A)') ''

   ! Set up altitude range (90 km to 2000 km, logarithmic spacing works well)
   h_min = 90.0_dp
   h_max = 2000.0_dp
   ! Use logarithmic spacing for better resolution at low altitudes
   do i = 1, n_altitudes
      altitudes(i) = h_min * (h_max/h_min)**((i-1.0_dp)/(n_altitudes-1.0_dp))
   end do

   ! Define epochs spanning different solar activity periods
   ! Solar Cycle 22: Sept 1986 - May 1996, peak July 1989
   epochs(1) = 46796.0_dp   ! Jan 1, 1987 (solar minimum)
   epochs(2) = 47892.0_dp   ! Jan 1, 1990 (solar maximum)
   epochs(3) = 48988.0_dp   ! Jan 1, 1993 (moderate/declining)

   epoch_labels(1) = '1987 (min)'
   epoch_labels(2) = '1990 (max)'
   epoch_labels(3) = '1993 (moderate)'

   ! Define colors for each epoch
   colors(1) = 'green'
   colors(2) = 'red'
   colors(3) = 'orange'

   ! Fixed position for all calculations
   ! Use equatorial position at noon local solar time
   geo_lat = 0.0_dp  ! Equator
   sun_vector = [1.0_dp, 0.0_dp, 0.0_dp]  ! Sun on +X axis

   write(*,'(A)') 'Computing density profiles...'
   write(*,'(A)') ''

   ! Loop over epochs
   do j = 1, n_epochs
      write(*,'(A,I1,A,A,A,F10.1)') '  Epoch ', j, ': ', trim(epoch_labels(j)), &
                                     ' (MJD ', epochs(j), ')'

      ! Loop over altitudes for this epoch
      do i = 1, n_altitudes
         height = altitudes(i)
         r_mag = EARTH_RADIUS + height

         ! Position on equator, sun overhead
         position = [r_mag, 0.0_dp, 0.0_dp]

         ! Calculate density
         density = atm%density(height, position, sun_vector, geo_lat, epochs(j))
         densities(i, j) = density
      end do
   end do

   write(*,'(A)') ''
   write(*,'(A)') 'Creating plot...'

   ! Create plot using pyplot-fortran
   call plt%initialize(grid=.true., xlabel='Altitude (km)', &
                       ylabel='Density (kg/m$^3$)', &
                       title='Jacchia-Roberts Density Profiles - Solar Cycle 22', &
                       legend=.true., figsize=[20,7], &
                        font_size = 20,&
                        axes_labelsize = 20,&
                        xtick_labelsize = 20,&
                        ytick_labelsize = 20,&
                        legend_fontsize = 20,&
                       axis_equal=.false., tight_layout=.true.)

   ! Add each epoch as a separate line
   do j = 1, n_epochs
      call plt%add_plot(altitudes, densities(:,j), &
                        label=trim(epoch_labels(j)), &
                        linestyle='-', linewidth=3, yscale='log',&
                        xlim=[0.0_dp, 2000.0_dp])
   end do

   ! Save the plot
   call plt%savefig('density_profiles.png', pyfile='plot_density_profiles.py')

   write(*,'(A)') ''
   write(*,'(A)') 'Plot saved to: density_profiles.png'
   write(*,'(A)') 'Python script: plot_density_profiles.py'
   write(*,'(A)') ''

   ! Print some statistics
   write(*,'(A)') 'Density Statistics:'
   write(*,'(A)') '------------------------------------------------------'
   write(*,'(A)') '  Epoch              Alt (km)   Density (kg/m³)'
   write(*,'(A)') '------------------------------------------------------'
   do j = 1, n_epochs
      ! Show density at 400 km (typical ISS altitude)
      do i = 1, n_altitudes
         if (altitudes(i) >= 400.0_dp) then
            write(*,'(2X,A,F8.1,2X,ES12.5)') trim(epoch_labels(j)), &
                                             altitudes(i), densities(i,j)
            exit
         end if
      end do
   end do
   write(*,'(A)') '------------------------------------------------------'
   write(*,'(A)') ''

   ! Clean up
   call atm%destroy()

   write(*,'(A)') 'Test completed successfully!'
   write(*,'(A)') ''

end program test_density_profiles
