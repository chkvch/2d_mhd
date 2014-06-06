module fluid

	use omp_lib
	
	implicit none
	
	integer :: nz, nm, nx ! number of zones, number of fourier modes, number of x zones
	integer :: output_interval, photo_interval, photo_count = 0, num_levels, png_count
	integer :: ktrack, kbdy, xwin_id1, xwin_id2,xwin_id3,xwin_id4, xwin_id5,png_id, png_interval
	integer :: step, steps_to_take, initial_condition_type, plot_type, num_mag_contours
	double precision, dimension(:), pointer :: z, sub, dia, sup, work1, work2, temp_at_x, x ! nz
	double precision, dimension(:,:), pointer :: temp, omega, stream, &
		temp_old, omega_old, stream_old, temp_zz, omega_zz, stream_zz, temp_z, omega_z, stream_z, &
		mag, mag_old, mag_z, mag_zz, current(:,:), current_z(:,:) ! nz, nm
	double precision, dimension(:,:), pointer :: temp_xspace, stream_xspace, omega_xspace, &
			cos_ncx, sin_ncx, mag_xspace, temp_zz_xspace, bx_xspace,bz_xspace, &
			current_xspace,current_z_xspace,mag_z_xspace,mag_zz_xspace
	double precision, dimension(:,:,:), pointer :: tdot, omegadot, magdot ! nz, nm, 2 (t_old and t)
	double precision :: pi = 4.*atan(1.), dz, dt, time, ra, pr, a, c, dx, tau, t_bdy, dt_alfven, &
		roberts, chandra, tmax_anywhere, &
		dt_thermal_diffusion,dt_gravity_wave,dt_mag_diffusion,dt_diffusion_gravity, dt_factor
	logical :: do_pgplot_xwin, do_load_state = .false., do_save_state = .true., include_background_temperature, &
		do_pgplot_pngs, include_thermal_forcing, timestep_cfl, pure_hydro = .false.
	character(len=100) :: infile, outfile
		
	contains
		
	subroutine read_inlist
		integer :: unit = 5
		
		namelist /mesh/ nz, nm, nx, dt_factor
		namelist /params/ steps_to_take, ra, pr, a, tau, initial_condition_type, &
			include_thermal_forcing, roberts, chandra, pure_hydro
		namelist /output/ do_pgplot_xwin, output_interval, photo_interval, &
			include_background_temperature, num_levels, do_pgplot_pngs,png_interval, &
			num_mag_contours,ktrack
		
		open(unit,file='inlist',action='read',delim='quote',status='old')
		read(unit,nml=mesh)
		read(unit,nml=params)
		read(unit,nml=output)
		close(unit)
	end subroutine read_inlist
	
	subroutine initialize_stable_above_unstable(k_bdy,t_bdy)
		integer, intent(in) :: k_bdy
		double precision, intent(in) :: t_bdy
		integer :: k
		write(*,*) 'set up stable gradient above convection zone'
		write(*,*)
		! also have forcing term in linear contribution to Tdot in evolve
		temp(k_bdy,0) = t_bdy
		do k=1,k_bdy ! lower (unstable) region
			temp(k,0) = 1d0 - (1-temp(k_bdy,0))/z(k_bdy) * z(k) ! linear decrease in lower (convective) zone
! 			temp(k,n_perturb) = sin(pi*z(k)/z(k_bdy)) ! initial perturbation
		end do
		do k=k_bdy+1,nz ! upper (stable) region
			temp(k,0) = t_bdy + (1-temp(k_bdy,0))/(1-z(k_bdy)) * ( z(k) - z(k_bdy) )
		end do
	end subroutine initialize_stable_above_unstable
	
	subroutine initialize_unstable_above_stable(k_bdy,t_bdy)
		integer, intent(in) :: k_bdy
		double precision, intent(in) :: t_bdy
		integer :: k
		double precision :: zp
		write(*,*) 'set up stable gradient beneath convection zone'
		write(*,*)
		! set up stable below kbdy, unstable above
		! also have forcing term in linear contribution to Tdot in evolve
		temp(:,0) = temp(:,0) / z(k_bdy) ! temperature goes to 1 at k_bdy
		temp(nz,0) = t_bdy ! temperature at top boundary
		do k=k_bdy,nz ! upper (unstable) region
			zp = 1.0 - ( (1-z(k))/(1-z(k_bdy)) )
			temp(k,0) = (1.0-zp) + zp * temp(nz,0) ! linear decrease in upper (convective) zone
! 			temp(k,n_perturb) = 0.3 * sin(pi*zp) ! perturbation in n=1
		end do
	end subroutine initialize_unstable_above_stable
	
	subroutine initialize_unstable
		integer :: k
		do k=1,nz
			temp(k,0) = 1d0 - z(k) ! linear initial background temperature (convectively unstable)
		end do
		temp(nz,:) = 0d0
	end subroutine initialize_unstable

	subroutine initialize_stable
		integer :: k
		do k=1,nz
			temp(k,0) = z(k) ! linearly INCREASING background (convectively stable)
		end do
	end subroutine initialize_stable
	
	subroutine sine_perturbation(amplitude,n_perturb,kmin,kmax)
		double precision, intent(in) :: amplitude
		integer, intent(in) :: n_perturb, kmin, kmax
		integer :: k
		do k=kmin,kmax
			temp(k,n_perturb) = temp(k,n_perturb) + amplitude * sin(pi*(z(k)-z(kmin))/z(kmax))
		end do
	end subroutine sine_perturbation
		
			
	
	subroutine initialize
		
		integer :: k, ierr, pgopen, j, n, nthreads
		allocate(z(nz), sub(nz), dia(nz), sup(nz), work1(nz), work2(nz), x(nx), temp_at_x(nz))
		allocate(temp(nz,0:nm),omega(nz,0:nm),stream(nz,0:nm), &
			temp_z(nz,0:nm), omega_z(nz,0:nm), stream_z(nz,0:nm), &
			temp_zz(nz,0:nm), omega_zz(nz,0:nm), stream_zz(nz,0:nm), &
			temp_old(nz,0:nm), omega_old(nz,0:nm), stream_old(nz,0:nm), &
			mag(nz,0:nm), mag_z(nz,0:nm), mag_zz(nz,0:nm), current(nz,0:nm), current_z(nz,0:nm))
		allocate(omegadot(nz,0:nm,2),tdot(nz,0:nm,2),magdot(nz,0:nm,2))
		allocate(temp_xspace(nx,nz),stream_xspace(nx,nz),omega_xspace(nx,nz),mag_xspace(nx,nz), &
			temp_zz_xspace(nx,nz),bx_xspace(nx,nz),bz_xspace(nx,nz),&
			current_xspace(nx,nz),current_z_xspace(nx,nz),mag_zz_xspace(nx,nz),mag_z_xspace(nx,nz))
		allocate(cos_ncx(nm,nx),sin_ncx(nm,nx))
	
 	nthreads=OMP_GET_MAX_THREADS()
 	write(*,*)
 	write(*,*) 'nthreads ', nthreads
 	write(*,*)
	call omp_set_num_threads(nthreads)
						
		if (do_load_state) then
			call load_state
		else ! get initial state from inlist
			step = 0
			time = 0d0
			
			dz = 1. / (nz-1)
			dx = a / (nx-1)

			do j=1,nx ! initialize x grid (will only transform to do plots)
				x(j) = dx * (j-1)
			end do
			
			! zero fluid velocity at t=0
	 		stream(:,:) = 0d0
	 		omega(:,:) = 0d0
			
			! purely vertical field at t=0 (second paragraph p. 185 of Glatz)
			mag(:,:) = 0d0
					
			! initialize z grid
			temp(:,:) = 0d0
			do k=1,nz
				z(k) = dz * (k-1)
			end do
			
			if (initial_condition_type.eq.0) then
				call initialize_stable_above_unstable(nz/2,1d-1) ! k_bdy, temp_bdy
				include_thermal_forcing = .true.
			else if (initial_condition_type.eq.1) then
				call initialize_unstable_above_stable(nz * 3 / 4 + 1,3d-1) ! k_bdy, temp_bdy
				include_thermal_forcing = .true.
			else
				call initialize_unstable
				include_thermal_forcing = .false.
			end if

			call sine_perturbation(1d-2,1,1,nz) ! amplitude, mode number, k_min, k_max
			call sine_perturbation(1d-2,8,1,nz) ! amplitude, mode number, k_min, k_max
			
			! initialize a spike in vorticity at midpoint in z
! 			do n=2,nm,2
! 				omega(50,n) = -1. * (-1.0)**(n/2) * n / nm
! 				write(*,*) 'n, omega(50,n) ', n, omega(50,n)
! 			end do

			tdot(:,:,1) = 0d0 ! initial "old" tdot
			tdot(:,:,2) = 0d0 ! initial "new" tdot -- will calculate in evolve loop
			omegadot(:,:,1) = 0d0 ! initial "old" omegadot
			omegadot(:,:,2) = 0d0 ! initial "new" omegadot
			magdot(:,:,1) = 0d0
			magdot(:,:,2) = 0d0

			png_count = 0
			
		end if ! load state or set initial conditions

		c = pi / a	
		dz = 1. / (nz-1)
		
		! spatial grid in z
		do k=1, nz
			z(k) = dz * (k-1)
		end do

		! timestep
		dt_thermal_diffusion = dz**2/4 ! thermal diffusion time over a cell
		dt_gravity_wave = 2*pi/(50*sqrt(ra*pr)) ! gravity wave time
		dt_mag_diffusion = dz**2*roberts/4 ! magnetic diffusion time
		dt_diffusion_gravity = minval( (/ dt_thermal_diffusion, dt_gravity_wave, dt_mag_diffusion  /) )

		dt = dt_diffusion_gravity * dt_factor
		timestep_cfl = .false.

		write(*,*) 'initial dt = ', dt
		
		work1(:) = 0d0
		work2(:) = 0d0
		sub(:) = -1./dz**2
		sup(:) = -1./dz**2
		! streamfunction = 0 at k=1 and k=nz
		sub(nz) = 0d0
		sup(1) = 0d0

		! initialize derivatives, etc -- to be calculated initially in evolve anyway
		temp_z(:,:) = 0d0
		omega_z(:,:) = 0d0
		stream_z(:,:) = 0d0
		temp_zz(:,:) = 0d0
		omega_zz(:,:) = 0d0
		stream_zz(:,:) = 0d0
		temp_old(:,:) = 0d0
		omega_old(:,:) = 0d0
		stream_old(:,:) = 0d0
		mag_z(:,:) = 0d0
		mag_zz(:,:) = 0d0

		! compute sines and cosines on x grid to speed up transform later
		do n=1,nm
			do j=1,nx
				cos_ncx(n,j) = cos(n*c*x(j))
				sin_ncx(n,j) = sin(n*c*x(j))
			end do
		end do
		
		if (do_pgplot_xwin) then
			xwin_id1 = pgopen('/xwin')
	        if (xwin_id1 .le. 0) stop 'pgopen failed to initialize device /xwin'
			
			call init_colors

			if(.not.pure_hydro) then
				xwin_id2 = pgopen('/xwin')
				xwin_id3 = pgopen('/xwin')
! 				xwin_id4 = pgopen('/xwin')
! 				xwin_id5 = pgopen('/xwin')
			else if(pure_hydro) then
				xwin_id2=pgopen('/xwin')
			end if
			
			call init_colors
		end if
		
	end subroutine initialize
	
	subroutine cleanup
		deallocate(z, sub, dia, sup, work1, work2, stream, omega, temp, &
		temp_old, omega_old, stream_old, omegadot, tdot, temp_zz, omega_zz, stream_zz, &
		omega_z, temp_z, stream_z, temp_at_x, temp_xspace, omega_xspace, x, cos_ncx, sin_ncx,&
		current,mag,mag_z,mag_zz,mag_xspace,temp_zz_xspace,bx_xspace,bz_xspace,current_xspace,&
		current_z_xspace,mag_zz_xspace,mag_z_xspace)
		if (do_pgplot_xwin) call pgend()
	end subroutine cleanup
	
	subroutine tridiag(nz, rhs, sol, sub, dia, sup, work1, work2)
		! tridiagonal matrix solver
		double precision, dimension(:) :: rhs, sol, sub, dia, sup, work1, work2
		integer, intent(in) :: nz
		integer :: i
				
		work1(1) = 1. / dia(1)
		work2(1) = sup(1) * work1(1)
		do i=2, nz-1
			work1(i) = 1. / (dia(i) - sub(i) * work2(i-1))
			work2(i) = sup(i) * work1(i)
		end do
		work1(nz) = 1. / (dia(nz)-sub(nz) * work2(nz-1))
		
		sol(1) = rhs(1) * work1(1)
		do i=2, nz
			sol(i) = (rhs(i) - sub(i) * sol(i-1)) * work1(i)
		end do
		do i=nz-1, 1, -1
			sol(i) = sol(i) - work2(i) * sol(i+1)
		end do
		
		return
	end subroutine tridiag

	subroutine report
		implicit none 
		integer :: n
		
! 		write(*,*) 'max temperature anywhere', tmax_anywhere
		write(*,'(8a16)') 'step', 'time', 'ktrack', 'nz', 'nm', 'nx', 'dt', 'cfl timestep?'
		write(*,'(i16,f16.8,4i16,es16.4,l16)') step, time, ktrack, nz, nm, nx, dt, timestep_cfl
		
		write(*,'(a3,7a16)') 'n', 'temperature', 'omega', 'stream', 'mag', &
			'magdot_now', 'magdot_old', 'current'
		do n=0,5
			write(*,'(i3,7es16.6)') n, &
				temp(ktrack,n), omega(ktrack,n), stream(ktrack,n), mag(ktrack,n), &
				magdot(ktrack,n,2), magdot(ktrack,n,1), current(ktrack,n)
		end do
		write(*,*) '...snipped...'
		do n=nm-5,nm
			write(*,'(i3,7es16.6)') n, &
				temp(ktrack,n), omega(ktrack,n), stream(ktrack,n), mag(ktrack,n), &
				magdot(ktrack,n,2), magdot(ktrack,n,1), current(ktrack,n)
		end do
! 		write(*,*)
! 		write(*,*) 'max temp @ z=1 ', maxval(temp(nz,:))
		write(*,*)
	end subroutine report
			
	subroutine load_state
		integer :: k, n, j
		integer :: nil
		
		open(15,file=infile)
		read(15,*) ! skip header
		read(15,*) ! skip header
		read(15,'(4i20,7es20.8,i20)') nz, nm, nx, step, time, ra, pr, &
			roberts, chandra, a, dt, photo_count
		read(15,*)
		read(15,*)
		do k=1,nz
			do n=0,nm
				read(15,'(2i5, 7es24.10)') nil, nil, temp(k,n), omega(k,n), stream(k,n), &
					mag(k,n), tdot(k,n,1), omegadot(k,n,1), magdot(k,n,1)
			end do
		end do
		
		close(15)

		! reconstruct x grid. could just invent a new one
		dx = a / (nx-1)
		do j=1,nx
			x(j) = dx * (j-1)
		end do
		
		write(*,*)
		write(*,*) 'loaded initial state from ', infile
		write(*,*)
		photo_count = photo_count + 1
	end subroutine load_state
	
	subroutine save_state
		integer :: k, n
				
		if (photo_count.lt.10)	then
			write(outfile,'(a10,i1)') 'photos/000', photo_count
		else if (photo_count.lt.100) then
			write(outfile,'(a9,i2)') 'photos/00', photo_count
		else if (photo_count.lt.1000)	then
			write(outfile,'(a8,i3)') 'photos/0', photo_count
		else
			write(outfile,'(a7,i4)') 'photos/', photo_count
		end if
		
		open(18,file='photos/photos.index',action='write',access='append')
!		write(18,'(2a16)') 'step', 'photo number'
		write(18,'(2i12)') step, photo_count
		close(18)

		open(16,file=trim(outfile),action='write',status='replace')
		write(16,*) '# saved model from flow. first lines are global info, rest is hydro variables for all k and n.'
		write(16,'(13a20)') 'nz', 'nm', 'nx', 'step', 'time', 'Ra', 'Pr', &
			'q', 'Q', 'a', 'dt', 'photo_count', 'png_count'
		write(16,'(4i20,7es20.8,2i20)') nz, nm, nx, step, time, ra, pr, &
			roberts, chandra, a, dt, photo_count, png_count
		write(16,*)
		write(16,'(2a5, 7a24)') 'k', 'n', 'temp', 'omega', 'stream', 'mag', &
			'tdot_old', 'omegadot_old', 'magdot_old'
		do k=1,nz
			do n=0,nm
				write(16,'(2i5, 7es24.10E3)') k, n, temp(k,n), omega(k,n), stream(k,n), &
					mag(k,n), tdot(k,n,1), omegadot(k,n,1), magdot(k,n,1)
			end do
		end do
		close(16)

		write(*,*) 'wrote ', outfile
		photo_count = photo_count + 1

	end subroutine save_state
	
	
	
	
	
	
	
	subroutine compute_interior_spatial_derivatives(k)
		integer, intent(in) :: k
		integer :: n
		! central finite differences
!$OMP PARALLEL DO PRIVATE(n)
		do n=0,nm
			temp_zz(k,n) = (temp(k+1,n) - 2 * temp(k,n) + temp(k-1,n)) / dz**2
			omega_zz(k,n) = (omega(k+1,n) - 2 * omega(k,n) + omega(k-1,n)) / dz**2
			stream_zz(k,n) = (stream(k+1,n) - 2 * stream(k,n) + stream(k-1,n)) / dz**2				
		
			temp_z(k,n) = (temp(k+1,n)-temp(k-1,n)) / 2 / dz
			omega_z(k,n) = (omega(k+1,n)-omega(k-1,n)) / 2 / dz
			stream_z(k,n) = (stream(k+1,n)-stream(k-1,n)) / 2 / dz
			
			mag_z(k,n) = (mag(k+1,n)-mag(k-1,n)) / 2 / dz
			mag_zz(k,n) = (mag(k+1,n) - 2 * mag(k,n) + mag(k-1,n)) / dz**2
		end do
!$OMP END PARALLEL DO
	end subroutine compute_interior_spatial_derivatives

	subroutine compute_boundary_spatial_derivatives(n)
		integer, intent(in) :: n
		! centered finite differences, using ghost points (eqs. 11.37-40)
		mag_zz(1,n) = 2d0 * (mag(2,n) - mag(1,n)) / dz**2
		mag_zz(nz,n) = 2d0 * (mag(nz-1,n) - mag(nz,n)) / dz**2
		stream_z(1,n) = stream(2,n) / dz
		stream_z(nz,n) = -1d0 * stream(nz-1,n) / dz
		
		! boundary conditions (eq. 11.24)
		mag_z(1,n) = 0d0
		mag_z(nz,n) = 0d0
	end subroutine compute_boundary_spatial_derivatives

	subroutine compute_current(k,n)
		integer, intent(in) :: k, n
		! eq. 11.33 for J_n. J(x,z,t) is expanded in sines in x
		current(k,n) = -1d0 * (mag_zz(k,n) - mag(k,n) * (n * c)**2d0)
	end subroutine compute_current
	
	subroutine compute_current_z(k,n)
		integer, intent(in) :: k,n
		current_z(k,n) = (current(k+1,n) - current(k-1,n)) / 2d0 / dz
	end subroutine compute_current_z
			
	subroutine add_thermal_forcing(k) ! only works for stable above unstable, as written
		integer, intent(in) :: k
		double precision :: term
		if (k.gt.kbdy) then
			term = -1d0 * (temp(k,0) - t_bdy - (1-t_bdy)/(1-z(kbdy))*(z(k)-z(kbdy))) / tau
			tdot(k,0,2) = tdot(k,0,2) + term  ! thermal forcing term
		end if
		write(*,*) 'doing thermal forcing at zone', k
	end subroutine add_thermal_forcing
	
	subroutine add_linear_terms(k,n)
		integer, intent(in) :: k, n
		! linear parts of the time derivatives (eqs. 3.3 and 3.4)
		tdot(k,n,2) = temp_zz(k,n) - (n * c)**2*temp(k,n)
		! tdot(k,n,2) = tdot(k,n,2) + c * pi * stream(k,n) ! linearized advection term, needed for linear case
		omegadot(k,n,2) = ra * pr * n * c * temp(k,n) &
			+ pr * (omega_zz(k,n) - (n * c)**2 * omega(k,n)) ! velocity diffusion term
	end subroutine add_linear_terms
	
	subroutine add_linear_magnetic_terms(k,n)
		integer, intent(in) :: k, n
		! eq. 11.25 for dA_n/dt
		! diffusion of the magnetic vector potential
		! eta -> 1/roberts
		magdot(k,n,2) = magdot(k,n,2) + 1d0 / roberts * (mag_zz(k,n) - (n * c)**2d0 * mag(k,n))
		
		! linear part of the advection of the vector potential  (second paragraph page 186)
		if(.not.pure_hydro) magdot(k,n,2) = magdot(k,n,2) + stream_z(k,n)
		
		! linear part of the Lorentz torque (first paragraph page 186 in Glatzmaier)
		if((k.gt.1).and.(k.lt.nz)) then ! otherwise no such thing as current_z
			omegadot(k,n,2) = omegadot(k,n,2) + chandra * pr / roberts * current_z(k,n)
		end if
	end subroutine add_linear_magnetic_terms
	
	subroutine galerkin_advection(k,n,factor,timederiv,advector,advectee) 
		! will not include an np=0 term (give Tdot one of those externally)
		double precision, intent(inout) :: timederiv(nz,nm)
		double precision, intent(in) :: advector(nz,nm), advectee(nz,nm)
		integer, intent(in) :: k,n
		double precision, intent(in) :: factor
		integer :: np ! n prime
		
		! eq. 4.4
		
		return
		
	end subroutine
			
	subroutine add_advection_terms(k)
		integer, intent(in) :: k
		integer :: n, np
		! update tdot for n=0, npp=np
		do np=1,nm
			tdot(k,0,2) = tdot(k,0,2) - c / 2 * np * &
				(stream_z(k,np) * temp(k,np) + stream(k,np) * temp_z(k,np))
		end do
!$OMP PARALLEL DO PRIVATE(n,np)	      
		do n=1,nm
			! do the np=0 case separately (only Tdot gets a contribution, for npp=n)
			tdot(k,n,2) = tdot(k,n,2) - n * c * stream(k,n) * temp_z(k,0)
			!       update tdot and omegadot for all npp, now n>0 and np>0
			do np=1,nm	! three values of npp contribute: n-np, n+np, np-n.
	
			! npp = n - np
			   if((n-np.ge.1).and.(n-np.le.nm)) tdot(k,n,2) = tdot(k,n,2) - c / 2 * ( &
			   -1d0 * np * stream_z(k,n-np) * temp(k,np) + (n-np) * stream(k,n-np) * temp_z(k,np))
			! npp = n + np
			   if((n+np.ge.1).and.(n+np.le.nm)) tdot(k,n,2) = tdot(k,n,2) - c / 2 * ( &
			   np * stream_z(k,n+np) * temp(k,np) + (n+np) * stream(k,n+np) * temp_z(k,np))
			! npp = np - n
			   if((np-n.ge.1).and.(np-n.le.nm)) tdot(k,n,2) = tdot(k,n,2) - c / 2 * ( &
			   np * stream_z(k,np-n) * temp(k,np) + (np-n) * stream(k,np-n) * temp_z(k,np))
		   
			! npp = n - np						
			   if((n-np.ge.1).and.(n-np.le.nm)) then 
				   omegadot(k,n,2) = omegadot(k,n,2) - c / 2 * & ! advection of vorticity
				   (-1d0 * np * stream_z(k,n-np) * omega(k,np) + (n-np) * stream(k,n-np) * omega_z(k,np)) &
				   + c / 2 * chandra * pr / roberts * & ! Lorentz torque
				   (-1d0 * np * mag_z(k,n-np) * current(k,np) + (n-np) * mag(k,n-np) * current_z(k,np)) 
			   end if
			! npp = n + np
			   if((n+np.ge.1).and.(n+np.le.nm)) then
				   omegadot(k,n,2) = omegadot(k,n,2) - c / 2 * & ! advection of vorticity
				   (-1d0 * (np * stream_z(k,n+np) * omega(k,np) + (n+np) * stream(k,n+np) * omega_z(k,np))) &
				   + c / 2 * chandra * pr / roberts * & ! Lorentz torque
				   (-1d0 * (np * mag_z(k,n+np) * current(k,np) + (n+np) * mag(k,n+np) * current_z(k,np))) 
			   end if
			! npp = np - n
			   if((np-n.ge.1).and.(np-n.le.nm)) then
				   omegadot(k,n,2) = omegadot(k,n,2) - c / 2 * & ! advection of vorticity
				   (np * stream_z(k,np-n) * omega(k,np) + (np-n) * stream(k,np-n) * omega_z(k,np)) &
				   + c / 2 * chandra * pr / roberts * & ! Lorentz torque
				   (np * mag_z(k,np-n) * current(k,np) + (np-n) * mag(k,np-n) * current_z(k,np))
			   end if
		    
		   end do		! loop over np
		end do		! loop over n
!$OMP END PARALLEL DO
		
	end subroutine add_advection_terms		
		
	subroutine add_mag_advection_terms(k,n)
		integer, intent(in) :: k,n
		integer :: np
				
		! following eq. 4.4 but with omega(k,np) -> mag(k,np) for n > 0
	      
		!       update magdot for all npp, now n>0 and np>0
!$OMP PARALLEL DO PRIVATE(np)
		do np=1,nm	! three values of npp contribute: n-np, n+np, np-n.
		   
		! npp = n - np						
		   if((n-np.ge.1).and.(n-np.le.nm)) then 
			   magdot(k,n,2) = magdot(k,n,2) &
			   - c / 2 * & ! advection of induced vector potential
			   (-1d0 * np * stream_z(k,n-np) * mag(k,np) + (n-np) * stream(k,n-np) * mag_z(k,np))		  
		   end if
		! npp = n + np
		   if((n+np.ge.1).and.(n+np.le.nm)) then
			   magdot(k,n,2) = magdot(k,n,2) &
			   - c / 2 * & ! advection of induced vector potential
			   (-1d0 * (np * stream_z(k,n+np) * mag(k,np) + (n+np) * stream(k,n+np) * mag_z(k,np)))
		   end if
		! npp = np - n
		   if((np-n.ge.1).and.(np-n.le.nm)) then
			   magdot(k,n,2) = magdot(k,n,2) &
			   - c / 2 * & ! advection of induced vector potential
			   (np * stream_z(k,np-n) * mag(k,np) + (np-n) * stream(k,np-n) * mag_z(k,np))
		   end if
	    
	   end do		! loop over np
!$OMP END PARALLEL DO
		
	end subroutine add_mag_advection_terms		
		
	subroutine transform_to_xspace
		integer :: k, j, n
		temp_xspace(:,:) = 0d0 ! nx, nz
		mag_xspace(:,:) = 0d0 ! nx, nz
		temp_zz_xspace(:,:) = 0d0 ! nx, nz
		tmax_anywhere = 0d0
		current_xspace(:,:) = 0d0
		current_z_xspace(:,:) = 0d0
		mag_zz_xspace(:,:) = 0d0
		mag_z_xspace(:,:) = 0d0
!$OMP PARALLEL DO PRIVATE(k,j,n)
		do k=1,nz
			do j=1,nx
				if (include_background_temperature) then
					temp_xspace(j,k) = temp_xspace(j,k) + temp(k,0)
				end if
				mag_xspace(j,k) = mag_xspace(j,k) + x(j) ! n=0 contribution B_0*x
				temp_zz_xspace(j,k) = temp_zz_xspace(j,k) + temp_zz(k,0)
				do n=1,nm
					current_xspace(j,k) = current_xspace(j,k) + current(k,n) * sin_ncx(n,j)
					current_z_xspace(j,k) = current_z_xspace(j,k) + current_z(k,n) * sin_ncx(n,j)
					temp_xspace(j,k) = temp_xspace(j,k) + temp(k,n) * cos_ncx(n,j)
 					temp_zz_xspace(j,k) = temp_zz_xspace(j,k) + temp_zz(k,n) * cos_ncx(n,j)
! 					omega_xspace(j,k) = omega_xspace(j,k) + omega(k,n) * sin_ncx(n,j)
					mag_xspace(j,k) = mag_xspace(j,k) + mag(k,n) * sin_ncx(n,j)
					mag_z_xspace(j,k) = mag_z_xspace(j,k) + mag_z(k,n) * sin_ncx(n,j)
					mag_zz_xspace(j,k) = mag_zz_xspace(j,k) + mag_zz(k,n) * sin_ncx(n,j)
	!	stream_xspace(j,k) = stream_xspace(j,k) + stream(k,n) * sin_ncx(n,j)
				end do
				tmax_anywhere = maxval( (/ tmax_anywhere, temp_xspace(j,k) /) )
			end do
		end do
!$OMP END PARALLEL DO
	end subroutine transform_to_xspace		
	
	subroutine set_dt_magnetic_cfl_condition
		! see eqs. 11.30 and 11.31 in Glatz
		integer :: k,j,n	
		double precision :: alfven_speed_squared, alfven_speed, max_alfven_speed
		bx_xspace(:,:) = 0d0 ! nx, nz 
		bz_xspace(:,:) = 1d0 ! nx, nz ! including background term here
		max_alfven_speed = 0d0

!$OMP PARALLEL DO PRIVATE(k,j,n)
		do k=1,nz
			do j=1,nx
				do n=1,nm
					bx_xspace(j,k) = bx_xspace(j,k) - mag_z(k,n) * sin_ncx(n,j)
					bz_xspace(j,k) = bz_xspace(j,k) + n * c * mag(k,n) * cos_ncx(n,j)
				end do
				alfven_speed_squared = 1d0 * chandra * pr / roberts * &
					(bx_xspace(j,k)**2 + bz_xspace(j,k)**2)
				alfven_speed = sqrt(alfven_speed_squared)
				if(alfven_speed.gt.max_alfven_speed) max_alfven_speed = alfven_speed
!				write(*,*) 'v2, v, maxv', alfven_speed_squared, alfven_speed, max_alfven_speed
			end do
		end do
!$OMP END PARALLEL DO
		! CFL condition
		dt_alfven = dz / max_alfven_speed
		if(mod(step,output_interval).eq.0) write(*,*) 'max alfven speed = ', max_alfven_speed, 'dt_alfven = ', dt_alfven
		if (dt_diffusion_gravity.gt.dt_alfven) then
			timestep_cfl = .true.
			dt = dt_alfven * dt_factor
		else
			timestep_cfl = .false.
			dt = dt_diffusion_gravity * dt_factor
		end if
	end subroutine set_dt_magnetic_cfl_condition		
		
		
		
		
		
		
		
		
		
		
		
		
	subroutine evolve
	integer :: k, n, np, j, ierr, perr
	
	do while (step .le. steps_to_take) ! step through time
	   
	! update tdot and omegadot (Eq. 3.3, 3.4, and 2.16)
	do k=2,nz-1

	! get temp_z, temp_zz, mag_z, mag_zz, omega_z, omega_zz, stream_z, stream_zz
	! centered finite differences
	call compute_interior_spatial_derivatives(k) ! does n=0,nm
							   
!$OMP PARALLEL DO PRIVATE(n,perr)		   
	do n=0,nm
		 
!$OMP CRITICAL
		if (isnan(temp(k,n))) then
			call report
			call save_one_png(temp_xspace,0,'T')
			stop 'got NaN in temperature'
		end if
!$OMP END CRITICAL
				
		call add_linear_terms(k,n) ! initializes tdot from temp and temp_zz, and omegadot from temp, omega, and omega_zz		

		if(.not.pure_hydro) then
			if(step.gt.1) then
				call compute_current(k,n) ! minus laplacian of mag ! needs mag and mag_zz
				call compute_current_z(k,n) ! centered difference - only works for an internal zone
			end if
			! otherwise leave J and J_z as zero. magdot just gets stream_z contribution

! 			if((k.eq.ktrack).and.(n.eq.1))write(*,*) &
! 				'current(ktrack,1), current_z(ktrack,1) = ', &
! 				current(ktrack,1), current_z(ktrack,1)

			! magnetic diffusion, and linear contributions to advection of mag and the lorentz torque
			call add_linear_magnetic_terms(k,n) ! needs mag, mag_zz, stream_z, current_z
			
! 			if((k.eq.ktrack).and.(n.eq.1))write(*,*) &
! 				'after add_linear_magnetic_terms on internal zones, magdot(ktrack,1,new), magdot(ktrack,1,old) = ', & 
! 				magdot(ktrack,1,2), magdot(ktrack,1,1)
! 				if((k.eq.ktrack).and.(n.eq.1))write(*,*)

			call add_mag_advection_terms(k,n) ! -[(v dot del)A]_n ! needs mag, mag_z, stream, stream_z. see eq. 11.25, cf. eq. 4.4
				
		end if

	end do ! loop over mode number
!$OMP END PARALLEL DO

	if (include_thermal_forcing) call add_thermal_forcing(k)

	! update tdot (adv. of temperature) and omegadot (adv. of omega and induced mag)
	! needs current, current_z, mag, mag_z for Lorentz torque
	call add_advection_terms(k) 

	call galerkin_advection(1,1,1d0,omegadot(:,:,2),stream,omega) 
			  
	end do ! loop over internal zones

	if(.not.pure_hydro) then
!$OMP PARALLEL DO PRIVATE(n)
	! deal with k=1 and k=nz cases with special care: one-side finite differences using ghost points
	do n=1,nm
		call compute_boundary_spatial_derivatives(n) ! mag_z, mag_zz, and stream_z at k=1 and k=nz
	end do
!$OMP END PARALLEL DO	
!$OMP PARALLEL DO PRIVATE(n)
	do n=1,nm
		call compute_current(1,n)
		call compute_current(nz,n)
! 		call compute_current_z(1,n)
! 		call compute_current_z(nz,n)
		call add_linear_magnetic_terms(1,n)
		call add_linear_magnetic_terms(nz,n)
		call add_mag_advection_terms(1,n)
		call add_mag_advection_terms(nz,n)
	end do
!$OMP END PARALLEL DO	
	end if	   		   
			   
	temp_old = temp
	omega_old = omega
	stream_old = stream				

	! some terminal and graphical output
	if (mod(step,output_interval).eq.0) then
	
		call transform_to_xspace
		call report
	
		if(do_pgplot_xwin.or.do_pgplot_pngs) then

			if(do_pgplot_xwin) then
				call pgslct(xwin_id1)
				call plot(temp_xspace,0,'T')
			end if
			if(do_pgplot_pngs.and.(mod(step,png_interval).eq.0)) &
				call save_one_png(temp_xspace,0,'T')
			
			if(.not.pure_hydro) then ! also plot mag
				if(do_pgplot_xwin) then
					call pgslct(xwin_id2)
					call plot(mag_xspace,1,'A')
				end if
				if(do_pgplot_pngs.and.(mod(step,png_interval).eq.0)) &
					call save_one_png(temp_xspace,1,'A')

					if(do_pgplot_xwin) then
						call pgslct(xwin_id3)
						call plot(current_xspace,0,'J')
					end if
					if(do_pgplot_pngs.and.(mod(step,png_interval).eq.0)) &
						call save_one_png(current_xspace,0,'J')
					
			end if
			if(do_pgplot_pngs.and.(mod(step,png_interval).eq.0)) png_count = png_count + 1
		end if
		
	end if	! show terminal and graphical output

	! now with time derivs in hand, update temp and omega (Eq. 2.18)
	! and mag
	do k=1,nz
		do n=0,nm
			temp(k,n) = temp(k,n) + dt * 0.5 * (3 * tdot(k,n,2) - tdot(k,n,1))
			omega(k,n) = omega(k,n) + dt  * 0.5 * (3 * omegadot(k,n,2) - omegadot(k,n,1))
			if((.not.pure_hydro).and.(step.ge.1)) &
				mag(k,n) = mag(k,n) + dt * 0.5 * (3 * magdot(k,n,2) - magdot(k,n,1))
		end do
	end do


	! boundary conditions
	if (initial_condition_type.eq.0) then
		temp(1,0) = 1d0
		temp(nz,0) = 1d0
	else
		temp(nz,0) = 0d0
	end if
	temp(nz,1:nm) = 0d0 ! perturbations must vanish on boundaries
	   
	! contruct dia (a function of n only)
	! call tridiag with omega as rhs to get stream=sol
	! laplacian(psi_n) + omega_n = 0
	dia(1) = 1
	dia(nz) = 1

	do n=1,nm
		dia(2:nz-1) = (c * n)**2 + 2. / dz**2
		call tridiag(nz,omega(:,n),stream(:,n),sub,dia,sup,work1,work2)
	end do

	   
	! advance time
	tdot(:,:,1) = tdot(:,:,2)
	tdot(:,:,2) = 0d0
	omegadot(:,:,1) = omegadot(:,:,2)
	omegadot(:,:,2) = 0d0
	if((.not.pure_hydro).and.(step.ge.1)) then
		magdot(:,:,1) = magdot(:,:,2)
		magdot(:,:,2) = 0d0
	end if

	time = time + dt
	step = step + 1

	if((.not.pure_hydro).and.(mod(step,1).eq.0))call set_dt_magnetic_cfl_condition ! calculates Bx and Bz, maximum Alfven speed anywhere, and sets new dt
	   
	if (mod(step,photo_interval).eq.0) then
		call save_state
	end if
	   	   
	end do
	
!       call report_final
	write(*,*) 'hit maximum step ', steps_to_take
	  
	end subroutine evolve
		
		
		
		
		
		
		
		
		
	subroutine plot_line
        integer :: i
		character(100) :: title
		real :: ymin, ymax, xvar(nz), y1(nz), y2(nz), y3(nz), yvar(nz)
		
		write(title,('(a5,i8,a10,es16.8,a6,es16.8)')) 'step ', step, 'time ', time, 'Ra', ra
		
		xvar = real(z)
		
		yvar = real(temp_xspace(50,:))
		
		ymax = 1.05 * maxval(yvar)
		ymin = 1.05 * minval(yvar)


		call pgbbuf
		call pgsci(1)
		CALL PGENV(minval(xvar),maxval(xvar),ymin,ymax,0,1)
		
	
		call pgask(.false.) ! whether or not pgplot pauses and prompts

        CALL pglab('z', 'temperature',title)
		call pgsci(11) ! blue
		call pgpt(nz,xvar,yvar,9) ! number of points, x points, y points, marker style				
		
		call pgslw(4) ! ~default thickness
		call pgebuf
	end subroutine plot_line
	
	subroutine init_colors
		integer :: i
		real :: r,g,b, norm
		do i=1,num_levels
			r = ((1.0*i-1)/(num_levels-1)*2) - 1
			b = (1-(1.0*i-1)/(num_levels-1)*2)
			if (i.le.num_levels/2) then
				r=0
			else if (i.ge.num_levels/2) then
				b=0
			end if
			g = 1 - r - b !sin(pi*(i-1)/num_levels)**2
! 			write(*,*)i, r, g, b
			call pgscr(15+i,r,g,b)
		end do
	end subroutine init_colors
	
	subroutine save_one_png(array_to_plot,plot_type,plot_name)
		integer, intent(in) :: plot_type ! 1 for mag (just draw contours), otherwise full colormap
		double precision, intent(in) :: array_to_plot(nx,nz)
		character, intent(in) :: plot_name
		integer :: pgopen
		character(len=256) :: png_outfile

		if (png_count.lt.10)	then
			write(png_outfile,'(a10,i1)') '000', png_count
		else if (png_count.lt.100) then
			write(png_outfile,'(a9,i2)') '00', png_count
		else if (png_count.lt.1000)	then
			write(png_outfile,'(a8,i3)') '0', png_count
		else
			write(png_outfile,'(i4)') png_count
		end if
		
		png_outfile = adjustl(png_outfile)
		png_outfile = 'png/' // trim(plot_name) // '_' // trim(png_outfile)

		png_id = pgopen(trim(png_outfile) // '.png/png')
        if (png_id .le. 0) then
			write(*,*) 'pgopen /png failed', png_id
			stop
		end if

 		call pgslct(png_id)
			
		call plot(array_to_plot,plot_type,plot_name)
	    call pgclos ! close this device instance
		
	    png_id = 0
		write(*,*) 'wrote ' // trim(png_outfile) // '.png'
		
	    if (do_pgplot_xwin) call pgslct(xwin_id1) ! reselect first xwindow
	end subroutine save_one_png

	subroutine plot(array_to_plot,plot_type,plot_name)
		integer, intent(in) :: plot_type ! 1 for mag (just draw contours), otherwise full colormap
		double precision, intent(in) :: array_to_plot(nx,nz)
		character, intent(in) :: plot_name
        integer :: i, j, k, ierr
		character(100) :: title, annot
		real :: xvar(nx), xmin, xmax, yvar(nz), ymin, ymax, zvar(nz,nx), zmin, zmax, alfa
		real :: transform(6)
		real, pointer :: z_contours(:) 
		logical :: dbg = .false.
								
		xvar = real(x)
		yvar = real(z)
		zvar = real(array_to_plot)
				
		call init_colors
		
		zmin = 0d0
		zmax = 0d0
		do k=1,nz
			do j=1,nx
				if (zvar(k,j).gt.zmax) zmax = zvar(k,j)
				if (zvar(k,j).lt.zmin) zmin = zvar(k,j)
			end do
		end do
			
		if (dbg) write(*,*) 'zmin, zmax', zmin, zmax
		
		write(annot,('(a5,i8,a7,es12.4,a7,es12.4,a7,es12.4)')) 'step ', step, 'time ', time, 'zmin', zmin, 'zmax ', zmax		
			
		ymax = 1d0
		ymin = 0d0

		xmax = a
		xmin = 0d0
		
		! the world coords of array point A(I,J) are given by
		! X = TR(1) + TR(2)*I + TR(3)*J
        ! Y = TR(4) + TR(5)*I + TR(6)*J		
		transform(1) = -1d0 * dx
		transform(2) = dx
		transform(3) = 0d0
		transform(4) = -1d0 * dz
		transform(5) = 0d0
		transform(6) = dz
			
		call pgbbuf ! begin buffer of pgplot commands
		call pgsci(1) ! white for frame, labels, etc
		CALL pgenv(xmin,xmax,ymin,ymax,0,1)
		
		call pgask(.false.) ! whether or not pgplot pauses and prompts

		! contours for drawing magnetic field lines

		call pgscir(16,16+num_levels-1)
		if (plot_type.eq.1) then
			allocate(z_contours(num_mag_contours))
			do i=1,num_mag_contours
				alfa = 1d0 * (i-1) / num_mag_contours
				z_contours(i) = (1-alfa) * zmin + alfa * zmax
			end do
			call pgcont(zvar,nx,nz,1,nx,2,nz,z_contours,num_mag_contours,transform)
		else
			call pgimag(zvar,nx,nz,1,nx,1,nz,zmin,zmax,transform)
		end if
		call pgsch(1.0)
        CALL pglab('x', 'z','')
		call pgmtext('t',0.5,0.5,0.5,annot)
		call pgsch(1.5)
		call pgmtext('t',1.5,0.5,0.5,trim(plot_name))			
		call pgsch(1.0)			

		
!		call pgslw(4) ! ~default thickness
		call pgebuf ! dump buffer

		
	end subroutine plot

end module fluid

program main

	use fluid
	
 	use omp_lib

	implicit none	

	character(100) :: action
	character(100) :: param
	character(100) :: png_count_c
	integer :: photo_num

	call getarg(1,action)
	call getarg(2,param)
	call getarg(3,png_count_c)
		
	if (action .eq. 're') then
		read(param,*) photo_num
		read(png_count_c,*) png_count
		do_load_state = .true.
		if (photo_num.lt.10) write(infile,'(a10,i1)') 'photos/000', photo_num
		if ((photo_num.ge.10).and.(photo_num.lt.100)) write(infile,'(a9,i2)') 'photos/00', photo_num
		if ((photo_num.ge.100).and.(photo_num.lt.1000)) write(infile,'(a8,i3)') 'photos/0', photo_num
		if (photo_num.ge.1000) write(infile,'(a7,i4)') 'photos/', photo_num
	else if (action .eq. 'tau') then
		read(param,*) tau
		write(*,*) 'using tau = ', tau
! 		do_pgplot_xwin = .false.
	end if
	if (action.eq.'quiet') do_pgplot_xwin = .false.
	
	call read_inlist
	call initialize
	call save_state ! always save initial state in photos/0000
	call evolve
	call cleanup

end program main
