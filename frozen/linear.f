module fluid

	integer, parameter :: nz = 100, nm = 5 ! number of zones, number of fourier modes
	integer :: ktrack = int(nz/2.) ! track the midpoint in z
	integer :: step=0, steps_to_take = 20000
	double precision, dimension(:), pointer :: z, sub, dia, sup, work1, work2 ! nz
	double precision, dimension(:,:), pointer :: stream, omega, temp, &
		temp_old, omega_old, stream_old, temp_zz, omega_zz ! nz, nm
	double precision, dimension(:,:,:), pointer :: tdot, omegadot ! nz, nm, 2 (t_old and t)
	double precision :: pi = 4.*atan(1.), dz, dt, t=0d0, ra, pr, a, c
	character(len=100) :: j
	logical :: do_pgplot = .true., do_linear_only = .false.
		
	contains
	
	subroutine initialize
		! does what it sounds like
		integer :: k, ierr, pgopen
		allocate(z(nz), sub(nz), dia(nz), sup(nz), work1(nz), work2(nz))
		allocate(stream(nz,nm),omega(nz,nm),temp(nz,nm),temp_old(nz,nm), &
			temp_zz(nz,nm), omega_zz(nz,nm), omega_old(nz,nm), stream_old(nz,nm))
		allocate(omegadot(nz,nm,2),tdot(nz,nm,2))
		
		ra = 100 ! g * alpha * DeltaT * D**3 / (nu * kappa)
		pr = 1d0 ! nu/kappa
		a = 1.618 ! aspect ratio (width/height)
		c = pi / a
		
		dz = 1. / (nz-1)
		
		! timestep must resolve a thermal time and a viscous diffusion time
		dt = 0.9 * dz**2 / 4 / pr

		! zero fluid velocity at t=0
		stream(:,:) = 0d0
		omega(:,:) = 0d0
		
		do k=1,nz
			z(k) = (k-1) * dz
			temp(k,:) = sin(pi * z(k))
		end do
		sub(:) = -1./dz**2
		sup(:) = -1./dz**2
		! streamfunction = 0 at k=1 and k=nz
		sub(nz) = 0d0
		sup(1) = 0d0
		
		temp_zz(:,:) = 0d0
		omega_zz(:,:) = 0d0
		
		tdot(:,:,1) = 0d0 ! initial "old" tdot
		omegadot(:,:,1) = 0d0 ! initial "old" omegadot

		if (do_pgplot) then
			ierr = pgopen('/xwin') !pgbeg(0,'?',1,1)
	        if (ierr .ne. 1) then
				write(*,*) 'pgopen failed in initialize'
				stop
			end if
		end if
		
	end subroutine initialize
	
	subroutine cleanup
		deallocate(z, sub, dia, sup, work1, work2, stream, omega, temp, &
		temp_old, omega_old, stream_old, omegadot, tdot, temp_zz, omega_zz)
		if (do_pgplot) call pgend()
	end subroutine cleanup
	
	subroutine tridiag(nz, rhs, sol, sub, dia, sup, work1, work2)
		! tridiagonal matrix solver
		double precision, dimension(:) :: rhs, sol, sub, dia, sup, work1, work2
		integer :: i
				
		work1(1) = 1. / dia(1)
		work2(1) = sup(1) * work1(1)
		do i=1, nz-1
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
	
    subroutine init_random_seed()
      use iso_fortran_env, only: int64
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid
      integer(int64) :: t
    
      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24_int64 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
         end if
         pid = getpid()
         t = ieor(t, int(pid, kind(t)))
         do i = 1, n
            seed(i) = lcg(t)
         end do
      end if
      call random_seed(put=seed)
    contains
      ! This simple PRNG might not be good enough for real work, but is
      ! sufficient for seeding a better PRNG.
      function lcg(s)
        integer :: lcg
        integer(int64) :: s
        if (s == 0) then
           s = 104729
        else
           s = mod(s, 4294967296_int64)
        end if
        s = mod(s * 279470273_int64, 4294967291_int64)
        lcg = int(mod(s, int(huge(0), int64)), kind(0))
      end function lcg
    end subroutine init_random_seed	
	subroutine test_tridiag
		integer :: j, test_nz = 10000
		double precision, dimension(:), pointer :: rhs, sol, expected_sol, wk1, wk2, &
			dia, sub, sup
		
		allocate(rhs(test_nz),sol(test_nz),expected_sol(test_nz), &
			wk1(test_nz), wk2(test_nz), dia(test_nz), sub(test_nz), sup(test_nz))
		
! 			dia = (/ 0.57141397, 0.10439855, 0.37844050, 0.30001649, 0.48309448 /)
! 			sub = (/ 0.00000000, 0.04690249, 0.08006328, 0.09211752, 0.01855902 /)
! 			sup = (/ 0.09923425, 0.05383660, 0.01701892, 0.02812829, 0.00000000 /)
! 
! 			expected_sol = (/ 0.33035731, 0.62151472, 0.33130685, 0.73153078, 0.14220734 /)
! 
! 			rhs = (/ 0.25044633, 0.09821625, 0.18759030, 0.25399052, 0.08227607 /)
			
! 			call init_random_seed()

			call init_random_seed
			call random_number(dia)
			call random_number(sub)
			call random_number(sup)
			sub = sub * 0.1
			sup = sup * 0.1
			call random_number(expected_sol)
			
			rhs(:) = 0d0
			
			rhs(1) = dia(1) * expected_sol(1) + sup(1) * expected_sol(2)
			do j=2,test_nz-1
				rhs(j) = rhs(j) + dia(j) * expected_sol(j) &
					+ sub(j) * expected_sol(j-1) &
					+ sup(j) * expected_sol(j+1)
			end do
			rhs(test_nz) = dia(test_nz) * expected_sol(test_nz) &
				+ sub(test_nz) * expected_sol(test_nz-1)	
				
				sol(:) = 0d0

! 			write(*,'(a,5es16.8)') 'dia:', dia
! 			write(*,'(a,5es16.8)') 'sub:', sub
! 			write(*,'(a,5es16.8)') 'sup:', sup
! 			write(*,*)
! 			write(*,'(a,5es16.8)') 'rhs:', rhs
		
		
			call tridiag(test_nz,rhs,sol,sub,dia,sup,wk1,wk2)

			write(*,*) 
! 			write(*,'(a,5es16.8,a)') 'sol:', sol, ' (routine)'
! 			write(*,'(a,5es16.8,a)') 'sol:', expected_sol, ' (expected)'
			write(*,'(a,es16.8)') 'maximum absolute difference: ', maxval(expected_sol - sol)
		
		deallocate(dia,sub,sup,rhs,sol,expected_sol)
		
		
	end subroutine test_tridiag
		
	subroutine report
		write(*,*)
		write(*,'(2a16,a16,a16)') 'step', 'time', 'ktrack', 'nz'
		write(*,'(i16,f16.8,2i16)') step, t, ktrack, nz
		
		write(*,'(a3,6a16)') 'n', 'temperature', 'omega', 'stream', &
			'dln_temp', 'dln_omega', 'dln_stream'
		do n=1,nm
			write(*,'(i3,7es16.8)') n, &
				temp(ktrack,n), omega(ktrack,n), stream(ktrack,n), &
				log(abs(temp(ktrack,n))) - log(abs(temp_old(ktrack,n))), &
				log(abs(omega(ktrack,n))) - log(abs(omega_old(ktrack,n))), &
				log(abs(stream(ktrack,n))) - log(abs(stream_old(ktrack,n)))
		end do
		write(*,*)
	end subroutine report
	
	
	subroutine report_final
		double precision :: dlnT(nm), dlnO(nm), dlnP(nm)
		dlnT = log(abs(temp(ktrack,:))) - log(abs(temp_old(ktrack,:)))
		dlnO = log(abs(omega(ktrack,:))) - log(abs(omega_old(ktrack,:)))
		dlnP = log(abs(stream(ktrack,:))) - log(abs(stream_old(ktrack,:)))
		open(unit=17,file='linear_outcomes.out',access='append')
! 		write(17,'(11a16)') 'Ra', 'Pr', 'dln_temp_n1', 'dln_omega_n1', 'dln_stream_n1', &
! 			'dln_temp_n2', 'dln_omega_n2', 'dln_stream_n2', &
! 			'dln_temp_n3', 'dln_omega_n3', 'dln_stream_n3'
		write(17,'(11es16.8)') ra, pr, dlnT(1), dlnO(1), dlnP(1), &
			dlnT(2), dlnO(2), dlnP(2), dlnT(3), dlnO(3), dlnP(3)
		close(17)
	end subroutine report_final
	
	subroutine compute_temp_advection
		return
	end subroutine compute_temp_advection
	
	subroutine compute_omega_advection
		return
	end subroutine compute_omega_advection
		
	subroutine evolve
		integer :: k, n
		
		do while (step .le. steps_to_take) ! step through time
		
			! update tdot and omegadot (Eq. 3.3, 3.4, and 2.16)
			do k=2,nz-1
				do n=1,nm
					if (isnan(temp(k,n))) stop 'got NaN in temperature'
					
					! central finite differences to get second derivatives wrt z
					temp_zz(k,n) = (temp(k+1,n) - 2 * temp(k,n) + temp(k-1,n)) / dz**2
					omega_zz(k,n) = (omega(k+1,n) - 2 * omega(k,n) + omega(k-1,n)) / dz**2
															
					! eqs. 3.3 and 3.4
					tdot(k,n,2) = temp_zz(k,n) - (n * c)**2*temp(k,n)
					! add linearized advection term if just doing linear case
					if (do_linear_only) tdot(k,n,2) = tdot(k,n,2) + n * c * stream(k,n)

					omegadot(k,n,2) = ra * pr * n * c * temp(k,n) &
						+ pr * (omega_zz(k,n) - (n * c)**2 * omega(k,n))

				end do
			end do
! 			write(*,'(3es16.8)') ra, pr, c
! 			write(*,'(8es16.8)') omegadot(:,1,2)
! 			write(*,'(8es16.8)') temp(:,1)
! 			write(*,*)

			if (mod(step,500).eq.0) then
				call report
				if(do_pgplot) call plot
			end if			
					
			temp_old = temp
			omega_old = omega
			stream_old = stream				
			! update temp and omega (Eq. 2.18)
			do k=2,nz-1
				do n=1,nm
					temp(k,n) = temp(k,n) + dt / 2 * (3 * tdot(k,n,2) - tdot(k,n,1))
					omega(k,n) = omega(k,n) + dt  / 2 * (3 * omegadot(k,n,2) - omegadot(k,n,1))
				end do
			end do
		
			! contruct dia
			! call tridiag with omega as rhs to get stream=sol
			! laplacian(psi_n) + omega_n = 0
			do n=1,nm
				dia(1) = 1
				dia(nz) = 1
				dia(2:nz-1) = (c * n)**2 + 2. / dz**2
								
! 				if (n.eq.1) then
! 					write(*,'(5es16.8)') dia(1), dia(ktrack), dia(nz), sub(1), sup(nz)
! 				end if
				
				! tridiag takes (nz, rhs, sol, sub, dia, sup, work1, work2)
				call tridiag(nz,omega(:,n),stream(:,n),sub,dia,sup,work1,work2)
			end do
					
			! advance time
			tdot(:,:,1) = tdot(:,:,2)
			t = t + dt
			step = step + 1
			
			! check that all fourier cofficients vanish on boundaries
! 			if (mod(step,100).eq.0) then
! 				write(*,*) 'am I zero?'
! 				do n=1,nm
! 					write(*,'(6es16.8)') temp(1,n), temp(nz,n), omega(1,n), omega(nz,n), stream(1,n), stream(nz,n)
! 				end do
! 			end if
						
		end do
		
! 		call report_final
	end subroutine evolve
		
	subroutine plot
        integer :: i
		character(100) :: title
		real :: ymin, ymax, xvar(nz), y1(nz), y2(nz)
		
		write(title,('(a5,i8,a10,es16.8)')) 'step ', step, 'time ', t
		
		xvar = real(z)
 		y1 = real(temp(:,1)) ! temp(:,1)
		y2 = real(temp_zz(:,1))
		
		ymax = 1.05 * maxval(y1)
		ymin = 1.05 * minval(y1)

! 		call pgsvp(0.0,1.0,0.5,1.0) ! set viewport

		call pgbbuf
		call pgsci(1)
		CALL PGENV(minval(xvar),maxval(xvar),ymin,ymax,0,1)
        CALL pglab('z', 'temperature',title)
	
		call pgask(.false.) ! whether or not pgplot pauses and prompts

		call pgsci(11) ! blue
		call pgpt(nz,xvar,y1,9) ! number of points, x points, y points, marker style
		call pgsci(13) ! pinkish
		call pgslw(6) ! thicker line
		call pgline(nz,xvar,y1)
		
! 		call pgslw(4) ! ~default thickness
! 		call pgsci(4)
! 		call pgpt(nz,xvar,y2,8) ! number of points, x points, y points, marker style
! 		call pgsci(2)
! 		call pgslw(6) ! thicker line
! 		call pgline(nz,xvar,y2)
		
		call pgslw(4) ! ~default thickness
		call pgebuf
	end subroutine plot

end module fluid

program main

	use fluid
	
	character(100) :: action
	
	call getarg(1,action)	
		
	if (action.eq.'test') then
		call test_tridiag
	else
		if (action.eq.'quiet') do_pgplot = .false.
		call initialize
		call evolve
		call cleanup
	end if

end program main
