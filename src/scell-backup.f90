! To simulate off-lattice cells as spheres, following Drasdo & Hohme
!
! Biomechanics, includung cell motion, is treated by the methods described in:
!
! Metropolis N, Rosenbluth A W, Rosenbluth M N, Teller A H and Teller E 1953 Equation of state calculations by fast computing
! machines J. Chem. Phys. 21 1087–92
!
! Drasdo D, Kree R and McCaskill J S 1995 Monte-carlo approach to tissue-cell populations Phys. Rev. E 52 6635–57
!
! D. Drasdo and S. Hoehme, A single-cell based model to tumor growth in-vitro: Monolayers and spheroids. Phys. Biol. 2:133 (2005).
! 
!	The state S; of a cell i is described by its position, its age class, and its shape. The cell position is given
!	by the vector r, of its center of mass. For an I cell, the shape is characterized by its radius R, whereas for
!	an M cell it is the vector d(i), joining the centers of the corresponding dumbbell which fixes its shape. As the
!	radii R(1,i) = R(2,i) = R(i), of the spheres making up the dumbbell are equal and furthermore determined by the
!	length of d(i), we may also describe the shape by R(i) and a unit vector d(i), which fixes the orientation of an M cell in
!	space.
!	During a single updating step, only one randomly chosen cell is considered as active and tries to perform one
!	of the following actions: (i) action_I (migration, growth) for I cells, (ii) action_M (migration, rotation, deformation)
!	for M cells.
!	A single migration trial consists of a shift r(i) -> r(i) + dr(i), of cell i in a random direction with a step length
!	which is chosen at random out of an interval [O, dr_max].
!	A rotation trial is defined analogously. During a growth step, the radius of an I cell is increased by a random
!	amount dR in [O, dR_max] whereas during a deformation step the length of d(i) is increased and at the same time
!	the dumbbell radius is decreased such that the volume of the cell is kept constant. All these trials are described in
!	detail in the Appendix.
!	Whether cell i actually performs its chosen trial or not depends on its interaction with other cells. Following
!	the discussion of the preceding subsection, we model the interaction between two cells i and j by potentials
!	V(S(i), S(j)). Note that these interaction potentials may be different for different types of actions. We will discuss
!	the biological basis of this complication below.
!	The cell dynamics is based on the assumption that a trial is accepted according to the standard prescription
!	of the Metropolis algorithm [49], i.e., an action a is always accepted if the change of the state S(i) -> S^(i) of cell i
!	induced by the action, leads to a decrease of its interaction potential with other cells:
!       E(S(i)) = Sum {V(S(i),S(j))}, 
!   whereas it is only accepted with probability 
!       exp(—[E(S^(i))—E(S(i))]) if E(S(i)) is increased by the action.
!
! Plan:
! First test motion of a collection of cells at various states (fixed), i.e. a mix of spheres and dumbells.

module scell

use global
implicit none

contains

!-----------------------------------------------------------------------------------------
! d~ = d - (R1+R2)
! d_detach is the value of d~ at which V -> 0, i.e. it depends on R1+R2
!-----------------------------------------------------------------------------------------
subroutine setup(ncpu, infile, outfile, ok)
integer :: ncpu
character*(128) :: infile, outfile
logical :: ok
type(cell_type), pointer :: cp
real(REAL_KIND) :: k_v, v, R0

ok = .true.
#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif
Mnodes = ncpu
seed = [1234,5678]

call rng_initialisation

nsteps = 0
R0 = 5
alpha_v = 0.25
epsilon = 7.5
es_e = 1
shift = -6
Dfactor = 0.01
sqr_es_e = sqrt(es_e)
dr_grow = 0.002
dr_max = 100*dr_grow	! 0.2
ntries = 10
k_v = 2/alpha_v - sqr_es_e + sqrt(es_e - shift/epsilon)
k_detach = k_v*alpha_v/2
!write(*,'(a,f8.4)') 'k_detach: ',k_detach

allocate(cell_list(1000))
ncells = 8
cp=>cell_list(1)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [0,0,0]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0
cp=>cell_list(2)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [8,0,0]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0
cp=>cell_list(3)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [0,8,0]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0
cp=>cell_list(4)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [8,8,0]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0
cp=>cell_list(5)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [0,0,9]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0
cp=>cell_list(6)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [8,0,9]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0
cp=>cell_list(7)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [0,8,9]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0
cp=>cell_list(8)
cp%Iphase = .true.
cp%radius(1) = R0
cp%centre(1,:) = [8,8,9]
allocate(cp%nbrlist(100,2))
cp%nbrs = 0

call setup_nbrlists
call testrun
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testrun
integer :: it, kcell, nits, k, i
integer :: kpar=0
real(REAL_KIND) :: R1, c1(3), E0, E1, dE, dr(3), p, Emin, drmin(3)
type(cell_type), pointer :: cp1

call testresults
nits = 1000
do it = 1,nits
	if (mod(it,100) == 0) call setup_nbrlists
	do k = 1,ncells
		kcell = mod(it+k,ncells) + 1
		cp1 => cell_list(kcell)
		if (cp1%Iphase) then
			R1 = cp1%radius(1)
			c1 = cp1%centre(1,:)
		endif
		write(*,'(a,i4,a,3f6.2,a,f6.2)') 'kcell: ',kcell,' centre: ',c1,' R: ',R1
		E0 = Efun(cp1,R1,c1)
		if (E0 >= BIG/2) then
			write(*,*) 'Bad E0: ',it,E0
			stop
		endif
		Emin = 100*BIG
		do i = 1,ntries
			call get_random_dr(dr)
			c1 = cp1%centre(1,:) + dr
			E1 = Efun(cp1,R1,c1)
			if (E1 < Emin) then
				Emin = E1
				drmin = dr
			endif
!			write(*,'(a,3f6.2,2e12.3)') 'trial: ',c1,E1
!			if (E1 < BIG/2) exit
		enddo
		dE = Emin - E0
!		write(*,'(a,f8.3)') 'dE: ',dE
		if (dE >= 0) then
			p = Dfactor*exp(-dE)
			if (p > par_uni(kpar)) cycle
		endif
		cp1%centre(1,:) = cp1%centre(1,:) + drmin
!		write(*,'(a,i4,3f6.2)') 'kcell, centre: ',kcell,c1
	enddo
	if (it <= 1000) call grower
enddo
call testresults
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testresults
integer :: k1, k2
real(REAL_KIND) :: R1, R2, c1(3), c2(3), v(3), d2, d

write(*,*) 'testresults'
do k1 = 1,ncells
	R1 = cell_list(k1)%radius(1)
	c1 = cell_list(k1)%centre(1,:)
	write(*,'(a,i4,3f6.2)') 'k1, c: ',k1,c1
	do k2 = k1+1, ncells
		R2 = cell_list(k2)%radius(1)
		c2 = cell_list(k2)%centre(1,:)
		v = c2 - c1
		d2 = dot_product(v,v)
		d = sqrt(d2)
		write(*,'(2i4,f8.2)') k1,k2,d
	enddo
enddo

end subroutine
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine grower
integer :: kcell
real(REAL_KIND) :: dr

dr = dr_grow
do kcell = 1,ncells
	cell_list(kcell)%radius(1) = cell_list(kcell)%radius(1) + dr
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function Efun(cp1,R1,c1) result(E)
type(cell_type), pointer :: cp1
real(REAL_KIND) :: R1, c1(3)
type(cell_type), pointer :: cp2
real(REAL_KIND) :: R2, c2(3)
integer :: k, knbr
logical :: incontact

E = 0
do k = 1,cp1%nbrs
	knbr = cp1%nbrlist(k,1)
	incontact = (cp1%nbrlist(k,2) == 1)
	cp2 => cell_list(knbr)
	if (cp2%Iphase) then
		R2 = cp2%radius(1)
		c2 = cp2%centre(1,:)
	endif
!	write(*,'(a,i4,3f6.2)') 'knbr: ',knbr,c2
	E = E + Vfun(R1,c1,R2,c2,incontact)
	if (incontact) then
		cp1%nbrlist(k,2) = 1
	else
		cp1%nbrlist(k,2) = 0
	endif
enddo
end function

!-----------------------------------------------------------------------------------------
! dd is d~
! Note that incontact may be changed here, both ways, implementing the hysteresis loop
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function Vfun(R1,c1,R2,c2,incontact)
real(REAL_KIND) :: R1, R2, c1(3), c2(3)
logical :: incontact
real(REAL_KIND) :: v(3), d, dd, Rsum, delta	!,Rc, d1, Ac

Vfun = 0
v = c2 - c1
d = sqrt(dot_product(v,v))
Rsum = R1+R2
dd = d - Rsum
delta = alpha_v*Rsum
if (d > k_detach*Rsum) then
	incontact = .false.
elseif (.not. incontact .and. dd < 0) then
	incontact = .true.
endif
! find Rc = radius of contact circle, Ac = area of contact circle
!if (d >= R1 + R2) then
!	Rc = 0
!else
!	d1 = 0.5*(d + (R1**2 - R2**2)/d)
!	write(*,*) 'd1: ',d1
!	Rc = sqrt(R1**2 - d1**2)
!endif
!write(*,*) 'Rc: ',Rc
!Ac = PI*Rc**2
!write(*,*) 'Ac: ',Ac

if (dd < -delta) then
	Vfun = BIG
elseif (.not.incontact) then
	Vfun = 0
else
	Vfun = epsilon*((2*dd/delta + sqr_es_e)**2 - es_e) + shift
endif
!write(*,'(a,2f6.2)') 'd,dd: ',d,dd

end function

!-----------------------------------------------------------------------------------------
! Dumb version to start
!-----------------------------------------------------------------------------------------
subroutine setup_nbrlists
type(cell_type), pointer :: cp1, cp2
real(REAL_KIND) :: R1, c1(3), R2, c2(3), v(3), d2, d, Rsum, dlim, dfactor
integer :: kcell, k2, k, nbrs, nbrlist(100,2)

do kcell = 1,ncells
	cp1 => cell_list(kcell)
	R1 = cp1%radius(1)
	c1 = cp1%centre(1,:)
	nbrs = 0
	do k2 = 1,ncells
		if (k2 == kcell) cycle
		cp2 => cell_list(k2)
		R2 = cp2%radius(1)
		c2 = cp2%centre(1,:)
		Rsum = R1 + R2
		dlim = 1.5*Rsum
		v = c2 - c1
		d2 = dot_product(v,v)
		d = sqrt(d2)
		if (d < dlim) then
			dfactor = 1
			do k = 1,cp1%nbrs
				if (cp1%nbrlist(k,1) == k2) then
					if (cp1%nbrlist(k,2) == 1) dfactor = k_detach
					exit
				endif
			enddo
			nbrs = nbrs + 1
			nbrlist(nbrs,1) = k2
			if (d < dfactor*Rsum) then
				nbrlist(nbrs,2) = 1
			else
				nbrlist(nbrs,2) = 0
			endif
		endif
	enddo
	cp1%nbrs = nbrs
	cp1%nbrlist(1:nbrs,:) = nbrlist(1:nbrs,:)
	write(*,*) 'kcell: ',kcell,nbrs
	write(*,*) nbrlist(1:nbrs,:)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! This needs to be given radial symmetry
!-----------------------------------------------------------------------------------------
subroutine get_random_dr(dr)
real(REAL_KIND) :: dr(3)
integer :: kpar=0

dr(1) = dr_max*2*(par_uni(kpar) - 0.5)
dr(2) = dr_max*2*(par_uni(kpar) - 0.5)
dr(3) = dr_max*2*(par_uni(kpar) - 0.5)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation
integer, allocatable :: zig_seed(:)
integer :: i, n, R
integer :: kpar = 0
integer :: npar, grainsize = 32

npar = Mnodes
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.

!n = 0
!do i = 1,1000000000
!	R = par_shr3(kpar)
!	if (R == -2147483648) n = n+1
!enddo
!write(*,*) 'n = ',n
!stop
end subroutine



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res

res=0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!subroutine Execute() BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
integer :: i, res
logical :: ok, isopen

infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

awp_0%is_open = .false.
awp_1%is_open = .false.
par_zig_init = .false.
logfile = 'scell.log'
inquire(unit=nflog,OPENED=isopen)
if (.not.isopen) then
    open(nflog,file=logfile,status='replace')
endif
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
if (use_TCP) then
	write(nflog,*) 'call connecter'
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call Setup(ncpu,infile,outfile,ok)
if (ok) then
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
	res = 0
else
	call logger('=== Setup failed ===')
	res = 1
	stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from spheroid_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
else
	call logger('  === Execution failed ===')
endif
!write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

!close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
!call logger('deallocated all arrays')

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine


end module