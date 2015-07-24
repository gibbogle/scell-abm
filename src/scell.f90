!---------------------------------------------------------------------------------------------------------------
! Force-based simulation
!
! UNITS: (OLD)
!    distance um
!    time     seconds
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!---------------------------------------------------------------------------------------------------------------
module scell

use global
use fmotion
use nbr
use chemokine
use react_diff

#include "../src/version.h"

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok

call RngInitialisation
! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends 
! it will still be possible to view the cell distributions and chemokine concentration fields.
if (allocated(cell_list)) deallocate(cell_list)
if (allocated(grid)) deallocate(grid)
if (allocated(perm_index)) deallocate(perm_index)
if (allocated(Cextra)) deallocate(Cextra)
if (allocated(Caverage)) deallocate(Caverage)
if (allocated(Cflux)) deallocate(Cflux)
if (allocated(Cflux_prev)) deallocate(Cflux_prev)
if (allocated(stencil)) deallocate(stencil)
call logger('did deallocation')
allocate(cell_list(MAX_NLIST))
allocate(grid(NX,NY,NZ))
allocate(perm_index(MAX_NLIST))
allocate(Cextra(NX,NY,NZ,NCONST))
allocate(Caverage(NX,NY,NZ,NCONST))
allocate(Cflux(NX,NY,NZ,NCONST))
allocate(Cflux_prev(NX,NY,NZ,NCONST))
allocate(stencil(NX,NY,NZ))

ncells = 0
nlist = 0
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
integer :: iV_depend, iV_random, Nmm3, itestcase, ictype, ishow_progeny
integer :: iuse_oxygen, iuse_glucose, iuse_tracer, iconstant
real(REAL_KIND) :: days, percent, fluid_fraction, medium_volume0, d_layer, sigma, Vsite_cm3

ok = .true.
!chemo(:)%used = .false.

open(nfcell,file=inputfile,status='old')

read(nfcell,*) gui_run_version				! program run version number
read(nfcell,*) dll_run_version				! DLL run version number
read(nfcell,*) NX							! size of fine grid
read(nfcell,*) NXB							! size of coarse grid
read(nfcell,*) DELTA_X						! grid size (um)
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median			! hours
read(nfcell,*) divide_time_shape
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) a_separation
read(nfcell,*) a_force
read(nfcell,*) c_force
read(nfcell,*) x0_force
read(nfcell,*) x1_force
read(nfcell,*) kdrag
read(nfcell,*) frandom
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
read(nfcell,*) Nmm3							! number of cells/mm^3
read(nfcell,*) fluid_fraction				! fraction of the (non-necrotic) tumour that is fluid
read(nfcell,*) medium_volume0				! initial total volume (medium + spheroid) (cm^3)
read(nfcell,*) d_layer						! thickness of the unstirred layer around the spheroid (cm)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) ANOXIA_THRESHOLD			    ! O2 threshold for anoxia (uM)
read(nfcell,*) anoxia_tag_hours				! hypoxic time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_input					! for GUI just a placeholder for ncpu, used only when execute parameter ncpu = 0
read(nfcell,*) Ncelltypes					! maximum number of cell types in the spheroid
do ictype = 1,Ncelltypes
	read(nfcell,*) percent
	celltype_fraction(ictype) = percent/100
enddo
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) show_progeny                 ! if != 0, the number of the cell to show descendents of

read(nfcell,*) iuse_oxygen		! chemo(OXYGEN)%used
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%medium_diff_coef
read(nfcell,*) chemo(OXYGEN)%membrane_diff_in
Vsite_cm3 = 2.0e-9		! this is from spheroid-abm, and the scaling here is for consistency with spheroid-abm
chemo(OXYGEN)%membrane_diff_in = chemo(OXYGEN)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(OXYGEN)%membrane_diff_out = chemo(OXYGEN)%membrane_diff_in
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) iconstant
chemo(OXYGEN)%constant = (iconstant == 1)
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
chemo(OXYGEN)%max_cell_rate = chemo(OXYGEN)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
read(nfcell,*) chemo(OXYGEN)%MM_C0
read(nfcell,*) chemo(OXYGEN)%Hill_N
read(nfcell,*) iuse_glucose		!chemo(GLUCOSE)%used
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%medium_diff_coef
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_in
chemo(GLUCOSE)%membrane_diff_in = chemo(GLUCOSE)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(GLUCOSE)%membrane_diff_out = chemo(GLUCOSE)%membrane_diff_in
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
read(nfcell,*) iconstant
chemo(GLUCOSE)%constant = (iconstant == 1)
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate
chemo(GLUCOSE)%max_cell_rate = chemo(GLUCOSE)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
read(nfcell,*) chemo(GLUCOSE)%MM_C0
read(nfcell,*) chemo(GLUCOSE)%Hill_N
read(nfcell,*) iuse_tracer		!chemo(TRACER)%used
read(nfcell,*) chemo(TRACER)%diff_coef
read(nfcell,*) chemo(TRACER)%medium_diff_coef
read(nfcell,*) chemo(TRACER)%membrane_diff_in
chemo(TRACER)%membrane_diff_in = chemo(TRACER)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(TRACER)%membrane_diff_out = chemo(TRACER)%membrane_diff_in
read(nfcell,*) chemo(TRACER)%bdry_conc
read(nfcell,*) iconstant
chemo(TRACER)%constant = (iconstant == 1)
read(nfcell,*) chemo(TRACER)%max_cell_rate
read(nfcell,*) chemo(TRACER)%MM_C0
read(nfcell,*) chemo(TRACER)%Hill_N

close(nfcell)

if (chemo(OXYGEN)%Hill_N /= 1 .and. chemo(OXYGEN)%Hill_N /= 2) then
	call logger('Error: OXYGEN_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
if (chemo(GLUCOSE)%Hill_N /= 1 .and. chemo(GLUCOSE)%Hill_N /= 2) then
	call logger('Error: GLUCOSE_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
DELTA_X = 1.0e-4*DELTA_X	! um -> cm
dxf = DELTA_X
dxb = NRF*dxf
NY = NX
NZ = NX
NYB = NXB
NZB = NXB
Kdrag = 1.0e5*Kdrag
MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
ANOXIA_THRESHOLD = ANOXIA_THRESHOLD/1000			! uM -> mM
!O2cutoff = O2cutoff/1000							! uM -> mM
chemo(OXYGEN)%used = (iuse_oxygen == 1)
chemo(GLUCOSE)%used = (iuse_glucose == 1)
chemo(TRACER)%used = (iuse_tracer == 1)
chemo(OXYGEN)%MM_C0 = chemo(OXYGEN)%MM_C0/1000		! uM -> mM
chemo(GLUCOSE)%MM_C0 = chemo(GLUCOSE)%MM_C0/1000	! uM -> mM
t_anoxic_limit = 60*60*anoxia_tag_hours				! hours -> seconds
anoxia_death_delay = 60*60*anoxia_death_hours		! hours -> seconds
nsteps = days*24*3600./DELTA_T
write(logmsg,*) 'nsteps: ',nsteps
call logger(logmsg)

divide_dist%class = LOGNORMAL_DIST
divide_time_median = 60*60*divide_time_median		! hours -> seconds
sigma = log(divide_time_shape)
!divide_dist%p1 = log(divide_time_mean/exp(sigma*sigma/2))	
divide_dist%p1 = log(divide_time_median)	
divide_dist%p2 = sigma
divide_time_mean = exp(divide_dist%p1 + 0.5*divide_dist%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,2e12.4)') 'shape, sigma: ',divide_time_shape,sigma
call logger(logmsg)
write(logmsg,'(a,2e12.4)') 'Median, mean divide time: ',divide_time_median/3600,divide_time_mean/3600
call logger(logmsg)
use_V_dependence = (iV_depend == 1)
!randomise_initial_volume = (iV_random == 1)
!use_extracellular_O2 = (iuse_extra == 1)
!randomise_initial_volume = (iV_random == 1)
!relax = (iuse_relax == 1)
!use_parallel = (iuse_par_relax == 1)


end subroutine

!-----------------------------------------------------------------------------------------
! d~ = d - (R1+R2)
! d_detach is the value of d~ at which V -> 0, i.e. it depends on R1+R2
!-----------------------------------------------------------------------------------------
subroutine setup(ncpu, infile, outfile, ok)
integer :: ncpu
character*(128) :: infile, outfile
logical :: ok
integer :: kcell
real(REAL_KIND) :: k_v

ok = .true.
initialized = .false.
par_zig_init = .false.

inputfile = infile
outputfile = outfile
call logger("ReadCellParams")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

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

NY = NX
NZ = NX
Ndim(1) = NX
Ndim(2) = NY
Ndim(3) = NZ
call ArrayInitialisation(ok)
if (.not.ok) return
call logger('did ArrayInitialisation')

!call create_stencil

call SetupChemo
call logger('did SetupChemo')
!call SetInitialConcs

! Assume that Raverage is for cells midway in growth, i.e. between 0.8 and 1.6, at 1.2
! as a fraction of the average cell volume, (or, equivalently, between 1.0 and 2.0, at 1.5)
! and that actual divide volume Vdivide is: Vdivide0-dVdivide < Vdivide < Vdivide0+dVdivide
! Note that 2.0/1.5 = 1.6/1.2
! Growth occurs at a constant rate for (divide_time - mitosis_duration)
! This is option A.
! Try a different option: B
Vdivide0 = (2.0/1.5)*(4.*PI/3.)*Raverage**3
dVdivide = 0.3*Vdivide0	
Rdivide0 = Raverage*(2.0/1.5)**(1./3.)
d_nbr_limit = 1.5*2*Rdivide0	! 1.5 is an arbitrary choice - was 1.2

if (MITOSIS_MODE == TERMINAL_MITOSIS) then
	mitosis_duration = 2*3600
	test_growthrate = Vdivide0/(2*(divide_time_median - mitosis_duration))	! um3/sec
else
	mitosis_duration = 0
	test_growthrate = Vdivide0/(2*divide_time_median)	! um3/sec
endif
write(logmsg,'(a,2f8.2,e12.3)') 'Raverage,Vdivide0,test_growthrate: ',Raverage,Vdivide0,test_growthrate
call logger(logmsg)

call setup_force_parameters

call make_jumpvec
call PlaceCells(ok)
call logger('did PlaceCells')
!do kcell = 1,nlist
!	write(nflog,'(i6,3f8.1)') kcell,1.0e4*cell_list(kcell)%centre(:,1)
!enddo
alpha_v = 0.25
epsilon = 7.5
es_e = 1
shift = -6
Dfactor = 1
sqr_es_e = sqrt(es_e)

ntries = 1
t_fmover = 0
ndtotal = 0
ndtotal_last = 0

k_v = 2/alpha_v - sqr_es_e + sqrt(es_e - shift/epsilon)
k_detach = k_v*alpha_v/2
!write(*,'(a,f8.4)') 'k_detach: ',k_detach

call setup_nbrlists
call logger('did setup_nbrlists')
call setup_react_diff
! For simple testing...
!cell_list(:)%Cin(OXYGEN) = chemo(OXYGEN)%bdry_conc
!cell_list(:)%Cin(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
call logger('did setup_react_diff')
call logger('completed Setup')
end subroutine



!-----------------------------------------------------------------------------------------
! Cells are initially placed in a regular rectangular grid pattern, spacing between cell
! centres equal to 2*Raverage
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: ix, iy, iz, kcell, site(3), irad, lastid
real(REAL_KIND) :: Radius, d, r2lim, r2, rad(3), rsite(3)
logical, allocatable :: occup(:,:,:)

blobcentre = DELTA_X*[NX/2,NY/2,NZ/2]
d = 2.0*Raverage
Radius = (3.0*initial_count/(4.0*PI))**(1./3.)	! scaled by /d
write(nflog,*) 'blobcentre, d: ',blobcentre,d,Radius
irad = Radius + 2
allocate(occup(-irad:irad,-irad:irad,-irad:irad))
occup = .false.
r2lim = 0.97*Radius*Radius
lastID = 0
ncells_mphase = 0
kcell = 0
do ix = -irad, irad
	do iy = -irad, irad
		do iz = -irad, irad
			site = (/ix,iy,iz/)
			rad = site
			r2 = dot_product(rad,rad)
			if (r2 > r2lim) cycle
			kcell = kcell + 1
			rsite = blobcentre + d*site
			call AddCell(kcell,rsite)
			occup(ix,iy,iz) = .true.
		enddo
	enddo
enddo
if (kcell > initial_count) then
	write(logmsg,*) 'Cell count already exceeds specified number: ',kcell,initial_count
	call logger(logmsg)
	ok = .false.
	return
endif
! Now add cells to make the count up to the specified initial_count
if (kcell < initial_count) then
	call AddBdryCells(kcell, occup, irad, d, blobcentre)
	kcell = initial_count
endif
deallocate(occup)
nlist = kcell
ncells = kcell
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(kcell, rsite)
integer :: kcell
real(REAL_KIND) :: rsite(3)
integer :: kpar = 0
real(REAL_KIND) :: r(3), c(3), R1, R2
type(cell_type), pointer :: cp
	
cp => cell_list(kcell)
cp%state = ALIVE
!			if (MITOSIS_MODE == TERMINAL_MITOSIS) then
cp%Iphase = .true.
cp%nspheres = 1
!			else
!				cp%Iphase = .false.
!				cp%nspheres = 2
!				ncells_mphase = ncells_mphase + 1
!			endif
			
!cp%radius(1) = Raverage
!cp%V = (4.*PI/3.)*Raverage**3						! need to randomise
cp%V_divide = Vdivide0
cp%V = (0.5 + 0.49*par_uni(kpar))*cp%V_divide
cp%radius(1) = (3*cp%V/(4*PI))**(1./3.)
cp%centre(:,1) = rsite
cp%site = rsite/DELTA_X + 1
cp%d = 0
cp%birthtime = 0
cp%growthrate = test_growthrate
!cp2%V_divide = get_divide_volume()
cp%d_divide = (3*cp%V_divide/PI)**(1./3.)
cp%mitosis = 0
call get_random_vector3(r)	! set initial axis direction
cp%d = 0.1*small_d
c = cp%centre(:,1)
cp%centre(:,1) = c + (cp%d/2)*r
cp%centre(:,2) = c - (cp%d/2)*r
cp%nbrs = 0
cp%Cex = Cextra(1,1,1,:)	! initially the concentrations are the same everywhere
cp%Cin = cp%Cex
end subroutine

!--------------------------------------------------------------------------------
! Add cells at the boundary to bring the total count from k up to initial_count
! (1) Make a list of all boundary sites (sites in contact with an OUTSIDE site)
! (2) Iteratively traverse the list to select the adjacent OUTSIDE site closest 
! to the centre.
!--------------------------------------------------------------------------------
subroutine AddBdryCells(klast, occup, irad, d, centre)
integer :: klast, irad
logical :: occup(-irad:irad,-irad:irad,-irad:irad)
real(REAL_KIND) :: d, centre(3)
integer :: ix, iy, iz
integer :: kcell, i, kb, site(3), nbsite(3), nbt, kbmin, imin
integer, allocatable :: sitelist(:,:)
real(REAL_KIND) :: r2, r2min, rad(3), rsite(3)

nbt = 0
do ix = -irad, irad
	do iy = -irad, irad
		do iz = -irad, irad
			if (.not.occup(ix,iy,iz)) cycle
			site = [ix,iy,iz]
			do i = 1,27
				if (i == 14) cycle
				nbsite = site + jumpvec(:,i)
				if (.not.occup(nbsite(1),nbsite(2),nbsite(3))) then
					nbt = nbt+1
					exit
				endif
			enddo
		enddo
	enddo
enddo

allocate(sitelist(3,nbt))

nbt = 0
do ix = -irad, irad
	do iy = -irad, irad
		do iz = -irad, irad
			if (.not.occup(ix,iy,iz)) cycle
			site = [ix,iy,iz]
			do i = 1,27
				if (i == 14) cycle
				nbsite = site + jumpvec(:,i)
				if (.not.occup(nbsite(1),nbsite(2),nbsite(3))) then
					nbt = nbt+1
					sitelist(:,nbt) = site
					exit
				endif
			enddo
		enddo
	enddo
enddo	

do kcell = klast+1,initial_count
	r2min = 1.0e10
	do kb = 1,nbt
		site = sitelist(:,kb)
		do i = 1,27
			if (i == 14) cycle
			nbsite = site + jumpvec(:,i)
			if (.not.occup(nbsite(1),nbsite(2),nbsite(3))) then
				rad = d*site
				r2 = dot_product(rad,rad)
				if (r2 < r2min) then
					kbmin = kb
					imin = i
					r2min = r2
				endif
			endif
		enddo
	enddo
	site = sitelist(:,kbmin) + jumpvec(:,imin)
	rsite = centre + d*site
	call AddCell(kcell,rsite)
	occup(site(1),site(2),site(3)) = .true.
enddo

deallocate(sitelist)
		
end subroutine

!-----------------------------------------------------------------------------------------
! Simulate through a full time step: DELTA_T
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: kcell, hour, kpar=0
real(REAL_KIND) :: dt
integer :: nchemo, i, k, nit, nt_diff, it_diff, ncells0
integer :: nshow = 100
integer :: Nhop
logical :: ok, done, changed
logical :: dbug

!Nhop = 10*(30/DELTA_T)
Nhop = 1
istep = istep + 1
dbug = .false.
if (Ncells == 0) then
	call logger('Ncells = 0')
!    res = 3
!    return
endif
!dt = DELTA_T/NT_CONC
! the idea is to accumulate time steps until DELTA_T is reached 
t_fmover = 0
nit = 0
done = .false.
call make_perm_index(ok)
if (.not.ok) then
	call logger('make_perm_index error')
	res = 4
	return
endif
!if (ncells > nshow) write(*,*) 'start moving'
do while (.not.done)
	nit = nit + 1
!	call mover(ok)
!	if (ncells > nshow) write(*,*) 'moving'
	call fmover(dt,done,ok)
	if (.not.ok) then
		call logger('mover error')
		res = 1
		return
	endif
	t_fmover = t_fmover + dt
	ncells0 = ncells
!	if (ncells >= nshow) write(*,'(a,2i6,3f8.1)') 'growing: istep, ndt, dt, delta_tmove: ',istep,ndt,dt,delta_tmove,t_fmover
	call grower(dt,changed,ok)
	if (.not.ok) then
		call logger('grower error')
		res = 2
		return
	endif
!	if (ncells > ncells0) then
	if (changed) then
		call make_perm_index(ok)
	endif
enddo

!if (mod(istep,Nhop) == 0) then
!	! determine cell death and tagging for death
!	call setup_grid_cells
!	call update_all_nbrlists
!	if (Ncells > 0) then
!		write(logmsg,'(a,5i8)') 'istep,Ncells,Nsteps,Nit,ndt: ',istep,Ncells,Nsteps,Nit,ndt
!		call logger(logmsg)
!	endif
!endif

! Reaction-diffusion system
! Assuming DELTA_T = 600 ...
if (ncells < 2000) then
	nt_diff = 1
elseif (ncells < 3000) then
	nt_diff = 2
elseif (ncells < 4000) then
	nt_diff = 3
elseif (ncells < 5000) then
	nt_diff = 4
elseif (ncells < 7000) then
	nt_diff = 5
elseif (ncells < 10000) then
	nt_diff = 6
else
	nt_diff = 7
endif
nt_diff = 1
dt = DELTA_T/nt_diff
do it_diff = 1,nt_diff
	call setup_grid_cells
!	if (ncells >= nshow) write(nflog,*) 'start setup_nbrlists'
	call update_all_nbrlists
	if (Ncells > 0) then
		write(logmsg,'(a,5i8)') 'istep,Ncells,Nsteps,Nit,ndt: ',istep,Ncells,Nsteps,Nit,ndt
		call logger(logmsg)
	endif
	call diff_solver(dt)
enddo
res = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine make_perm_index(ok)
logical :: ok
integer :: np, kcell, kpar=0

np = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	np = np + 1
	perm_index(np) = kcell
enddo
if (np /= ncells) then
	write(logmsg,*) 'Error: make_perm_index: np /= Ncells: ',np,ncells
	call logger(logmsg)
	ok = .false.
	return
endif
if (use_permute) then
	call permute(perm_index,np,kpar)
endif
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! In the Iphase, the cell grows as a sphere.
! When V > Vdivide, the cell enters Mphase
! In the Mphase, V is constant, d increases and R decreases (R1 = R2)
! The level of mitosis grows at a constant rate, and d is proportional to mitosis
! When d > 2R cell division is completed, and two cells replace one.
!-----------------------------------------------------------------------------------------
subroutine grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: k, kcell
type(cell_type), pointer :: cp
real(REAL_KIND) :: r(3), c(3), rad, tnow, d_desired, tgrowth, c_rate, r_mean
integer :: ndivide, divide_list(1000)

ok = .true.
changed = .false.
if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then
	tgrowth = divide_time_mean
else
	tgrowth = divide_time_mean - mitosis_duration
endif
c_rate = log(2.0)/tgrowth		! Note: to randomise divide time need to use random number, not mean!
r_mean = Vdivide0/(2*tgrowth)

ndivide = 0
tnow = istep*DELTA_T + t_fmover
do k = 1,ncells
	kcell = perm_index(k)
	cp => cell_list(kcell)
	if (cp%Cex(OXYGEN) < ANOXIA_THRESHOLD) then
		cp%state = DEAD
		Ncells = Ncells - 1
		changed = .true.
		write(*,*) 'Cell died: ',kcell
		cycle
	endif
!	write(*,*) 'test_growthrate, dVdt: ',test_growthrate*dt, cp%dVdt
!	if (cp%state == DEAD) cycle
	if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then	! always Mphase, growing and separating
		! need to set initial mitosis axis and d_divide at birth
!		cp%V = cp%V + cp%growthrate*DELTA_T
		call growcell(cp,dt,c_rate,r_mean)
		cp%mitosis = cp%V/cp%V_divide
		cp%d = cp%mitosis*cp%d_divide
		r = cp%centre(:,2) - cp%centre(:,1)
		r = r/sqrt(dot_product(r,r))	! axis direction unit vector
		c = (cp%centre(:,1) + cp%centre(:,2))/2
		cp%centre(:,1) = c - (cp%d/2)*r
		cp%centre(:,2) = c + (cp%d/2)*r
		call cubic_solver(cp%d,cp%V,rad)
		cp%radius = rad
		if (cp%d > 2*rad) then			! time for completion of cell division
			ndivide = ndivide + 1
			divide_list(ndivide) = kcell
		endif
		write(nflog,'(a,i4,3f8.3)') 'V,d,rad: ',kcell,cp%V,cp%d,rad
	elseif (MITOSIS_MODE == TERMINAL_MITOSIS) then
		if (cp%Iphase) then
!			cp%V = cp%V + cp%growthrate*dt
			call growcell(cp,dt,c_rate,r_mean)
			cp%radius(1) = (3*cp%V/(4*PI))**(1./3.)
			if (cp%V > cp%V_divide) then
				cp%Iphase = .false.
                cp%nspheres = 2
				ncells_mphase = ncells_mphase + 1
				call get_random_vector3(r)	! set initial axis direction
				cp%d = 0.1*small_d
				c = cp%centre(:,1)
				cp%centre(:,1) = c + (cp%d/2)*r
				cp%centre(:,2) = c - (cp%d/2)*r
				cp%mitosis = 0
				cp%t_start_mitosis = tnow
				cp%d_divide = 2.0**(2./3)*cp%radius(1)
!				write(*,*) 'start mitosis'
			endif
		else
			cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
			d_desired = max(cp%mitosis*cp%d_divide,small_d)
			r = cp%centre(:,2) - cp%centre(:,1)
			cp%d = sqrt(dot_product(r,r))
			r = r/cp%d	! axis direction
			c = (cp%centre(:,1) + cp%centre(:,2))/2
			cp%site = c/DELTA_X + 1
!			cp%centre(1,:) = c - (cp%d/2)*r		! For fmover we do not need to set centre positions
!			cp%centre(2,:) = c + (cp%d/2)*r
			call cubic_solver(d_desired,cp%V,rad)
!			write(*,*) 'mitosis,d,d_divide,V,rad: ',cp%mitosis,cp%d,cp%d_divide,cp%V,rad
	if (rad > 10.e-4) then
		write(*,*) 'grower: big rad, d: ',kcell,rad,cp%d
		stop
	endif
			cp%radius = rad
			if (cp%d > 2*rad) then			! time for completion of cell division
				ndivide = ndivide + 1
				divide_list(ndivide) = kcell
			endif
		endif
	endif
	if (cp%radius(1) > 10.e-4) then
		write(*,*) 'grower: big radius: ',kcell,cp%radius
		stop
	endif
enddo
!if (ncells > 3000 .and. ndivide > 0) write(*,*) 'dividing'
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	call divider(kcell)
enddo
!if (ndivide > 0) then
!	if (ncells > 3000) write(*,*) 'setup_nbrlists'
!    call setup_nbrlists
!endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine growcell(cp, dt, c_rate, r_mean)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt, c_rate, r_mean
real(REAL_KIND) :: Cin_0(NCONST), Cex_0(NCONST)		! at some point NCONST -> MAX_CHEMO
real(REAL_KIND) :: dVdt,  Vin_0, dV, metab
integer :: C_option = 1	! we must use this

Cin_0 = cp%Cin
metab = O2_metab(Cin_0(OXYGEN))	! Note that currently growth depends only on O2
if (use_V_dependence) then
	dVdt = c_rate*metab*cp%V/(Vdivide0/2)
else
	dVdt = r_mean*metab
endif
!if (suppress_growth) then	! for checking solvers
!	dVdt = 0
!endif
!site = cell_list(kcell)%site
Cex_0 = cp%Cex
cp%dVdt = dVdt
cp%growthrate = dVdt
Vin_0 = cp%V
dV = dVdt*dt
cp%V = Vin_0 + dV
if (C_option == 1) then
	! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
	cp%Cin = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
!	Cextra = (Vex_0*Cex_0 - dV*Cex_0)/(Vex_0 - dV)	! = Cex_0
elseif (C_option == 2) then
	! Calculation based on change in volumes without mass transfer of constituents
	cp%Cin = Vin_0*Cin_0/(Vin_0 + dV)
!	Cextra = Vex_0*Cex_0/(Vex_0 - dV)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! A single cell is replaced by two.
!-----------------------------------------------------------------------------------------
subroutine divider(kcell1)
integer :: kcell1
integer :: kcell2, nbrs0
real(REAL_KIND) :: r(3), c(3)
type(cell_type), pointer :: cp1, cp2

!write(*,*) 'divider:'
!write(logmsg,*) 'divider: ',kcell1
!call logger(logmsg)
tnow = istep*DELTA_T
cp1 => cell_list(kcell1)
nlist = nlist + 1
ncells = ncells + 1
ncells_mphase = ncells_mphase - 1
kcell2 = nlist
!write(*,*) 'divider: ',kcell1,kcell2
cp2 => cell_list(kcell2)
cp1%state = ALIVE
cp1%V = cp1%V/2
cp1%site = cp1%centre(:,1)/DELTA_X + 1
cp1%d = 0
cp1%birthtime = tnow
cp1%V_divide = get_divide_volume()
cp1%d_divide = (3*cp1%V_divide/PI)**(1./3.)
cp1%mitosis = 0
nbrs0 = cp1%nbrs
cp1%nbrs = nbrs0 + 1
cp1%nbrlist(cp1%nbrs)%indx = kcell2
cp1%nbrlist(cp1%nbrs)%contact = .false.
cp1%nbrlist(cp1%nbrs)%contact(1,1) = .true.
cp2%state = ALIVE
cp2%V = cp1%V
cp2%radius(1) = cp1%radius(2)
cp2%centre(:,1) = cp1%centre(:,2)
cp2%site = cp2%centre(:,1)/DELTA_X + 1
cp2%d = 0
cp2%birthtime = tnow
cp2%growthrate = cp1%growthrate		!!!!!!!!!!!!! for testing, keep this constant !!!!!!!!!!!!!
cp2%V_divide = get_divide_volume()
cp2%d_divide = (3*cp2%V_divide/PI)**(1./3.)
cp2%mitosis = 0

cp2%Cin = cp1%Cin
cp2%Cex = cp1%Cex

!allocate(cp2%nbrlist(MAX_NBRS))
! Need to set up neighbour lists for both cp1 and cp2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< DO THIS
! Use neighbour lists of neighbours
! First simple test:
cp2%nbrlist(1:nbrs0) = cp1%nbrlist(1:nbrs0)
cp2%nbrs = nbrs0 + 1
cp2%nbrlist(cp2%nbrs)%indx = kcell1
cp2%nbrlist(cp2%nbrs)%contact = .false.
cp2%nbrlist(cp2%nbrs)%contact(1,1) = .true.
if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then
	cp1%Iphase = .false.
	cp2%Iphase = .false.
	call get_random_vector3(r)	! set initial axis direction
	cp1%d = 0.1*small_d
	c = cp1%centre(:,1)
	cp1%centre(:,1) = c + (cp1%d/2)*r
	cp1%centre(:,2) = c - (cp1%d/2)*r
	cp2%d = 0.1*small_d
	c = cp2%centre(:,1)
	cp2%centre(:,1) = c + (cp2%d/2)*r
	cp2%centre(:,2) = c - (cp2%d/2)*r
else
	cp1%Iphase = .true.
    cp1%nspheres = 1
	cp2%Iphase = .true.
    cp2%nspheres = 1
endif
call update_nbrlist(kcell1)
call update_nbrlist(kcell2)
! Note: any cell that has kcell1 in its nbrlist now needs to add kcell2 to the nbrlist.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function get_divide_volume() result(vol)
real(REAL_KIND) :: vol
integer :: kpar = 0
real(REAL_KIND) :: U

U = par_uni(kpar)
vol = Vdivide0 + dVdivide*(2*U-1)
end function

!-----------------------------------------------------------------------------------------
! Cap volume: http://mathworld.wolfram.com/SphericalCap.html
! http://en.wikipedia.org/wiki/Cubic_function
! We find x1, the first real root
! dc = centre separation distance (um)
! V = constant volume (um3)
!-----------------------------------------------------------------------------------------
subroutine cubic_solver(dc,V,R)
real(REAL_KIND) :: dc, V, R
real(REAL_KIND) :: a, b, c, d
real(REAL_KIND) :: DD0, DD1, CC, u1, x1

a = 2.0
b = 1.5*dc
c = 0
d = -(dc**3/8 + 3*V/(2*PI))

DD0 = b**2 - 3*a*c
DD1 = 2*b**3 - 9*a*b*c + 27*a**2*d
CC = (DD1 + sqrt(DD1**2 - 4*DD0**3))/2
if (CC < 0) then
	CC = -(-CC)**(1.d0/3)
else
	CC = CC**(1.d0/3)
endif

u1 = 1
R = -(1/(3*a))*(b + u1*CC + DD0/(u1*CC))

!write(*,*) 'R: ',R
!write(*,*) 'cubic: ',a*R**3 + b*R**2 + c*R + d
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

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim, NY_dim, NZ_dim, nsteps_dim, deltat, maxchemo, nextra, cused, dfraction, deltax) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,nsteps_dim, maxchemo, nextra
real(c_double) :: deltat, dfraction, deltax
logical(c_bool) :: cused(*)
integer :: ichemo
integer :: N_EXTRA=1

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
nsteps_dim = nsteps
deltat = DELTA_T
deltax = DELTA_X
maxchemo = MAX_CHEMO
nextra = N_EXTRA
do ichemo = 1,MAX_CHEMO
	cused(ichemo+1) = .true.	!chemo(ichemo)%used
enddo
cused(1) = .true.			! CFSE
cused(MAX_CHEMO+2) = .true.	! Growth rate
dfraction = 1.0	!2*cell_radius/DELTA_X
end subroutine

!--------------------------------------------------------------------------------
! TC = tumour cell
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nTC_list
type(celldata_type) :: TC_list(*)
integer :: kcell, nspheres, is
type(cell_type), pointer :: cp

!write(nflog,*) 'get_scene'
nTC_list = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	if (cp%Iphase) then
		nspheres = 1
	else
		nspheres = 2
	endif
	do is = 1,nspheres
		nTC_list = nTC_list + 1
		TC_list(nTC_list)%tag = nTC_list
		TC_list(nTC_list)%radius = 10000*cp%radius(is)		! cm -> um
		TC_list(nTC_list)%centre = 10000*cp%centre(:,is)	! cm -> um
		TC_list(nTC_list)%celltype = 1
		TC_list(nTC_list)%highlight = 0
!		write(nflog,'(2i6,4e12.3)') nTC_list,TC_list(nTC_list)%tag,TC_list(nTC_list)%radius,TC_list(nTC_list)%centre
	enddo
enddo
!write(nflog,*) 'nTC_list: ',nTC_list
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_DLL_build_version(version_array,array_len) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_dll_build_version
use, intrinsic :: iso_c_binding
character(c_char) :: version_array(12)
integer(c_int) :: array_len
integer :: k

dll_version = DLL_BUILD_VERSION
gui_version = GUI_BUILD_VERSION
!write(nflog,*) 'get_DLL_build_version: ',dll_version
do k = 1,12
	version_array(k) = dll_version(k:k)
!	write(nflog,'(i2,a,a)') k,' ',version_array(k)
	if (version_array(k) == ' ') then
		version_array(k) = char(0)
		array_len = k
		exit
	endif
enddo
!write(nflog,*) 'array_len: ',array_len
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen,centre) BIND(C)
!subroutine Execute() BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
real(c_double) :: centre(*)
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
centre(1:3) = 10000*blobcentre
write(nflog,*) 'blobcentre: ',blobcentre
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
!call logger(logmsg)

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

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testresults
integer :: k1, k2
real(REAL_KIND) :: R1, R2, c1(3), c2(3), v(3), d2, d

write(*,*) 'testresults'
do k1 = 1,nlist
	if (cell_list(k1)%state == DEAD) cycle
	R1 = cell_list(k1)%radius(1)
	c1 = cell_list(k1)%centre(:,1)
	write(*,'(a,i4,3f6.2)') 'k1, c: ',k1,c1
	do k2 = k1+1, nlist
		if (cell_list(k2)%state == DEAD) cycle
		R2 = cell_list(k2)%radius(1)
		c2 = cell_list(k2)%centre(:,1)
		v = c2 - c1
		d2 = dot_product(v,v)
		d = sqrt(d2)
		write(*,'(2i4,f8.2)') k1,k2,d
	enddo
enddo

end subroutine


end module