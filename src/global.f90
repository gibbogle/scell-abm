module global
use real_kind_mod
use omp_lib
use par_zig_mod
use winsock
use, intrinsic :: ISO_C_BINDING

implicit none

#include "itsol_interface.f90"

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting) 
integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4
integer, parameter :: MAX_CELLTYPES = 4
integer, parameter :: MAX_DRUGTYPES = 2

integer, parameter :: X_AXIS = 1
integer, parameter :: Y_AXIS = 2
integer, parameter :: Z_AXIS = 3

integer, parameter :: TPZ_CLASS = 1
integer, parameter :: DNB_CLASS = 2
integer, parameter :: DRUG_EVENT = 1
integer, parameter :: RADIATION_EVENT = 2
integer, parameter :: MEDIUM_EVENT = 3


!real(REAL_KIND), parameter :: dx = 3.0e-3			! (cm) = 30 um
integer, parameter :: NRF = 4
!real(REAL_KIND), parameter :: dxb = NRF*dx			! (cm) = 120 um

real(REAL_KIND), parameter :: Vsite_cm3 = 2.0e-9	! from spheroid-abm, to scale Kin, Kout
real(REAL_KIND), parameter :: PI = 4*atan(1.d0)
real(REAL_KIND), parameter :: BIG = 1.0d10
real(REAL_KIND), parameter :: um3_cm3 = 1.0e-12
real(REAL_KIND), parameter :: small_d = 0.1e-4		! 0.1 um -> cm
real(REAL_KIND), parameter :: CFSE_std = 0.05
real(REAL_KIND), parameter :: Vcell_pL = 1.0e09

integer, parameter :: ALIVE = 0
integer, parameter :: DEAD = 1
integer, parameter :: TERMINAL_MITOSIS   = 1
integer, parameter :: CONTINUOUS_MITOSIS = 2
!integer, parameter :: MAX_CHEMO = 5
integer, parameter :: CFSE = 0
integer, parameter :: OXYGEN = 1
integer, parameter :: GLUCOSE = 2
integer, parameter :: TRACER = 3
integer, parameter :: DRUG_A = 4
integer, parameter :: TPZ_DRUG = DRUG_A
integer, parameter :: TPZ_DRUG_METAB_1 = TPZ_DRUG + 1
integer, parameter :: TPZ_DRUG_METAB_2 = TPZ_DRUG + 2
integer, parameter :: DRUG_B = DRUG_A + 3
integer, parameter :: DNB_DRUG = DRUG_B
integer, parameter :: DNB_DRUG_METAB_1 = DNB_DRUG + 1
integer, parameter :: DNB_DRUG_METAB_2 = DNB_DRUG + 2
integer, parameter :: MAX_CHEMO = DRUG_B + 2
integer, parameter :: GROWTH_RATE = MAX_CHEMO + 1	! (not used here, used in the GUI)
integer, parameter :: CELL_VOLUME = MAX_CHEMO + 2
integer, parameter :: O2_BY_VOL = MAX_CHEMO + 3

integer, parameter :: N_EXTRA = O2_BY_VOL - MAX_CHEMO + 1	! = 4 = total # of variables - MAX_CHEMO
integer, parameter :: NCONST = MAX_CHEMO

integer, parameter :: MITOSIS_MODE = TERMINAL_MITOSIS

type neighbour_type
	integer :: indx
!	logical :: incontact
	logical*1 :: contact(2,2)
end type

type :: grid_type
	integer :: nc
	integer :: cell(100)
end type

type dist_type
	integer :: class
	real(REAL_KIND) :: p1, p2, p3
end type

type :: stencil_type
	integer :: nindx
	integer :: indx(3,6)
end type

type, bind(C) :: celldata_type
	integer(c_int) :: tag
	real(c_double) :: radius
	real(c_double) :: centre(3)
	integer(c_int) :: celltype
	integer(c_int) :: status
end type

type, bind(C) :: fielddata_type
    integer(c_int) :: NX, NY, NZ, NCONST
    real(c_double) :: DX
    type(c_ptr) :: Conc_ptr   ! Caverage(NX,NY,NZ,NCONST)
    integer(c_int) :: ncells
    type(c_ptr) :: cell_ptr
end type


type cell_type
!	real(REAL_KIND) :: V_n          ! "normal" volume
!	integer :: varindex
	integer :: ID
	integer :: celltype
	integer :: state
	logical :: Iphase
    integer :: nspheres             ! =1 for Iphase, =2 for Mphase
    integer :: site(3)				! this is the gridcell = (i,j,k) then cell is in: i-i+1, j-j+1, k-k+1
	real(REAL_KIND) :: V			! actual volume (um3) -> cm3
	real(REAL_KIND) :: dVdt
	real(REAL_KIND) :: radius(2)	! sphere radii (um) -> cm
	real(REAL_KIND) :: centre(3,2)  ! sphere centre positions
	real(REAL_KIND) :: d			! centre separation distance (um) -> cm
	real(REAL_KIND) :: birthtime
	real(REAL_KIND) :: V_divide
	real(REAL_KIND) :: t_hypoxic
	real(REAL_KIND) :: t_anoxia_die
	real(REAL_KIND) :: t_start_mitosis
	real(REAL_KIND) :: mitosis		! level of mitosis (0 - 1)
	real(REAL_KIND) :: d_divide		! centre separation distance at the end of mitosis
	integer :: cnr(3,8)
	real(REAL_KIND) :: wt(8)
	
	integer :: nbrs
	type(neighbour_type) :: nbrlist(100)
	real(REAL_KIND) :: Cin(NCONST)
	real(REAL_KIND) :: Cex(NCONST)
	real(REAL_KIND) :: dCdt(NCONST)
	real(REAL_KIND) :: dMdt(NCONST)	! mumol/s
	real(REAL_KIND) :: CFSE
	logical :: radiation_tag, anoxia_tag	!, drugA_tag, drugB_tag
	logical :: drug_tag(MAX_DRUGTYPES)
	real(REAL_KIND) :: p_death
end type

type XYZ_type
    real(REAL_KIND) :: x, y, z
end type

type treatment_type
	integer :: ichemo
	integer :: n
!	character*(16) :: name
	real(REAL_KIND), allocatable :: tstart(:)
	real(REAL_KIND), allocatable :: tend(:)
	real(REAL_KIND), allocatable :: conc(:)
	real(REAL_KIND), allocatable :: dose(:)
	logical, allocatable :: started(:)
	logical, allocatable :: ended(:)
end type

type event_type
	integer :: etype
	real(REAL_KIND) :: time
	integer :: idrug			! DRUG
	integer :: ichemo			! DRUG CHEMO INDEX
	real(REAL_KIND) :: volume	! DRUG MEDIUM
	real(REAL_KIND) :: conc		! DRUG
	real(REAL_KIND) :: dose		! RADIATION
	logical :: done
end type	

type drug_type
	character*(3)   :: classname
	integer         :: drugclass
	character*(16)  :: name
	integer         :: nmetabolites
	logical         :: use_metabolites
	real(REAL_KIND) :: diff_coef(0:2)
	real(REAL_KIND) :: medium_diff_coef(0:2)
	real(REAL_KIND) :: membrane_diff_in(0:2)
	real(REAL_KIND) :: membrane_diff_out(0:2)
	real(REAL_KIND) :: halflife(0:2)
	logical         :: kills(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kmet0(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: C2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: KO2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Vmax(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Klesion(MAX_CELLTYPES,0:2)
	integer         :: kill_model(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kd(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_drug(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_duration(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_fraction(MAX_CELLTYPES,0:2)
	logical         :: sensitises(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_max(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: SER_KO2(MAX_CELLTYPES,0:2)
end type
type SN30K_type
	integer :: nmetabolites
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: medium_diff_coef
	real(REAL_KIND) :: membrane_diff
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: metabolite_halflife
	real(REAL_KIND) :: Kmet0(MAX_CELLTYPES)
!	real(REAL_KIND) :: C1(MAX_CELLTYPES)
	real(REAL_KIND) :: C2(MAX_CELLTYPES)
	real(REAL_KIND) :: KO2(MAX_CELLTYPES)
!	real(REAL_KIND) :: gamma(MAX_CELLTYPES)
	real(REAL_KIND) :: Klesion(MAX_CELLTYPES)
	integer         :: kill_model(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_O2(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_drug(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_duration(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_fraction(MAX_CELLTYPES)
	real(REAL_KIND) :: Kd(MAX_CELLTYPES)
end type

type TPZ_type
	character*(16) :: name
	integer :: nmetabolites
	logical :: use_metabolites
	real(REAL_KIND) :: diff_coef(0:2)
	real(REAL_KIND) :: medium_diff_coef(0:2)
	real(REAL_KIND) :: membrane_diff_in(0:2)
	real(REAL_KIND) :: membrane_diff_out(0:2)
	real(REAL_KIND) :: halflife(0:2)
	real(REAL_KIND) :: Kmet0(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: C2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: KO2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Vmax(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Klesion(MAX_CELLTYPES,0:2)
	integer         :: kill_model(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_O2(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_drug(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_duration(MAX_CELLTYPES)
	real(REAL_KIND) :: kill_fraction(MAX_CELLTYPES)
	real(REAL_KIND) :: Kd(MAX_CELLTYPES)
end type


type DNB_type
	character*(16) :: name
	integer :: nmetabolites
	logical :: use_metabolites
	real(REAL_KIND) :: diff_coef(0:2)
	real(REAL_KIND) :: medium_diff_coef(0:2)
	real(REAL_KIND) :: membrane_diff_in(0:2)
	real(REAL_KIND) :: membrane_diff_out(0:2)
	real(REAL_KIND) :: halflife(0:2)
	real(REAL_KIND) :: Kmet0(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: C2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: KO2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Vmax(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Km(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Klesion(MAX_CELLTYPES,0:2)
	integer         :: kill_model(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_O2(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_drug(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_duration(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: kill_fraction(MAX_CELLTYPES,0:2)
	real(REAL_KIND) :: Kd(MAX_CELLTYPES,0:2)
end type

type LQ_type
	real(REAL_KIND) :: OER_am, OER_bm
	real(REAL_KIND) :: alpha_H, beta_H
	real(REAL_KIND) :: K_ms
	real(REAL_KIND) :: death_prob
end type

integer, parameter :: nflog=10, nfin=11, nfout=12, nfres=13, nfcell=14
integer, parameter :: MAX_NLIST = 100000
integer, parameter :: MAX_NBRS = 100
integer, parameter :: ndt_max = 30
real(REAL_KIND), parameter :: Raverage = 0.64e-3	! as in spheroid-abm.  was: 5.0e-4*1.5**(1./3)		! average cells radius (um) -> cm

character*(128) :: inputfile
character*(128) :: outputfile
character*(12) :: dll_version, dll_run_version
character*(12) :: gui_version, gui_run_version

integer :: Mnodes, ncpu_input, ncells, ncells_mphase, nlist, nsteps, nevents
integer :: Ndrugs_used
integer :: NX, NY, NZ, NXB, NYB, NZB
integer :: Ndim(3)
integer :: nt_saveprofiledata, it_saveprofiledata
integer :: NT_CONC, NT_GUI_OUT, initial_count, ntries, Ncelltypes, Ncells_type(MAX_CELLTYPES)
integer :: istep, ndt, ndtotal, ndtotal_last, ichemo_curr
integer :: seed(2)
integer :: jumpvec(3,27)
real(REAL_KIND) :: DELTA_X, DELTA_T, tnow, t_simulation, dt_saveprofiledata
real(REAL_KIND) :: blobcentre(3)	! blob centre
real(REAL_KIND) :: epsilon, es_e, sqr_es_e, shift, Dfactor
real(REAL_KIND) :: alpha_v, k_detach
real(REAL_KIND) :: dr_mitosis, mitosis_hours, mitosis_duration
real(REAL_KIND) :: test_growthrate
real(REAL_KIND) :: Vdivide0, dVdivide, Rdivide0, MM_THRESHOLD, medium_volume0, total_volume
real(REAL_KIND) :: t_anoxic_limit, anoxia_death_delay, ANOXIA_THRESHOLD, anoxia_tag_hours, anoxia_death_hours
real(REAL_KIND) :: divide_time_median(MAX_CELLTYPES), divide_time_shape(MAX_CELLTYPES), divide_time_mean(MAX_CELLTYPES), celltype_fraction(MAX_CELLTYPES)
type(dist_type) :: divide_dist(MAX_CELLTYPES)
real(REAL_KIND) :: d_nbr_limit

logical :: use_dropper
integer :: Ndrop
real(REAL_KIND) :: alpha_shape, beta_shape	! squashed sphere shape parameters
real(REAL_KIND) :: adrop, bdrop, cdrop		! drop shape transformation parameters
integer :: zmin     						! drop lower bound at drop time = lower limit of blob thereafter
logical :: is_dropped

integer :: Nradiation_tag(MAX_CELLTYPES), Nanoxia_tag(MAX_CELLTYPES)
!integer :: NdrugA_tag(MAX_CELLTYPES), NdrugB_tag(MAX_CELLTYPES), Nradiation_tag(MAX_CELLTYPES), Nanoxia_tag(MAX_CELLTYPES)
integer :: Ndrug_tag(MAX_DRUGTYPES,MAX_CELLTYPES)
!integer :: NdrugA_dead(MAX_CELLTYPES), NdrugB_dead(MAX_CELLTYPES)
integer :: Nradiation_dead(MAX_CELLTYPES), Nanoxia_dead(MAX_CELLTYPES)
integer :: Ndrug_dead(MAX_DRUGTYPES,MAX_CELLTYPES)
real(REAL_KIND) :: O2cutoff(3)
real(REAL_KIND) :: growthcutoff(3)


! From react_diff
real(REAL_KIND) :: dxf, dxb, dx3, dxb3, Rcell, Vcell
real(REAL_KIND), allocatable :: Vin(:,:,:), dVindt(:,:,:)
integer, allocatable :: ngcells(:,:,:)
real(REAL_KIND), allocatable :: afull(:,:)

real(REAL_KIND) :: a_separation, kdrag, frandom
real(REAL_KIND) :: a_force, b_force, c_force, x0_force, x1_force, xcross1_force, xcross2_force
real(REAL_KIND) :: t_fmover, delta_tmove, dt_min, delta_min, delta_max

real(REAL_KIND) :: C_O2_bdry
!real(REAL_KIND) :: Kmemb(NCONST)
real(REAL_KIND) :: total_dMdt

integer :: show_progeny
logical :: use_V_dependence, randomise_initial_volume
logical :: simulation_start, par_zig_init, initialized

TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send

character*(128) :: logfile
character*(2048) :: logmsg
character*(128) :: profiledatafilebase

type(cell_type), allocatable, target :: cell_list(:)
type(grid_type), allocatable, target :: grid(:,:,:)
integer, allocatable :: perm_index(:)
type(stencil_type), allocatable :: stencil(:,:,:)
real(REAL_KIND), allocatable, target :: Cextra_all(:,:,:,:)
real(REAL_KIND), allocatable, target :: Caverage(:,:,:,:)
real(REAL_KIND), allocatable, target :: Cflux(:,:,:,:)
real(REAL_KIND), allocatable, target :: Cflux_prev(:,:,:,:)

type(treatment_type), allocatable :: protocol(:)
type(event_type), allocatable :: event(:)


type(TPZ_type) :: TPZ
type(DNB_type) :: DNB
type(drug_type), allocatable, target :: drug(:)
type(LQ_type) :: LQ(MAX_CELLTYPES)

logical :: use_events = .true.
logical :: use_radiation, use_treatment
logical :: saveprofiledata
logical :: use_death = .true.
logical :: use_extracellular_O2 = .false.
logical :: use_migration = .false.
logical :: suppress_growth = .false.
logical :: use_hysteresis = .false.
logical :: use_permute = .false.
logical :: use_SS = .false.
logical :: use_integration = .true.
logical :: use_packer = .true.
logical :: dbug = .false.

integer :: kcell_now

!dec$ attributes dllexport :: nsteps, DELTA_T

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: logfile_isopen
character*(1) :: LF = char(94)

error = 0
inquire(unit=nflog,OPENED=logfile_isopen)
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    elseif (logfile_isopen) then
        write(nflog,*) trim(msg)
    else
        write(99,*) trim(msg)
    endif
else
	write(*,*) trim(msg)
endif
if (logfile_isopen) then
	write(nflog,'(a,a)') 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!-----------------------------------------------------------------------------------------
! Estimate blob radius by determining the range in the three axis directions
! This is clearly inaccurate when the blob shape deviates significantly from a sphere,
! and especially when cell killing causes the blob to break up.
!-----------------------------------------------------------------------------------------
function getRadius() result(R)
real(REAL_KIND) :: R
real(REAL_KIND) :: cntr(3), rng(3), diam

call getBlobCentreRange(cntr,rng, R)
end function

!-----------------------------------------------------------------------------------------
! Estimate blob centre, range in each axis direction, and radius.
!-----------------------------------------------------------------------------------------
subroutine getBlobCentreRange(cntr,rng,radius)
real(REAL_KIND) :: cntr(3), rng(3), radius
real(REAL_KIND) :: rmin(3), rmax(3), diam
integer :: kcell, k
type(cell_type), pointer :: cp

rmin = 1.0e10
rmax = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	do k = 1,cp%nspheres
		rmin = min(rmin,cp%centre(:,k)-cp%radius(k))
		rmax = max(rmax,cp%centre(:,k)+cp%radius(k))
	enddo
enddo
cntr = (rmax + rmin)/2
rng = rmax - rmin
diam = ((rng(1)**3 + rng(2)**3 + rng(3)**3)/3)**(1./3.)
radius = diam/2
end subroutine

!-----------------------------------------------------------------------------------------
! Find the cells that are at the min and max locations in the X, Y and Z directions
!-----------------------------------------------------------------------------------------
subroutine getBdryCells(bcells)
integer :: bcells(3,2)
real(REAL_KIND) :: rmin(3), rmax(3)
integer :: kcell, k, i
type(cell_type), pointer :: cp

rmin = 1.0e10
rmax = -1.0e10
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	do k = 1,cp%nspheres
		do i = 1,3
			if (cp%centre(i,k) < rmin(i)) then
				rmin(i) = cp%centre(i,k)
				bcells(i,1) = kcell
			endif
			if (cp%centre(i,k) > rmax(i)) then
				rmax(i) = cp%centre(i,k)
				bcells(i,2) = kcell
			endif
		enddo
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
!if (Mnodes == 1) return
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif

call logger('did omp_initialisation')
!call test_omp1

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RngInitialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
write(logmsg,*) 'npar = ',npar,seed
call logger(logmsg)
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
if (R == -2147483648) R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))
end function

!-----------------------------------------------------------------------------------------
! This needs to be given radial symmetry
!-----------------------------------------------------------------------------------------
subroutine get_random_dr(dr)
real(REAL_KIND) :: dr(3)
integer :: kpar=0

dr(1) = 2*(par_uni(kpar) - 0.5)
dr(2) = 2*(par_uni(kpar) - 0.5)
dr(3) = 2*(par_uni(kpar) - 0.5)
end subroutine

!-----------------------------------------------------------------------------------------
! Returns a unit vector with random 3D direction
!-----------------------------------------------------------------------------------------
subroutine get_random_vector3(v)
real(REAL_KIND) :: v(3)
real(REAL_KIND) :: R1, R2, s, a
integer :: kpar=0

R1 = par_uni(kpar)
R2 = par_uni(kpar)
s = sqrt(R2*(1-R2))
a = 2*PI*R1
v(1) = 2*cos(a)*s
v(2) = 2*sin(a)*s
v(3) = 1 - 2*R2
end subroutine

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(REAL_KIND) :: p(:)
integer :: k
real(REAL_KIND) :: R, psum

R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
end function
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!subroutine get_random_drot(axis,drot)
!real(REAL_KIND) :: drot
!integer :: axis
!integer :: kpar=0
!
!axis = random_int(1,3,kpar)
!drot = 2*(par_uni(kpar) - 0.5)
!end subroutine
!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!----------------------------------------------------------------------------------------
! R = A x B 
!----------------------------------------------------------------------------------------
subroutine cross_product(A, B, R)
real(REAL_KIND) :: A(3), B(3), R(3)

R(1) = A(2)*B(3) - A(3)*B(2)
R(2) = A(3)*B(1) - A(1)*B(3)
R(3) = A(1)*B(2) - A(2)*B(1)
end subroutine

!-----------------------------------------------------------------------------------------------------
! Rotate v0 about unit vector vN through angle to get v
!-----------------------------------------------------------------------------------------------------
subroutine Rotate(v0, vN, angle, v)
real(REAL_KIND) :: v0(3), vN(3), angle, v(3)
real(REAL_KIND) :: cosa, sina, d
type(XYZ_type) :: q1, q2, u

cosa = cos(angle)
sina = sin(angle)

q1%x = v0(1)
q1%y = v0(2)
q1%z = v0(3)
u%x = vN(1)
u%y = vN(2)
u%z = vN(3)

d = sqrt(u%y*u%y + u%z*u%z)

! Step 2 
if (d /= 0) then
    q2%x = q1%x
    q2%y = q1%y * u%z / d - q1%z * u%y / d
    q2%z = q1%y * u%y / d + q1%z * u%z / d
else
  q2 = q1
endif

! Step 3
q1%x = q2%x * d - q2%z * u%x
q1%y = q2%y
q1%z = q2%x * u%x + q2%z * d

! Step 4
q2%x = q1%x * cosa - q1%y * sina
q2%y = q1%x * sina + q1%y * cosa
q2%z = q1%z

! Inverse of step 3
q1%x =   q2%x * d + q2%z * u%x
q1%y =   q2%y
q1%z = - q2%x * u%x + q2%z * d

! Inverse of step 2 
if (d /= 0) then
    v(1) =   q1%x
    v(2) =   q1%y * u%z / d + q1%z * u%y / d
    v(3) = - q1%y * u%y / d + q1%z * u%z / d
else
    v(1) = q1%x
    v(2) = q1%y
    v(3) = q1%z
endif
end subroutine

!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
subroutine test_rotate
real(REAL_KIND) :: v0(3), vN(3), angle, v(3), r(3), F(3), M(3), mamp

v0(1) = 1   ! cell axis
v0(2) = 0
v0(3) = 0
F(1) = 0    ! force
F(2) = 0
F(3) = -1
r(1) = 3  ! offset of point of contact from centre
r(2) = 0
r(3) = 0
call cross_product(r,F,M)
mamp = sqrt(dot_product(M,M))   ! moment amplitude
vN = M/mamp     ! unit vector of moment axis
angle = DELTA_T*0.01*mamp
call rotate(v0,vN,angle,v)
write(nflog,'(a,3f8.4)') 'rotate: ',v0
write(nflog,'(a,3f8.4)') 'about:  ',vN
write(nflog,'(a,3f8.4)') 'angle:  ',angle
write(nflog,'(a,3f8.4)') 'yields: ',v

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine make_jumpvec
integer :: k,ix,iy,iz

k = 0
do ix = -1,1
	do iy = -1,1
		do iz = -1,1
			k = k+1
			jumpvec(:,k) = (/ix,iy,iz/)
		enddo
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!--------------------------------------------------------------------------------
SUBROUTINE qsort(a, n, t)
IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL(REAL_KIND), INTENT(INOUT)    :: a(n)
INTEGER, INTENT(INOUT) :: t(n)

!     Local Variables

INTEGER                :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL(REAL_KIND)        :: w, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 CONTINUE
l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20 CONTINUE
i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

DO
  DO
    IF (a(i).LT.x) THEN                ! Search from lower end
      i = i + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  DO
    IF (x.LT.a(j)) THEN                ! Search from upper end
      j = j - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (i.LE.j) THEN                     ! Swap positions i & j
    w = a(i)
    ww = t(i)
    a(i) = a(j)
    t(i) = t(j)
    a(j) = w
    t(j) = ww
    i = i + 1
    j = j - 1
    IF (i.GT.j) EXIT
  ELSE
    EXIT
  END IF
END DO

IF (j-l.GE.r-i) THEN
  IF (l.LT.j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i.LT.r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF

IF (l.LT.r) GO TO 20
IF (s.NE.0) GO TO 10

RETURN
END SUBROUTINE qsort

end module