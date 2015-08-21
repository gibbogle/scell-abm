! To test solution of reaction-diffusion on a rectangular grid, using fgmres with ITSOL ILUK preconditioning
! This version uses two rectangular grids, coarse and fine.
! The fine grid is large enough to contain the spheroid blob completely, and has grid resolution
! such that a grid-cell can hold a few tumour cells, e.g. dxf = 20 -> 30 um
! The fine grid is embedded in a coarse grid, which has a resolution dxb that is a multiple of dxf,
! e.g. dxb = 4*dxf. Coarse grid points coincide with points in the fine grid, and in particular
! the points on the boundary of the fine grid coincides with a line of the course grid.
! The size of the coarse grid is set to match the medium volume, if not the shape.
! 
! The idea behind the solution method is that solving on the coarse grid provides the boundary
! values for solution on the fine grid.
!
! Quasi-steady-state method:
! -------------------------
! Time-dependence is handled in the coarse solution, using the IMEX 2-SBDF scheme.  The most recent 
! rates of uptake-secretion of constituents by cells are used to compute grid flux F().
! After the Cave values on the coarse grid have been determined, the boundary concentrations on
! the fine grid are estimated by interpolation, and the steady-state solution on the fine grid
! is computed.  It is this solution that provides the flux values for the next coarse grid time step.
! 
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!
! Parallel issues:
! https://software.intel.com/en-us/articles/threading-fortran-applications-for-parallel-performance-on-multi-core-systems/
!
module react_diff

use real_kind_mod
use global
use omp_lib
use sparse_map
use par_zig_mod
use chemokine
use continuum

use, intrinsic :: iso_c_binding

implicit none

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine setup_react_diff
integer :: ix, iy, iz, ic, maxnz, ichemo
real(REAL_KIND) :: C0
character*(10) :: emapfile, bmapfile
real(REAL_KIND), pointer :: Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
logical :: zero
logical :: ok

write(*,*) 'setup_react_diff: NX,NY,NZ: ',NX,NY,NZ
dxf = DELTA_X
dx3 = dxf*dxf*dxf
nrow = (NX-2)*(NY-2)*(NZ-1)		! embedded map, Dirichlet conditions on ix=1,NX, iy=1,NY, iz=NZ
maxnz = MAX_CHEMO*nrow
if (allocated(amap)) deallocate(amap)
if (allocated(ja)) deallocate(ja)
if (allocated(ia)) deallocate(ia)
allocate(amap(maxnz,0:3))
allocate(ja(maxnz))
allocate(ia(nrow+1))

dxb3 = dxb*dxb*dxb
nrow_b = NXB*NYB*NZB
maxnz = MAX_CHEMO*nrow_b
if (allocated(amap_b)) deallocate(amap_b)
if (allocated(ja_b)) deallocate(ja_b)
if (allocated(ia_b)) deallocate(ia_b)
allocate(amap_b(maxnz,0:3))
allocate(ja_b(maxnz))
allocate(ia_b(nrow_b+1))

emapfile = ''
write(emapfile,'(a,i2.0,a)') 'emap',NX,'.dat'
write(nflog,*) 'setup_react_diff: ',emapfile
write(*,*) 'make emapfile: ',emapfile
call make_sparse_emap(emapfile,.true.)
write(*,*) 'made emapfile: ',emapfile

bmapfile = ''
write(bmapfile,'(a,i2.0,a)') 'bmap',NXB,'.dat'
write(nflog,*) 'setup_react_diff: ',bmapfile
write(*,*) 'make bmapfile: ',bmapfile
call make_sparse_map(bmapfile,.false.)
write(*,*) 'made bmapfile: ',bmapfile

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (ichemo > TRACER) then
		chemo(ichemo)%bdry_conc = 0
	endif
	Cextra_all(:,:,:,ichemo) = chemo(ichemo)%bdry_conc
	Caverage(:,:,:,ichemo) = chemo(ichemo)%bdry_conc
	write(*,*) 'call update_Cex_Cin: ',ichemo
	call update_Cex_Cin_const(ichemo)	! initialise with SS values
	Fcurr => Cflux(:,:,:,ichemo)
	call getF_const(ichemo,Fcurr)	! estimates all fine grid pt fluxes Cflux(:,:,:,:) from Cextra(:,:,:,:)
! Sum fine grid fluxes to initialise Fcurr_b, Fprev_b
	C0 = chemo(ichemo)%bdry_conc
	Cprev => chemo(ichemo)%Cprev
	Fprev => chemo(ichemo)%Fprev
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
	Cave_b => chemo(ichemo)%Cave_b
	
	Cave_b = C0
	Cprev_b = C0
	Cprev = C0
	Fprev = Cflux(:,:,:,ichemo)
	call makeF_b(Fprev_b,Fprev,DELTA_T,zero)
	Fcurr_b = Fprev_b
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! The flux values on the coarse grid are derived from the values on the fine grid by
! appropriate summation.  It's important that the total flux is consistent.
!-------------------------------------------------------------------------------------------
subroutine makeF_b(F_b,F,dt,zero)
integer, parameter :: dN = NRF/2
real(REAL_KIND) :: F_b(:,:,:), F(:,:,:), dt
logical :: zero
real(REAL_KIND) :: Fsum_b, Fsum, Fmin
real(REAL_KIND) :: wt(-dN:dN,-dN:dN,-dN:dN)
real(REAL_KIND) :: xfac, yfac, zfac
integer :: idx, idy, idz, xb0, yb0, idxb, idyb, zb1
integer :: ixb, iyb, izb, ix0, iy0, iz0, ix, iy, iz
integer :: ichemo

xb0 = (NXB+1)/2			! these must all be integer - check NX,NY,NZ,NXB,NYB,NZB,NRF earlier
idxb = (NX-1)/(2*NRF)
yb0 = (NYB+1)/2
idyb = (NY-1)/(2*NRF)
zb1 = (NZ-1)/NRF + 1

do idx = -dN,dN
do idy = -dN,dN
do idz = -dN,dN
	if (idx > -dN .and. idx < dN) then
		xfac = 1
	else
		xfac = 0.5
	endif
	if (idy > -dN .and. idy < dN) then
		yfac = 1
	else
		yfac = 0.5
	endif
	if (idz > -dN .and. idz < dN) then
		zfac = 1
	else
		zfac = 0.5
	endif
	wt(idx,idy,idz) = xfac*yfac*zfac
enddo
enddo
enddo

Fsum_b = 0
F_b = 0
do ixb = xb0-idxb,xb0+idxb
	ix0 = (ixb - xb0)*NRF + (NX+1)/2
!	write(*,*) 'ix0: ',ix0
	do iyb = yb0-idyb,yb0+idyb
!		write(*,*) 'iyb: ',iyb
		iy0 = (iyb - yb0)*NRF + (NY+1)/2
!		write(*,*) 'iy0: ',iy0
		do izb = 1,zb1
			iz0 = (izb - 1)*NRF + 1
			do idx = -dN,dN
				ix = ix0 + idx
				if (ix < 1 .or. ix > NX) cycle
				do idy = -dN,dN
					iy = iy0 + idy
					if (iy < 1 .or. iy > NY) cycle
					do idz = -dN,dN
						iz = iz0 + idz
						if (iz < 1 .or. iz > NZ) cycle
						F_b(ixb,iyb,izb) = F_b(ixb,iyb,izb) + wt(idx,idy,idz)*F(ix,iy,iz)
						Fsum_b = Fsum_b + wt(idx,idy,idz)*F(ix,iy,iz)
					enddo
				enddo
			enddo
		enddo
	enddo
enddo
zero = (Fsum_b == 0)	
end subroutine

!-------------------------------------------------------------------------------------------
! This version is for the embedded grid
!-------------------------------------------------------------------------------------------
subroutine make_csr_SS(a, ichemo, Cave, Fcurr, rhs)
integer :: ichemo
real(REAL_KIND) :: a(:), Cave(:,:,:), Fcurr(:,:,:), rhs(:)
integer :: k, ix, iy, iz, krow, kcol, nc
integer :: nc_max = 10	! just a wild guess but not a bad one
real(REAL_KIND) :: Kdiff, Kr, Vex, Cbdry
logical, save :: first = .true.

krow = 0
do k = 1,nnz
	if (k == ia(krow+1)) krow = krow+1
	kcol = ja(k)
	a(k) = amap(k,0)
!	if (amap(k,2) == NY .and. kcol == krow-NY) then
!		if (a(k) /= -2) then
!			write(*,*) 'Error in OXYGEN bdry adjustment'
!			stop
!		endif
!		a(k) = -1
!	endif
enddo

do ix = 2,NX-1
	do iy = 2,NY-1
		do iz = 1,NZ-1
			krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
			! This is very crude!!!!!!!!!!!!!!! e.g. not centred on the grid pt
			nc = grid(ix,iy,iz)%nc
			Kdiff = chemo(ichemo)%medium_diff_coef
			Kdiff = Kdiff*(1 - chemo(ichemo)%diff_reduction_factor*min(nc,nc_max)/nc_max)
			Kr = 1/(dxf*Kdiff)
			rhs(krow) = -Kr*Fcurr(ix,iy,iz)
		enddo
	enddo
enddo

ix = 2
do iy = 2,NY-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(1,iy,iz)
		rhs(krow) = rhs(krow) + Cbdry		
	enddo
enddo
ix = NX-1
do iy = 2,NY-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(NX,iy,iz)
		rhs(krow) = rhs(krow) + Cbdry		
	enddo
enddo
iy = 2
do ix = 2,NX-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(ix,1,iz)
		rhs(krow) = rhs(krow) + Cbdry		
	enddo
enddo
iy = NY-1
do ix = 2,NX-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(ix,NY,iz)
		rhs(krow) = rhs(krow) + Cbdry		
	enddo
enddo
iz = NZ-1
do ix = 2,NX-1
	do iy = 2,NY-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(ix,iy,NZ)
		rhs(krow) = rhs(krow) + Cbdry		
	enddo
enddo
first = .false.
end subroutine

!-------------------------------------------------------------------------------------------
! This is the version for the coarse grid.
! Use variable numbering (ix,iy) -> k = (ix-1)*NY + iy
! Since the equation is now dM/dt = ... need to derive different expression for Kr
! Need to add treatment of top boundary for O2
! Need to add decay!  For now add it explicitly in solver
!-------------------------------------------------------------------------------------------
subroutine make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs, zero)
integer :: ichemo
real(REAL_KIND) :: dt, a_b(:), Cave_b(:,:,:), Cprev_b(:,:,:), Fcurr_b(:,:,:), Fprev_b(:,:,:), rhs(:)
logical :: zero
integer :: ixb, iyb, izb, k, i, krow, kcol, nc
real(REAL_KIND) :: Kdiff, Kr, Cbdry, Fsum
integer, parameter :: m = 3

zero = .true.
Kdiff = chemo(ichemo)%medium_diff_coef
Fsum = 0
krow = 0
do k = 1,nnz_b
	if (k == ia_b(krow+1)) krow = krow+1
	kcol = ja_b(k)
	if (amap_b(k,0) == 2*m) then
		Kr = dxb*dxb/Kdiff
		a_b(k) = 3*Kr/(2*dt) + 2*m	! ... note that Kdiff should depend on (ixb,iyb,izb), ultimately on # of cells
	else
		a_b(k) = amap_b(k,0)
	endif
	if (ichemo == OXYGEN) then
		if (amap_b(k,3) == NZB .and. kcol == krow-1) then		! check this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if (a_b(k) /= -2) then
				write(*,*) 'Error in OXYGEN bdry adjustment'
				stop
			endif
			a_b(k) = -1
		endif
	endif
enddo
!!$omp parallel do private(iyb, ixb, krow, Kr)
do izb = 1,NZB
	do iyb = 1,NYB
		do ixb = 1,NXB
			krow = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
			Kr = dxb*dxb/Kdiff
			rhs(krow) = Kr*((-2*Fcurr_b(ixb,iyb,izb) + Fprev_b(ixb,iyb,izb))/dxb3 + (1./(2*dt))*(4*Cave_b(ixb,iyb,izb) - Cprev_b(ixb,iyb,izb)))
			if (rhs(krow) /= 0) zero = .false.
		enddo
	enddo
enddo
!!$omp end parallel do
if (ichemo == OXYGEN) then
	Cbdry = chemo(ichemo)%bdry_conc
	izb = NZB
	do ixb = 1,NXB
		do iyb = 1,NYB
			krow = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
			rhs(krow) = rhs(krow) + Cbdry		
		enddo
	enddo
endif
end subroutine

!-------------------------------------------------------------------------------------------
! Update Cextra(:,:,:) and thereby cp%Cex(ichemo) and cp%Cin(ichemo) from Caverage(:,:,:,ichemo), 
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin_const(ichemo)
integer :: ichemo
integer :: kcell
real(REAL_KIND) :: Kin, Kout, dC, dCdt, dMdt
type(cell_type), pointer :: cp
real(REAL_KIND), pointer :: Cextra(:,:,:)
!write(*,*) 'update_Cex_Cin_const'

Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!

Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	call extra_concs_const(kcell, Cextra, cp%Cex(ichemo))
	cp%Cin(ichemo) = getCin_SS(kcell,ichemo,cp%V,cp%Cex(ichemo))
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! Update Cextra(:,:,:) and thereby cp%Cex(ichemo) and cp%Cin(ichemo) from Caverage(:,:,:,ichemo), 
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin_dCdt_const(ichemo, dt)
integer :: ichemo
real(REAL_KIND) :: dt
integer :: kcell, ix, iy, iz
real(REAL_KIND) :: alfa(3), Clast, Kin, Kout, dC, dCexdt, dMdt
type(cell_type), pointer :: cp
real(REAL_KIND), pointer :: Cextra(:,:,:)
real(REAL_KIND), pointer :: Cprev(:,:,:)

!write(*,*) 'update_Cex_Cin_dCdt_const: ',ichemo

Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!
Cprev => chemo(ichemo)%Cprev
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
!$omp parallel do private(cp, ix, iy, iz, alfa, Clast, dCexdt)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	call grid_interp(kcell, alfa)
	ix = cp%site(1)
	iy = cp%site(2)
	iz = cp%site(3)
	cp%Cex(ichemo) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1)
	Clast = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cprev(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cprev(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cprev(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cprev(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cprev(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cprev(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cprev(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cprev(ix+1,iy,iz+1)
    dCexdt = (cp%Cex(ichemo) - Clast)/dt
	cp%Cin(ichemo) = getCin(kcell,ichemo,cp%V,cp%Cex(ichemo),dCexdt)
enddo
!$omp end parallel do
end subroutine

!-------------------------------------------------------------------------------------------
! Estimate total flux values associated with each fine grid pt, consistent with the current
! average concentrations Cave in the neighbourhood of the grid pt.  The grid pt flux is
! determined from the uptake fluxes of nearby cells, and a positive computed flux F 
! corresponds to a reduction in grid pt concentration.
! In an iterative procedure, the average grid pt concentrations C are used to estimate
! Cex and Cin for each nearby cell, which then provides an estimate of the effective
! extracellular concentration associated with each grid pt.
! The membrane flux is the equilibrium value, i.e. 
! rate of consumption is balanced by the rate of mass transport across the membrane.  
! Note: It might make sense to compute all the cell concentrations at the same time,
! since the interpolations are identical for each constituent, at a given cell.
! The current best estimate of the extracellular concentrations at the fine grid pts
! is stored in Cextra(:,:,:,ichemo)
! This is initialised to Caverage(:,:,:,:)
! Not necessarily SS
!-------------------------------------------------------------------------------------------
subroutine getF_const(ichemo, Cflux_const)
integer :: ichemo
real(REAL_KIND) :: Cflux_const(:,:,:)
real(REAL_KIND) :: Kin, Kout
integer :: kcell
type(cell_type), pointer :: cp

!write(*,*) 'getF_const: ',ichemo,nlist
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out

! Compute cell fluxes cp%dMdt
!$omp parallel do private(cp)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cp%dMdt(ichemo) = Kin*cp%Cex(ichemo) - Kout*cp%Cin(ichemo)
enddo
!$omp end parallel do

! Estimate grid pt flux values F
call make_grid_flux(ichemo,Cflux_const)
end subroutine

!-------------------------------------------------------------------------------------------
! alfa is the amount of the previous flux
!-------------------------------------------------------------------------------------------
subroutine make_grid_flux(ichemo,Cflux_const)
integer :: ichemo
real(REAL_KIND) :: Cflux_const(:,:,:)
integer :: kcell, k, cnr(3,8)
type(cell_type), pointer :: cp
real(REAL_KIND) :: alfa = 0.7

Cflux_const = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cnr = cp%cnr
	do k = 1,8
		Cflux_const(cnr(1,k),cnr(2,k),cnr(3,k)) = Cflux_const(cnr(1,k),cnr(2,k),cnr(3,k)) + cp%dMdt(ichemo)*cp%wt(k)
	enddo
enddo
Cflux_const(:,:,:) = (1-alfa)*Cflux_const(:,:,:) + alfa*Cflux_prev(:,:,:,ichemo)
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine make_grid_flux_weights
integer :: ic, ix, iy, iz, kcell
type(cell_type), pointer :: cp

do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ix = cp%site(1)
	iy = cp%site(2)
	iz = cp%site(3)
	cp%cnr(:,1) = [ix, iy, iz]
	cp%cnr(:,2) = [ix, iy+1, iz]
	cp%cnr(:,3) = [ix, iy, iz+1]
	cp%cnr(:,4) = [ix, iy+1, iz+1]
	cp%cnr(:,5) = [ix+1, iy, iz]
	cp%cnr(:,6) = [ix+1, iy+1, iz]
	cp%cnr(:,7) = [ix+1, iy, iz+1]
	cp%cnr(:,8) = [ix+1, iy+1, iz+1]
	call grid_flux_weights(kcell, cp%cnr, cp%wt)
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! This valid for steady-state with any constituent, with decay, when Cex is known
! In fact it is currently valid for oxygen and glucose.
! This is effectively the same as getCin() with dCexdt = 0
!-------------------------------------------------------------------------------------------
function getCin_SS(kcell,ichemo, Vin, Cex) result(Cin)
integer :: kcell, ichemo
real(REAL_KIND) :: Vin, Cex, Cin, Cex_t, Cin_t
real(REAL_KIND) :: Kin, Kout, Kd, Kmax, VKdecay, C0, a, b, c, D, r(3)
integer :: i, n, ictyp, idrug, im
real(REAL_KIND) :: CO2, C_parent, C_metab1
real(REAL_KIND) :: C2_0, C2_1, C2_2, KO2_0, KO2_1, KO2_2, Kmet0_0, Kmet0_1, Kmet0_2
real(REAL_KIND) :: K1, K2, Km, Vmax
type(drug_type), pointer :: dp

!write(*,*) 'getCin_SS: istep,ichemo,Vin,Cex: ',istep,ichemo,Vin,Cex
!if (Cex < 0) then
!	write(logmsg,*) 'getCin_SS: istep,ichemo,Vin,Cex: ',istep,ichemo,Vin,Cex
!	call logger(logmsg)
!	stop
!endif
if (Cex <= 0) then
	Cex = 0
	Cin = 0
	return
endif
ictyp = cell_list(kcell)%celltype
CO2 = cell_list(kcell)%Cin(OXYGEN)
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
Kd = chemo(ichemo)%decay_rate
if (ichemo <= GLUCOSE) then
	Kmax = chemo(ichemo)%max_cell_rate
	VKdecay = Vin*Kd
	C0 = chemo(ichemo)%MM_C0
	if (chemo(ichemo)%Hill_N == 2) then
		b = C0*C0
		a = (Kmax + b*VKdecay - Kin*Cex)/(Kout+VKdecay)
		c = -b*Cex*Kin/(Kout + VKdecay)
		call cubic_roots(a,b,c,r,n)
		if (n == 1) then
			Cin = r(1)
		else
			n = 0
			do i = 1,3
				if (r(i) > 0) then
					n = n+1
					Cin = r(i)
				endif
			enddo
			if (n > 1) then
				write(*,*) 'getCin_SS: two roots > 0: ',r
				stop
			endif
		endif
	elseif (chemo(ichemo)%Hill_N == 1) then
		b = (Kmax + Kout*C0 + 2*VKdecay*C0 - Kin*Cex)/(Kout + VKdecay)
		c = -C0*Cex*Kin/(Kout + VKdecay)
		D = sqrt(b*b - 4*c)
		Cin = (D - b)/2
	endif
elseif (ichemo == TRACER) then
	
else	! parent drug or drug metabolite
	idrug = (ichemo - TRACER - 1)/3 + 1
	dp => drug(idrug)
	im = ichemo - TRACER - 1 - 3*(idrug-1)
	if (im == 0) then		! parent
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Km = dp%Km(ictyp,0)
		Vmax = dp%Vmax(ictyp,0)
		K1 = (1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0
		K2 = K1*Vmax/Kmet0_0
		if (K2 /= 0) then	!quadratic: a.x^2 + b.x + c = 0
			a = K1 + Kd + Kout/Vin
			b = a*Km + K2 - Kin*Cex/Vin
			c = -Kin*Cex*Km/Vin
			b = b/a
			c = c/a
			D = sqrt(b*b - 4*c)
			Cin = (D - b)/2
		else				! linear: a.x + b = 0
			a = K1 + Kd + Kout/Vin
			b = -Kin*Cex/Vin
			Cin = -b/a
		endif
		if (Cin > Cex) then
			write(*,'(a,2e12.3)') 'Cex,Cin: ',Cex,Cin
			stop
		endif
	elseif (im == 1) then	! metab1
		CO2 = cell_list(kcell)%Cin(OXYGEN)
		C_parent = cell_list(kcell)%Cin(ichemo-1)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = drug(idrug)%Kmet0(ictyp,1)
		Cin = ((1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0*C_parent + Kin*Cex/Vin) &
		     /((1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1 + Kd + Kout/Vin)
	elseif (im == 2) then	! metab2
		CO2 = cell_list(kcell)%Cin(OXYGEN)
		C_metab1 = cell_list(kcell)%Cin(ichemo-1)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = dp%Kmet0(ictyp,1)
		C2_2 = dp%C2(ictyp,2)
		KO2_2 = dp%KO2(ictyp,2)
		Kmet0_2 = dp%Kmet0(ictyp,2)
		Cin = ((1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1*C_metab1 + Kin*Cex/Vin) &
		     /((1 - C2_2 + C2_2*KO2_2/(KO2_2 + CO2))*Kmet0_2 + Kd + Kout/Vin)
	endif
endif
end function

!-------------------------------------------------------------------------------------------
! Use an approximation to dCin/dt derived from dCex/dt to arrive at a better estimate of Cin and flux
!-------------------------------------------------------------------------------------------
function getCin(kcell, ichemo, Vin, Cex, dCexdt) result(Cin)
integer :: kcell, ichemo
real(REAL_KIND) :: Vin, Cex, dCexdt, Cin
real(REAL_KIND) :: Kin, Kout, Kd, Kmax, VKdecay, dCdt, delta, C0, a, b, c, D, r(3)
integer :: i, n, ictyp, idrug, im
real(REAL_KIND) :: CO2, C_parent, C_metab1
real(REAL_KIND) :: C2_0, C2_1, C2_2, KO2_0, KO2_1, KO2_2, Kmet0_0, Kmet0_1, Kmet0_2
real(REAL_KIND) :: Km, Vmax, K1_0, K2_0, K1_1, K2_1, K1_2, K2_2
type(drug_type), pointer :: dp
type(cell_type), pointer :: cp


if (Cex < 0) then
	Cex = 0
	Cin = 0
	return
endif
cp => cell_list(kcell)
ictyp = cp%celltype
CO2 = cp%Cin(OXYGEN)
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
Kd = chemo(ichemo)%decay_rate
if (ichemo <= GLUCOSE) then
	delta = Vin*dCdt
	Kmax = chemo(ichemo)%max_cell_rate
	VKdecay = Vin*Kd
	C0 = chemo(ichemo)%MM_C0
	if (chemo(ichemo)%Hill_N == 2) then
		b = C0*C0
		a = (Kmax + b*VKdecay - (Kin*Cex-delta))/(Kout+VKdecay)
		c = -b*(Kin*Cex-delta)/(Kout + VKdecay)
		call cubic_roots(a,b,c,r,n)
		if (n == 1) then
			Cin = r(1)
		else
			n = 0
			do i = 1,3
				if (r(i) > 0) then
					n = n+1
					Cin = r(i)
				endif
			enddo
			if (n > 1) then
				write(*,*) 'getCin_SS: two roots > 0: ',r
				stop
			endif
		endif
	elseif (chemo(ichemo)%Hill_N == 1) then
		b = (Kmax + Kout*C0 + 2*VKdecay*C0 - (Kin*Cex-delta))/(Kout + VKdecay)
		c = -C0*(Kin*Cex-delta)/(Kout + VKdecay)
		D = sqrt(b*b - 4*c)
		Cin = (D - b)/2
	endif
elseif (ichemo == TRACER) then
	stop
else	! parent drug or drug metabolite
	idrug = (ichemo - TRACER - 1)/3 + 1
	dp => drug(idrug)
	im = ichemo - TRACER - 1 - 3*(idrug-1)
	if (im == 0) then		! parent
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Km = dp%Km(ictyp,0)
		Vmax = dp%Vmax(ictyp,0)
		K1_0 = (1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0
		K2_0 = K1_0*Vmax/Kmet0_0
		dCdt = (Kin*dCexdt/Vin)/(Kout/Vin + Kd + K1_0 + K2_0*Km/(Km + CO2)**2)
		if (K2_0 /= 0) then	!quadratic: a.x^2 + b.x + c = 0
			a = K1_0 + Kd + Kout/Vin
			b = a*Km + K2_0 - (Kin*Cex/Vin - dCdt)
			c = -Km*(Kin*Cex/Vin - dCdt)
			b = b/a
			c = c/a
			D = sqrt(b*b - 4*c)
			Cin = (D - b)/2
		else				! linear: a.x + b = 0
			a = K1_0 + Kd + Kout/Vin
			b = -(Kin*Cex/Vin - dCdt)
			Cin = -b/a
		endif
		Cin = max(Cin,0.0)
	elseif (im == 1) then	! metab1
		C_parent = cp%Cin(ichemo-1)
		C2_0 = dp%C2(ictyp,0)
		KO2_0 = dp%KO2(ictyp,0)
		Kmet0_0 = dp%Kmet0(ictyp,0)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = drug(idrug)%Kmet0(ictyp,1)
		K1_0 = (1 - C2_0 + C2_0*KO2_0/(KO2_0 + CO2))*Kmet0_0
		K1_1 = (1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1
		dCdt = (Kin*dCexdt/Vin + K1_0*cp%dCdt(ichemo-1))/(Kout/Vin + Kd + K1_1)
		Cin = (K1_0*C_parent + (Kin*Cex/Vin - dCdt))/(K1_1 + Kd + Kout/Vin)
		Cin = max(Cin,0.0)
	elseif (im == 2) then	! metab2
		C_metab1 = cp%Cin(ichemo-1)
		C2_1 = dp%C2(ictyp,1)
		KO2_1 = dp%KO2(ictyp,1)
		Kmet0_1 = dp%Kmet0(ictyp,1)
		C2_2 = dp%C2(ictyp,2)
		KO2_2 = dp%KO2(ictyp,2)
		Kmet0_2 = dp%Kmet0(ictyp,2)
		K1_1 = (1 - C2_1 + C2_1*KO2_1/(KO2_1 + CO2))*Kmet0_1
		K1_2 = (1 - C2_2 + C2_2*KO2_2/(KO2_2 + CO2))*Kmet0_2
		dCdt = (Kin*dCexdt/Vin + K1_1*cp%dCdt(ichemo-1))/(Kout/Vin + Kd + K1_2)
		Cin = (K1_1*C_metab1 + (Kin*Cex/Vin - dCdt))/(K1_2 + Kd + Kout/Vin)
		Cin = max(Cin,0.0)
	endif
endif
end function

!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
subroutine cubic_roots(a, b, c, r, n)
real(REAL_KIND) :: a, b, c, r(3)
integer :: n
real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB

QQ = (a*a - 3*b)/9
RR = (2*a*a*a - 9*a*b + 27*c)/54
Q3 = QQ*QQ*QQ
R2 = RR*RR
if (R2 < Q3) then
	n = 3
	theta = acos(RR/sqrt(Q3))
	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
else
	n = 1
	AA = -sign(1.d0,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
	if (AA == 0) then
		BB = 0
	else
		BB = QQ/AA
	endif
	r(1) = AA + BB - a/3
endif
end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine unpack_csr(a)
real(REAL_KIND) :: a(:)
integer :: k, irow, icol
real(REAL_KIND) :: asum

write(nflog,*) 'nrow: ',nrow
allocate(afull(nrow,nrow))
irow = 0
do k = 1,nnz
	if (k == ia(irow+1)) irow = irow+1
	icol = ja(k)
	afull(irow,icol) = a(k)
enddo
!do irow = 1,nrow
!	asum = 0
!	do icol = 1,nrow
!		asum = asum + afull(irow,icol)*9
!	enddo
!!	write(*,'(i6,2f8.3)') irow,asum,rhs(irow)
!enddo
!write(*,*) 'rhs:'
!write(*,'(25f7.2)') rhs(:)
end subroutine

!--------------------------------------------------------------------------------------
! Cex values on the boundary of the fine grid are interpolated from Cave_b on the coarse grid
! The coarse grid indices that coincide with the boundary of the fine grid are:
!	ixb = xb0-idxb, xb0+idxb = xb1, xb2
!	iyb = yb0-idyb, yb0+idyb = yb1, yb2
!	izb = zb1
! (see below)
! There are 5 faces.
!
! Interpolate over the whole fine grid?
!--------------------------------------------------------------------------------------
subroutine interpolate_Cave(Cave, Cave_b)
real(REAL_KIND) :: Cave(:,:,:), Cave_b(:,:,:)
integer :: idx, idy, idz, xb0, yb0, idxb, idyb, xb1, xb2, yb1, yb2, zb1, zb2
integer :: ixb, iyb, izb, ix0, iy0, iz0, ix, iy, iz, nsum
real(REAL_KIND) :: ax, ay, az, asum(2), Cnew
real(REAL_KIND) :: alfa = 0.3

xb0 = (NXB+1)/2			
idxb = (NX-1)/(2*NRF)
xb1 = xb0 - idxb
xb2 = xb0 + idxb
yb0 = (NYB+1)/2
idyb = (NY-1)/(2*NRF)
yb1 = yb0 - idyb
yb2 = yb0 + idyb
zb1 = 1
zb2 = (NZ-1)/NRF + 1

! Left and right faces: ixb = xb1, xb2  
!                       iyb = yb1..yb2
!						izb = zb1..zb2
asum = 0
nsum = 0
do iyb = yb1,yb2-1
	iy0 = (iyb - yb0)*NRF + (NY+1)/2
	do izb = zb1,zb2-1
		iz0 = (izb - 1)*NRF + 1
		do idy = 0,4
			iy = iy0 + idy
			ay = (4. - idy)/4.
			do idz = 0,4
				iz = iz0 + idz
				az = (4. - idz)/4.
				ixb = xb1
				ix = (ixb - xb0)*NRF + (NX+1)/2
				nsum = nsum + 1
				Cave(ix,iy,iz) =     ay*az*Cave_b(ixb, iyb,   izb) + &
								(1-ay)*az*Cave_b(ixb, iyb+1, izb) + &
								ay*(1-az)*Cave_b(ixb, iyb,   izb+1) + &
							(1-ay)*(1-az)*Cave_b(ixb, iyb+1, izb+1)
				asum(1) = asum(1) + Cave(ix,iy,iz)
				ixb = xb2
				ix = (ixb - xb0)*NRF + (NX+1)/2
				Cnew =     ay*az*Cave_b(ixb,iyb,  izb) + &
								(1-ay)*az*Cave_b(ixb,iyb+1,izb) + &
								ay*(1-az)*Cave_b(ixb,iyb,  izb+1) + &
							(1-ay)*(1-az)*Cave_b(ixb,iyb+1,izb+1)
				Cave(ix,iy,iz) = alfa*Cnew + (1-alfa)*Cave(ix,iy,iz)
				asum(2) = asum(2) + Cave(ix,iy,iz)
			enddo
		enddo
	enddo
enddo
!write(*,*) 'Average left and right faces: ',asum/nsum

! Back and front faces: iyb = yb1, yb2
!                       ixb = xb1..xb2
!						izb = zb1..zb2
asum = 0
nsum = 0
do ixb = xb1,xb2-1
	ix0 = (ixb - xb0)*NRF + (NX+1)/2
	do izb = zb1,zb2-1
		iz0 = (izb - 1)*NRF + 1
		do idx = 0,4
			ix = ix0 + idx
			ax = (4. - idx)/4.
			do idz = 0,4
				iz = iz0 + idz
				az = (4. - idz)/4.
				iyb = yb1
				iy = (iyb - yb0)*NRF + (NY+1)/2
				nsum = nsum + 1
				Cave(ix,iy,iz) =     ax*az*Cave_b(ixb,   iyb, izb) + &
								(1-ax)*az*Cave_b(ixb+1, iyb, izb) + &
								ax*(1-az)*Cave_b(ixb,   iyb, izb+1) + &
							(1-ax)*(1-az)*Cave_b(ixb+1, iyb, izb+1)
				asum(1) = asum(1) + Cave(ix,iy,iz)
				iyb = yb2
				iy = (iyb - yb0)*NRF + (NY+1)/2
				Cnew =     ax*az*Cave_b(ixb,   iyb, izb) + &
								(1-ax)*az*Cave_b(ixb+1, iyb, izb) + &
								ax*(1-az)*Cave_b(ixb,   iyb, izb+1) + &
							(1-ax)*(1-az)*Cave_b(ixb+1, iyb, izb+1)
				Cave(ix,iy,iz) = alfa*Cnew + (1-alfa)*Cave(ix,iy,iz)
				asum(2) = asum(2) + Cave(ix,iy,iz)
			enddo
		enddo
	enddo
enddo
!write(*,*) 'Average back and front faces: ',asum/nsum

! Top face: izb = zb2
!           ixb = xb1..xb2
!			iyb = yb1..yb2
nsum = 0
asum = 0
do ixb = xb1,xb2-1
	ix0 = (ixb - xb0)*NRF + (NX+1)/2
!	write(*,*) 'ix0: ',ix0
	do iyb = yb1,yb2-1
		iy0 = (iyb - yb0)*NRF + (NY+1)/2
!		write(*,*) 'iy0: ',iy0
		do idx = 0,4
			ix = ix0 + idx
			ax = (4. - idx)/4.
			do idy = 0,4
				iy = iy0 + idy
				ay = (4. - idy)/4.
				izb = zb2
				iz = (izb - 1)*NRF + 1
				nsum = nsum + 1
				Cnew =     ax*ay*Cave_b(ixb,   iyb,   izb) + &
								(1-ax)*ay*Cave_b(ixb+1, iyb,   izb) + &
								ax*(1-ay)*Cave_b(ixb,   iyb+1, izb) + &
							(1-ax)*(1-ay)*Cave_b(ixb+1, iyb+1, izb)
				Cave(ix,iy,iz) = alfa*Cnew + (1-alfa)*Cave(ix,iy,iz)
			enddo
		enddo
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------------
! See D:\ITSOL\tests\input
!
! Time for 20 iterations with NX=1000
! ILUtype=1 33 sec
! ILUtype=2 is hoplessly slow with the parameters specified
! ILUtype=3 56 sec
! ILUtype=4 197 sec
!
! Note: Solving for each constituent in parallel with OpenMP like this is possible
!       only if the constituent reactions are independent
! For drugs, the three related constituents are solved for sequentially in the fine grid.
! First drug, then metab1 (using drug results), then metab2 (using metab1 results).
!-------------------------------------------------------------------------------------- 
subroutine diff_solver(dt)
real(REAL_KIND) :: dt
integer :: i, k, k1, ix, iy, iz, irow, icol, kc, ic, icc, it
integer :: ixb, iyb, izb
integer :: ichemo, ierr, nfill, iters, maxits, im_krylov
real(REAL_KIND) :: R, tol, tol_b, asum, t, Vex_curr, Vex_next, Vin_curr, Vin_next, fdecay, Csum, dCsum, msum
real(REAL_KIND), allocatable :: x(:), rhs(:)
real(REAL_KIND), pointer :: Cave(:,:,:), Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
real(REAL_KIND), allocatable :: a(:), a_b(:)
real(REAL_KIND) :: alfa = 0.7
real(REAL_KIND) :: dCtol = 1.0e-4
integer :: ILUtype = 1
integer :: nfinemap, finemap(MAX_CHEMO)
integer :: im, im1, im2, ichemof
logical :: zeroF(MAX_CHEMO), zeroC(MAX_CHEMO)
logical :: done
logical :: do_fine = .true.
logical :: use_const = .true.

nfinemap = 0
do ic = 1,nchemo
	ichemo = chemomap(ic)
	if (ichemo <= TRACER) then
		nfinemap = nfinemap + 1
		finemap(nfinemap) = ichemo
	elseif (mod(ichemo-TRACER-1,3) == 0) then
		nfinemap = nfinemap + 1
		finemap(nfinemap) = ichemo		! idrug = (ichemo-TRACER-1)/3 + 1
	endif
enddo		

nfill = 1	! Level of fill for ILUK preconditioner
tol = 1.0d-6
tol_b = 1.0d-6
im_krylov = 60	! dimension of Krylov subspace in (outer) FGMRES
maxits = 100

! Compute all steady-state grid point fluxes in advance from Cextra(:,:,:,:): Cflux(:,:,:,:)

!$omp parallel do private(Cave, Fcurr, Cave_b, Cprev_b, Fprev_b, Fcurr_b, a_b, x, rhs, ix, iy, iz, ixb, iyb, izb, it, done, ichemo, icc, k, dCsum, msum, iters, ierr)
do ic = 1,nchemo
	ichemo = chemomap(ic)
	if (chemo(ichemo)%constant) cycle
	write(*,'(a,i2)') 'coarse grid: ichemo: ',ichemo
	write(nflog,'(a,i2)') 'coarse grid: ichemo: ',ichemo
	ichemo_curr = ichemo
	icc = ichemo - 1
	allocate(rhs(nrow_b))
	allocate(x(nrow_b))
	allocate(a_b(MAX_CHEMO*nrow_b))
	Fcurr => Cflux(:,:,:,ichemo)
	Cave_b => chemo(ichemo)%Cave_b
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
	Cave => Caverage(:,:,:,ichemo)
!	write(*,*) 'Cave_b:'
!	write(*,'(5e15.6)') Cave_b(NXB/2,NYB/2,:)
		
	Fprev_b = Fcurr_b
	call makeF_b(Fcurr_b, Fcurr, dt,zeroF(ichemo))
	call make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs, zeroC(ichemo))		! coarse grid

	! Solve Cave_b(t+dt) on coarse grid
	!----------------------------------
	call itsol_create_matrix(icc,nrow_b,nnz_b,a_b,ja_b,ia_b,ierr)
	!write(nflog,*) 'itsol_create_matrix: ierr: ',ierr
		
	if (ILUtype == 1) then
		call itsol_create_precond_ILUK(icc,nfill,ierr)
	!	write(nflog,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
	elseif (ILUtype == 2) then
		call itsol_create_precond_VBILUK(icc,nfill,ierr)
	!	write(*,*) 'itsol_create_precond_VBILUK: ierr: ',ierr 
	elseif (ILUtype == 3) then
		call itsol_create_precond_ILUT(icc,nfill,tol_b,ierr)
	!	write(*,*) 'itsol_create_precond_ILUT: ierr: ',ierr 
	elseif (ILUtype == 4) then
		call itsol_create_precond_ARMS(icc,nfill,tol_b,ierr)
	!	write(*,*) 'itsol_create_precond_ARMS: ierr: ',ierr 
	endif

	do izb = 1,NZB
		do iyb = 1,NYB
			do ixb = 1,NXB
				k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
				x(k) = Cave_b(ixb,iyb,izb)		! initial guess
			enddo
		enddo
	enddo
	if (.not.zeroC(ichemo)) then
	!	write(nflog,*) 'call itsol_solve_fgmr_ILU'
		call itsol_solve_fgmr_ILU(icc, rhs, x, im_krylov, maxits, tol_b, iters, ierr)
	!	write(nflog,*) 'itsol_solve_fgmr_ILU: Cave_b: ierr, iters: ',ierr,iters
	else
		write(*,*) 'no solve, zeroC: ',ichemo
	endif
	call itsol_free_precond_ILU(icc, ierr)
!	write(nflog,*) 'did itsol_free_precond_ILU'
	call itsol_free_matrix(icc, ierr)
!	write(nflog,*) 'did itsol_free_matrix'

	Cprev_b = Cave_b
	fdecay = 1 - chemo(ichemo)%decay_rate*dt
	msum = 0
	do izb = 1,NZB
		do iyb = 1,NYB
			do ixb = 1,NXB
				k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
				Cave_b(ixb,iyb,izb) = fdecay*x(k)
				msum = msum + x(k)*dxb3		! this sums the mass of constituent in mumols
				if (x(k) < 0) then
					write(nflog,*) 'Cave_b < 0: ',ixb,iyb,izb,x(k)
				endif 
			enddo
		enddo
	enddo
	! interpolate Cave_b on fine grid boundary
	call interpolate_Cave(Cave, Cave_b)
	deallocate(a_b, x, rhs)
enddo
!$omp end parallel do
	
!$omp parallel do private(Cave, Cprev, Fprev, Fcurr, a, x, rhs, ix, iy, iz, ixb, iyb, izb, it, done, ichemo, icc, k, Csum, dCsum, msum, iters, ichemof, im, im1, im2, ierr)
do ic = 1,nfinemap
	ichemof = finemap(ic)
	if (chemo(ichemof)%constant) cycle
	allocate(rhs(nrow))
	allocate(x(nrow))
	allocate(a(MAX_CHEMO*nrow))
	im1 = 0
	if (ichemof <= TRACER) then
		im2 = 0
	else
		im2 = 2
	endif
	do im = im1, im2
		ichemo = ichemof + im
		if (.not.chemo(ichemo)%present) cycle
		write(*,'(a,i2)') 'fine grid: ichemo: ',ichemo
		ichemo_curr = ichemo
		icc = ichemo - 1
!		write(*,'(a,i2)') 'ichemo: ',ichemo
		Cave => Caverage(:,:,:,ichemo)
		Cprev => chemo(ichemo)%Cprev
		Fprev => chemo(ichemo)%Fprev
		Fcurr => Cflux(:,:,:,ichemo)
		
		it = 0
		do
			it = it + 1
!			write(*,*) 'fine grid: ichemo,it: ',ichemo,it
			if (it == 50) stop
			call make_csr_SS(a, ichemo, Cave, Fcurr, rhs)	! fine grid - note: using the same flux values as the Cave_b solution!
			
			! Solve Cave(t+dt) steady-state on fine grid
			!-------------------------------------------
			call itsol_create_matrix(icc,nrow,nnz,a,ja,ia,ierr)
!			write(nflog,*) 'itsol_create_matrix: ierr: ',ierr
			if (ILUtype == 1) then
				call itsol_create_precond_ILUK(icc,nfill,ierr)
			!	write(*,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
			elseif (ILUtype == 2) then
				call itsol_create_precond_VBILUK(icc,nfill,ierr)
			!	write(*,*) 'itsol_create_precond_VBILUK: ierr: ',ierr 
			elseif (ILUtype == 3) then
				call itsol_create_precond_ILUT(icc,nfill,tol,ierr)
			!	write(*,*) 'itsol_create_precond_ILUT: ierr: ',ierr 
			elseif (ILUtype == 4) then
				call itsol_create_precond_ARMS(icc,nfill,tol,ierr)
			!	write(*,*) 'itsol_create_precond_ARMS: ierr: ',ierr 
			endif

			if (it == 1) then
				do iz = 1,NZ-1
					do iy = 2,NY-1
						do ix = 2,NX-1
							k = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
							x(k) = Cave(ix,iy,iz)		! initial guess
						enddo
					enddo
				enddo
			endif
			
			if (.not.zeroC(ichemo)) then
			!	write(nflog,*) 'call itsol_solve_fgmr_ILU'
				call itsol_solve_fgmr_ILU(icc,rhs, x, im_krylov, maxits, tol, iters, ierr)
			!	write(nflog,*) 'itsol_solve_fgmr_ILU: Cave: ierr, iters: ',ierr,iters
			else
				write(*,*) 'no solve, zeroC: ',ichemo
			endif
			call itsol_free_precond_ILU(icc,ierr)
			call itsol_free_matrix(icc,ierr)
			
			done = .false.
!			dCsum = 0
			Csum = 0
			do iz = 1,NZ-1
				do iy = 2,NY-1
					do ix = 2,NX-1
						k = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
						Cave(ix,iy,iz) = (1-alfa)*x(k) + alfa*Cave(ix,iy,iz)
						Csum = Csum + Cave(ix,iy,iz)
!						dCsum = dCsum + abs((Cave(ix,iy,iz) - Cprev(ix,iy,iz))/Cave(ix,iy,iz))
					enddo
				enddo
			enddo
!			dCsum = dCsum/(NX*NY*NZ)
!			done = (abs(dCsum) < dCtol)
			write(*,'(a,i2,3e12.3)') 'Csum,Cave,Cflux(17,17,17): ',ichemo,Csum,Cave(17,17,17),Cflux(17,17,17,ichemo)
			if (Csum > 0 .and. Csum < 1.0e-10) then
				write(*,*) 'Stop solving for: ',ichemo
				chemo(ichemo)%present = .false.
				! Probably should set all concentrations to 0
			endif
			Cflux_prev(:,:,:,ichemo) = Fcurr
			call update_Cex_Cin_dCdt_const(ichemo,dt)
			call getF_const(ichemo, Fcurr)
			Cprev = Cave
				
			done = .true.	! try just a single iteration

			if (done) exit
		enddo
			
	enddo
	deallocate(a, x, rhs)
	
enddo
!$omp end parallel do

end subroutine

end module

