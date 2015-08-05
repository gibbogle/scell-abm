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
integer :: ix, iy, iz, maxnz, ichemo
real(REAL_KIND) :: C0
character*(10) :: emapfile, bmapfile
real(REAL_KIND), pointer :: Cprev(:,:,:), Fprev(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
logical :: ok

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
call make_sparse_emap(emapfile,.true.)

bmapfile = ''
write(bmapfile,'(a,i2.0,a)') 'bmap',NXB,'.dat'
write(nflog,*) 'setup_react_diff: ',bmapfile
call make_sparse_map(bmapfile,.false.)

do ichemo = 1,2
	Cextra_all(:,:,:,ichemo) = chemo(ichemo)%bdry_conc
	Caverage(:,:,:,ichemo) = chemo(ichemo)%bdry_conc
enddo

call update_Cex_Cin(.true.,0.0d0)	! initialise with SS values
call getF_SS	! estimates all fine grid pt fluxes Cflux(:,:,:,:) from Cextra(:,:,:,:)
! Sum fine grid fluxes to initialise Fcurr_b, Fprev_b

do ichemo = 1,2
	C0 = chemo(ichemo)%bdry_conc
	Cprev => chemo(ichemo)%Cprev
	Fprev => chemo(ichemo)%Fprev
!	Fcurr => chemo(ichemo)%Fcurr
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
	Cave_b => chemo(ichemo)%Cave_b
	
	Cave_b = C0
	Cprev_b = C0
	Cprev = C0
	Fprev = Cflux(:,:,:,ichemo)
!	Fcurr = Fprev
	call makeF_b(Fprev_b,Fprev,DELTA_T)
	Fcurr_b = Fprev_b
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! The flux values on the coarse grid are derived from the values on the fine grid by
! appropriate summation.  It's important that the total flux is consistent.
!-------------------------------------------------------------------------------------------
subroutine makeF_b(F_b,F,dt)
integer, parameter :: dN = NRF/2
real(REAL_KIND) :: F_b(:,:,:), F(:,:,:), dt
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
!						if (izb == 5 .and. ixb == xb0 .and. F(ix,iy,iz) /= 0) then
!							write(*,'(6i6,f8.4,e12.3)') ixb,iyb,izb,ix,iy,iz,wt(dx,dy,dz),F(ix,iy,iz)
!						endif
						F_b(ixb,iyb,izb) = F_b(ixb,iyb,izb) + wt(idx,idy,idz)*F(ix,iy,iz)
						Fsum_b = Fsum_b + wt(idx,idy,idz)*F(ix,iy,iz)
!						if (ichemo_curr == OXYGEN .and. istep >= 850 .and. ixb==18.and.iyb==18.and.izb==5) then
!							write(*,'(4i6,4e12.3)') istep,ix,iy,iz,wt(idx,idy,idz),F(ix,iy,iz),wt(idx,idy,idz)*F(ix,iy,iz),Cextra(ix,iy,iz,OXYGEN)
!						endif
					enddo
				enddo
			enddo
		enddo
	enddo
enddo	
!write(*,'(a,e12.3)') 'Fsum_b: ',Fsum_b
!do izb = 1,zb1
!	write(*,*) 'izb: ',izb
!	do ixb = xb0-dxb,xb0+dxb
!		write(*,'(10e12.3)') (F_b(ixb,iyb,izb),iyb=yb0-dyb,yb0+dyb)
!	enddo
!enddo
!stop
!do ixb = 1,NXB
!	do iyb = 1,NYB
!		do izb = 1,NZB
!			if (F_b(ixb,iyb,izb) /= 0) then
!				write(*,'(a,3i4,e12.3)') 'F_b: ',ixb,iyb,izb,F_b(ixb,iyb,izb)
!			endif
!		enddo
!	enddo
!enddo
! Try limiting F_b - not thought through, transport by diffusion is the main source of O2 in a time step
!do ichemo = 1,NCONST
!	do ixb = 1,NXB
!		do iyb = 1,NYB
!			do izb = 1,NZB
!				if (F_b(ixb,iyb,izb) /= 0) then
!					Fmin = 0.5*chemo(ichemo)%Cave_b(ixb,iyb,izb)*dxb3/dt
!					F_b(ixb,iyb,izb) = min(Fmin,F_b(ixb,iyb,izb))
!				endif
!			enddo
!		enddo
!	enddo
!enddo
end subroutine

!-------------------------------------------------------------------------------------------
! This version is for the embedded grid
!-------------------------------------------------------------------------------------------
subroutine make_csr_SS(a, ichemo, Cave, Fcurr, rhs)
integer :: ichemo
real(REAL_KIND) :: a(:), Cave(:,:,:), Fcurr(:,:,:), rhs(:)
integer :: k, ix, iy, iz, krow, kcol, nc
integer :: nc_max = 10	! just a wild guess but not a bad one
real(REAL_KIND) :: Kd, Kr, Vex, Cbdry
logical, save :: first = .true.

!Kr = dx*dx/Kd

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
			Kd = chemo(ichemo)%medium_diff_coef
			Kd = Kd*(1 - chemo(ichemo)%diff_reduction_factor*min(nc,nc_max)/nc_max)
			Kr = 1/(dxf*Kd)
			rhs(krow) = -Kr*Fcurr(ix,iy,iz)
		enddo
	enddo
enddo
!write(*,*) 'rhs: '
!write(*,'(5e12.3)') rhs

!Cbdry = chemo(ichemo)%bdry_conc
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
!write(*,*) 'adjusted rhs: '
!write(*,'(5e12.3)') rhs
first = .false.
end subroutine

!-------------------------------------------------------------------------------------------
! This is the version for the coarse grid.
! Use variable numbering (ix,iy) -> k = (ix-1)*NY + iy
! Since the equation is now dM/dt = ... need to derive different expression for Kr
! Need to add treatment of top boundary for O2
!-------------------------------------------------------------------------------------------
subroutine make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs)
integer :: ichemo
real(REAL_KIND) :: dt, a_b(:), Cave_b(:,:,:), Cprev_b(:,:,:), Fcurr_b(:,:,:), Fprev_b(:,:,:), rhs(:)
integer :: ixb, iyb, izb, k, i, krow, kcol, nc
real(REAL_KIND) :: Kd, Kr, Cbdry, Fsum
integer, parameter :: m = 3

Kd = chemo(ichemo)%medium_diff_coef
Fsum = 0
krow = 0
do k = 1,nnz_b
	if (k == ia_b(krow+1)) krow = krow+1
	kcol = ja_b(k)
	if (amap_b(k,0) == 2*m) then
		Kr = dxb*dxb/Kd
		a_b(k) = 3*Kr/(2*dt) + 2*m	! ... note that Kd should depend on (ixb,iyb,izb), ultimately on # of cells
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
			Kr = dxb*dxb/Kd	
			rhs(krow) = Kr*((-2*Fcurr_b(ixb,iyb,izb) + Fprev_b(ixb,iyb,izb))/dxb3 + (1./(2*dt))*(4*Cave_b(ixb,iyb,izb) - Cprev_b(ixb,iyb,izb)))
!			rhs(krow) = Kr*(-0.1/nrow_b + (1./(2*dt))*(4*Cave_b(ixb,iyb,izb) - Cprev_b(ixb,iyb,izb)))
!			if (Fcurr_b(ixb,iyb,izb) > 0) then
!				Fsum = Fsum + Fcurr_b(ixb,iyb,izb)/dxb3
!				write(*,*) Fcurr_b(ixb,iyb,izb)/dxb3
!			endif
!			if (ixb == NXB/2 .and. iyb == NYB/2) then
!				write(*,'(i6,6e12.3)') krow,Kr,Fcurr_b(ixb,iyb,izb),Fprev_b(ixb,iyb,izb),Cave_b(ixb,iyb,izb),Cprev_b(ixb,iyb,izb),rhs(krow)
!			endif
		enddo
	enddo
enddo
!!$omp end parallel do
!write(*,*) 'Fsum: ',Fsum
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

!write(*,*) 'rhs: '
!write(*,'(5e12.3)') rhs(1:nrow_b)
!stop
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine gmres_solver(dt)
real(REAL_KIND) :: dt

!call make_csr_b(dt, rhs)
end subroutine

!-------------------------------------------------------------------------------------------
! For now assume that the rate of consumption (or production) of the constituent is not
! dependent on cell volume V.
! Also consumption of GLUCOSE does not depend on O2.
! The function returns the mass rate in mumol/s
!-------------------------------------------------------------------------------------------
function fmetab(ic, C, V) result (f)
integer :: ic
real(REAL_KIND) :: C, V, f
real(REAL_KIND) :: C0
integer :: N

if (ic == OXYGEN .or. ic == GLUCOSE) then
	N = chemo(ic)%Hill_N
	C0 = chemo(ic)%MM_C0
	f = chemo(ic)%max_cell_rate*C**N/(C0**N + C**N)
endif
end function

!-------------------------------------------------------------------------------------------
! If there is a cell, the flux F includes the rate of transport through the cell membrane
! into the extracellular space.
! Decay is in common.
! Note that in this simplified version we will only let ncells = 0 or 1.
!-------------------------------------------------------------------------------------------
!subroutine getF(ichemo,Cin,Cex,F)
!integer :: ichemo
!real(REAL_KIND) :: Cin(:,:,:), Cave(:,:,:), F(:,:,:)
!real(REAL_KIND) :: C, Vex, Fsum
!real(REAL_KIND) :: Kin, Kout, decayrate_ex
!integer :: ix, iy, iz
!
!write(*,*) 'Need to check getF!'
!stop
!Kin = chemo(ichemo)%membrane_diff_in
!Kout = chemo(ichemo)%membrane_diff_out
!decayrate_ex = 0
!!Fsum = 0
!do ix = 1,NX
!do iy = 1,NY
!do iz = 1,NZ
!	C = Cave(ix,iy,iz)
!	if (ngcells(ix,iy,iz) > 0) then
!		Vex = dx3 - Vin(ix,iy,iz)
!		F(ix,iy,iz) = ((Kout*Cin(ix,iy,iz) - Kin*C) - C*Vex*decayrate_ex + C*dVindt(ix,iy,iz))/Vex
!	else
!		Vex = dx3
!		F(ix,iy,iz) = - C*Vex*decayrate_ex
!	endif
!!	Fsum = Fsum + F(ix,iy,iz)
!enddo
!enddo
!enddo
!!write(*,*) 'Fsum: ',Fsum
!end subroutine

!-------------------------------------------------------------------------------------------
! Update Cextra(:,:,:,:) and cp%Cin(:) from Caverage(:,:,:,:), taking account of 
! cell locations in nbr_list
! loop ...
! From Cex estimate nearby cell concentrations cp%Cex and cp%Cin
! Improve estimate of Cex using cp%Cex, cp%Cin and cp%V for nearby cells
! ... until convergence
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin(SS,dt)
logical :: SS
real(REAL_KIND) :: dt
integer :: ichemo, kcell, ix, iy, iz, idx, idy, idz, ixx, iyy, izz, nc, k, it
real(REAL_KIND) :: x1, x2, y1, y2, z1, z2, xmax, ymax, zmax, c(3), Vin, dCdt
real(REAL_KIND) :: s, r, frac(3), dV, vsum, cvsum(NCONST), Kin(NCONST), Kout(NCONST)
logical, allocatable :: skip(:,:,:)
type(cell_type), pointer :: cp

!write(*,*) 'update_Cex_Cin'
Cextra_all = Caverage
goto 10
allocate(skip(NX,NY,NZ))
skip = .false.

! start loop...

do it = 1,1
! Interpolate cp%Cex and infer cp%Cin
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	call extra_concs(kcell,	cp%Cex)
	do ichemo = 1,2
		cp%Cin(ichemo) = getCin_SS(ichemo,cp%V,cp%Cex(ichemo))
	enddo
!	write(*,'(i6,4e14.5)') kcell,cp%Cin(1:2),cp%Cex(1:2)
enddo
! Need to find cells contributing to the volume associated with each grid pt., and
! each one's volume contribution.
! The boundaries of the grid pt (ix,iy,iz) "cube of influence" (dotted lines) are:
!	x in (ix-1)*dxf - dxf/2, (ix-1)*dxf + dxf/2
!	y in (iy-1)*dxf - dyf/2, (iy-1)*dyf + dyf/2
!	z in (iz-1)*dxf - dzf/2, (iz-1)*dzf + dzf/2
! (respecting grid boundaries)
! The candidate cells are those in the 8 surrounding grid cubes
! A cell is assumed to contribute if it overlaps the cube of influence when represented by
! an equivalent cube (i.e. with s^3 = cell volume).
xmax = (NX-1)*dxf
ymax = (NY-1)*dxf
zmax = (NZ-1)*dxf
do ix = 1,NX
	x1 = max(0.0,(ix-1.5)*dxf)
	x2 = min(xmax,x1 + dxf)
	do iy = 1,NY
		y1 = max(0.0,(iy-1.5)*dxf)
		y2 = min(ymax,y1 + dxf)
		do iz = 1,NZ
			z1 = max(0.0,(iz-1.5)*dxf)
			z2 = min(zmax,z1 + dxf)
			if (skip(ix,iy,iz)) cycle
			vsum = 0
			cvsum = 0
			
			do idx = -1,0
				ixx = ix + idx
				if (ixx < 1 .or. ixx == NX) cycle
				do idy = -1,0
					iyy = iy + idy
					if (iyy < 1 .or. iyy == NY) cycle
					do idz = -1,0
						izz = iz + idz
						if (izz < 1 .or. izz == NZ) cycle
						nc = grid(ixx,iyy,izz)%nc
						if (nc == 0) cycle
						do k = 1,nc
							kcell = grid(ixx,iyy,izz)%cell(k)
							cp => cell_list(kcell)
							if (cp%state == DEAD) cycle
							if (cp%nspheres == 1) then
								c = cp%centre(:,1)
							else
								c = 0.5*(cp%centre(:,1) + cp%centre(:,2))
							endif
							s = cp%V**(1./3.)
							r = s/2
							if ((c(1)+r < x1 .or. c(1)-r > x2) .or. (c(2)+r < y1 .or. c(2)-r > y2) .or. (c(3)+r < z1 .or. c(3)-r > z2)) cycle
							! There is overlap!
							if (c(1)-r < x1) then		! cut by x=x1
								frac(1) = (c(1)+r-x1)/s
							elseif (c(1)+r > x2) then	! 
								frac(1) = (x2-c(1)+r)/s
							else
								frac(1) = 1
							endif
							if (c(2)-r < y1) then
								frac(2) = (c(2)+r-y1)/s
							elseif (c(2)+r > y2) then
								frac(2) = (y2-c(2)+r)/s
							else
								frac(2) = 1
							endif
							if (c(3)-r < z1) then
								frac(3) = (c(3)+r-z1)/s
							elseif (c(3)+r > z2) then
								frac(3) = (z2-c(3)+r)/s
							else
								frac(3) = 1
							endif
							dV = cp%V*frac(1)*frac(2)*frac(3)	! this is the overlap volume
!							write(*,'(a,4e12.3)') 'frac: ',frac,frac(1)*frac(2)*frac(3)
!							write(*,'(4i4,3e12.3)') ixx,iyy,izz,kcell,dV/cp%V,cp%Cin(1:2)
							vsum = vsum + dV
							cvsum = cvsum + dV*cp%Cin
						enddo
					enddo
				enddo
			enddo
			if (vsum == 0) then
				skip(ix,iy,iz) = .true.
			else
				Cextra_all(ix,iy,iz,:) = (Caverage(ix,iy,iz,:)*dx3 - cvsum(:))/(dx3-vsum)
!				if (ix==17.and.iy==17.and.iz==17) then
!					write(*,'(3i4,5e12.3)') ix,iy,iz,Caverage(ix,iy,iz,1),cvsum(1),dx3,vsum,Cextra(ix,iy,iz,1)
!				endif
!				do ichemo = 1,NCONST
!					if (Cextra_all(ix,iy,iz,ichemo) < 0) then
!						write(*,'(a,4i4,5e12.3)') 'Cextra < 0: ',ix,iy,iz,ichemo,Caverage(ix,iy,iz,ichemo),dx3,vsum,cvsum(ichemo),Cextra(ix,iy,iz,ichemo)
!						Cextra(ix,iy,iz,ichemo) = 0
!						stop
!					endif 
!				enddo
			endif
		enddo
	enddo
enddo
enddo
! end loop

deallocate(skip)

10	continue

Kin = chemo(1:NCONST)%membrane_diff_in
Kout = chemo(1:NCONST)%membrane_diff_out
if (SS) then
	do kcell = 1,nlist
		cp => cell_list(kcell)
		if (cp%state == DEAD) cycle
		call extra_concs(kcell,	cp%Cex)
		do ichemo = 1,2
			cp%Cin(ichemo) = getCin_SS(ichemo,cp%V,cp%Cex(ichemo))
			dCdt = (Kin(ichemo)*cp%Cex(ichemo) - Kout(ichemo)*cp%Cin(ichemo) - fmetab(ichemo,cp%Cin(ichemo),cp%v))/cp%v
!			if (ichemo == 1 .and. kcell == 1) then
!				write(*,'(a,5e12.3)') 'SS: ',cp%Cex(1),cp%Cin(1),fmetab(ichemo,cp%Cin(ichemo),cp%V), cp%V
!				write(*,'(5e14.5)') Kin(ichemo)*cp%Cex(ichemo), - Kout(ichemo)*cp%Cin(ichemo), - fmetab(ichemo,cp%Cin(ichemo),cp%V),dCdt*cp%V, dCdt
!			endif
		enddo
	enddo
else
	do kcell = 1,nlist
		cp => cell_list(kcell)
		if (cp%state == DEAD) cycle
		call extra_concs(kcell,	cp%Cex)
		Vin = cp%V
		do ichemo = 1,2
			dCdt = (Kin(ichemo)*cp%Cex(ichemo) - Kout(ichemo)*cp%Cin(ichemo) - fmetab(ichemo,cp%Cin(ichemo),Vin))/Vin
!			if (ichemo == 1 .and. kcell == 1) then
!				write(*,'(a,5e12.3)') 'update_Cex_Cin: ',cp%Cex(1),cp%Cin(1),fmetab(ichemo,cp%Cin(ichemo),Vin), Vin
!				write(*,'(5e14.5)') Kin(ichemo)*cp%Cex(ichemo), -Kout(ichemo)*cp%Cin(ichemo), - fmetab(ichemo,cp%Cin(ichemo),Vin),dCdt*Vin, dCdt
!			endif
			cp%Cin(ichemo) = cp%Cin(ichemo) + dCdt*dt
		enddo
	enddo
	stop
endif							
end subroutine

!-------------------------------------------------------------------------------------------
! Update Cextra(:,:,:) and thereby cp%Cex(ichemo) and cp%Cin(ichemo) from Caverage(:,:,:,ichemo), 
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin_const(ichemo,dt)
integer :: ichemo
real(REAL_KIND) :: dt
integer :: kcell
real(REAL_KIND) :: Kin, Kout, dCdt, dMdt
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
	cp%Cin(ichemo) = getCin_SS(ichemo,cp%V,cp%Cex(ichemo))
	if (kcell == 1) then
		dMdt = -(Kin*cp%Cex(ichemo) - Kout*cp%Cin(ichemo))	! rate of uptake
		dCdt = (-dMdt - fmetab(ichemo,cp%Cin(ichemo),cp%V))/cp%V
!		write(*,'(a,5e12.3)') 'Cex,Cin,dC,dCdt: ',cp%Cex(ichemo),cp%Cin(ichemo),cp%Cex(ichemo)-cp%Cin(ichemo),dMdt,dCdt
	endif
enddo
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
!-------------------------------------------------------------------------------------------
subroutine getF_SS
!real(REAL_KIND) :: Cave(:,:,:), F(:,:,:)
real(REAL_KIND) :: Km, Cin, Fsum
real(REAL_KIND) :: Kin(NCONST), Kout(NCONST)
real(REAL_KIND), allocatable :: Cex(:,:,:)
type(cell_type), pointer :: cp
integer :: ix, iy, iz, kcell
real(REAL_KIND) :: C, Cthreshold = 0.01

Kin = chemo(:)%membrane_diff_in
Kout = chemo(:)%membrane_diff_out

! This is handled by update_Cex_Cin
! loop ...
! From Cex estimate nearby cell concentrations cp%Cex and cp%Cin
! Improve estimate of Cex using cp%Cex, cp%Cin and cp%V for nearby cells
! ... until convergence

! Compute SS cell fluxes cp%dMdt
Fsum = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cp%dMdt = Kin*cp%Cex - Kout*cp%Cin
!	write(*,'(a,i6,6e11.3)') 'getF_SS: ',kcell,cp%Cex(1),cp%Cin(1),cp%dMdt(1)
! Adjust for O2 only for now
!	C = cp%Cex(OXYGEN)
!	if (C < Cthreshold) then
!		cp%dMdt(OXYGEN) = (C/Cthreshold)*cp%dMdt(OXYGEN)
!	endif
enddo

! Estimate grid pt flux values F
call make_grid_flux1
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
!-------------------------------------------------------------------------------------------
subroutine getF_SS_const(ichemo, Cflux_const)
integer :: ichemo
real(REAL_KIND) :: Cflux_const(:,:,:)
real(REAL_KIND) :: Kin, Kout
integer :: kcell
type(cell_type), pointer :: cp

!write(*,*) 'getF_SS_const: ',ichemo,nlist
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out

! Compute SS cell fluxes cp%dMdt
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cp%dMdt(ichemo) = Kin*cp%Cex(ichemo) - Kout*cp%Cin(ichemo)
enddo

! Estimate grid pt flux values F
call make_grid_flux(ichemo,Cflux_const)
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine getF
integer :: ichemo, kcell
real(REAL_KIND) :: Km, Cin, Fsum
real(REAL_KIND) :: Kin(NCONST), Kout(NCONST)
real(REAL_KIND), allocatable :: Cex(:,:,:)
type(cell_type), pointer :: cp
integer :: ix, iy, iz
real(REAL_KIND) :: C, Cthreshold = 0.01

Kin = chemo(:)%membrane_diff_in
Kout = chemo(:)%membrane_diff_out

! Compute cell fluxes cp%dMdt
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cp%dMdt = Kin*cp%Cex - Kout*cp%Cin
! Adjust for O2 only for now
!	C = cp%Cex(OXYGEN)
!	if (C < Cthreshold) then
!		cp%dMdt(OXYGEN) = (C/Cthreshold)*cp%dMdt(OXYGEN)
!	endif
	write(*,'(a,i6,6e11.3)') 'getF: ',kcell,cp%Cex(1),cp%Cin(1),cp%dMdt(1)
enddo
! Estimate grid pt flux values F
call make_grid_flux1
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine make_grid_flux1
integer :: ic, ix, iy, iz, kcell, cnr(3,8), i
type(cell_type), pointer :: cp
real(REAL_KIND) :: alfa = 0.7

Cflux = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ix = cp%site(1)
	iy = cp%site(2)
	iz = cp%site(3)
	cnr(:,1) = [ix, iy, iz]
	cnr(:,2) = [ix, iy+1, iz]
	cnr(:,3) = [ix, iy, iz+1]
	cnr(:,4) = [ix, iy+1, iz+1]
	cnr(:,5) = [ix+1, iy, iz]
	cnr(:,6) = [ix+1, iy+1, iz]
	cnr(:,7) = [ix+1, iy, iz+1]
	cnr(:,8) = [ix+1, iy+1, iz+1]
	call grid_flux_contribution(kcell, cnr)
!	write(*,'(a,3i4,2e12.3)') 'grid_flux: ',ix,iy,iz,Cflux(ix,iy,iz,1:2)
enddo
Cflux(:,:,:,1:NCONST) = (1-alfa)*Cflux(:,:,:,1:NCONST) + alfa*Cflux_prev(:,:,:,1:NCONST)
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
!-------------------------------------------------------------------------------------------
function getCin_SS(ichemo, Vin, Cex) result(Cin)
integer :: ichemo
real(REAL_KIND) :: Vin, Cex, Cin, Cex_t, Cin_t
real(REAL_KIND) :: Kin, Kout, Kmax, VKdecay, C0, a, b, c, D, r(3)
integer :: i, n, isign

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
!if (Cex < 0) then
!	isign = -1
!	Cex_t = -Cex
!else
!	isign = 1
!	Cex_t = Cex
!endif
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
Kmax = chemo(ichemo)%max_cell_rate
VKdecay = Vin*chemo(ichemo)%decay_rate
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
!Cin = isign*Cin_t
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
!				asum(1) = asum(1) + Cave(ix,iy,iz)
!				write(*,'(3i6,f8.4)') ix,iy,iz,Cave(ix,iy,iz)
			enddo
		enddo
!		if (mod(istep,6) == 0) then
!			write(*,'(a,(10f7.4))') 'Top: ',Cave(:,NY/2,NZ)
!		endif
	enddo
enddo
!write(*,*) 'Average top face: ',asum(1)/nsum
!stop
end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine simulate_step1(dt)
real(REAL_KIND) :: dt
integer :: ix, iy, iz

istep = istep + 1
!write(*,'(a,i6,f8.1)') 'istep, hour: ',istep,istep*dt/3600
call diff_solver(dt)
! Currently not inferring Cex and Cin from Cave and Vin
!do iz = 1,NZ-1
!	do iy = 2,NY-1
!		do ix = 2,NX-1
!			if (ngcells(ix,iy,iz) > 0) then
!				Vin(ix,iy,iz) = Vin(ix,iy,iz) + dt*dVindt(ix,iy,iz)	! update Vin using dV/dt - constant for now
!			endif
!		enddo
!	enddo
!enddo
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
!-------------------------------------------------------------------------------------- 
subroutine diff_solver(dt)
real(REAL_KIND) :: dt
integer :: i, k, k1, ix, iy, iz, irow, icol, kc, ic, icc, it
integer :: ixb, iyb, izb
integer :: ichemo, ierr, nfill, iters, maxits, im_krylov
real(REAL_KIND) :: R, tol, tol_b, asum, t, Vex_curr, Vex_next, Vin_curr, Vin_next, dCsum, msum
real(REAL_KIND), allocatable :: x(:), rhs(:)
real(REAL_KIND), pointer :: Cave(:,:,:), Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
real(REAL_KIND), allocatable :: a(:), a_b(:)
real(REAL_KIND) :: alfa = 0.7
real(REAL_KIND) :: dCtol = 1.0e-4
integer :: ILUtype = 1
logical :: done
logical :: do_fine = .true.
logical :: use_const = .true.

nfill = 1	! Level of fill for ILUK preconditioner
tol = 1.0d-6
tol_b = 1.0d-6
im_krylov = 60	! dimension of Krylov subspace in (outer) FGMRES
maxits = 100

! Compute all steady-state grid point fluxes in advance from Cextra(:,:,:,:): Cflux(:,:,:,:)
! call getF_SS

! Parallel operation is currently not possible because Cflux and Cflux_prev have not been accounted for - see make_grid_flux
!$omp parallel do private(Cave, Cprev, Fprev, Fcurr, Cave_b, Cprev_b, Fprev_b, Fcurr_b, a, a_b, x, rhs, ix, iy, iz, ixb, iyb, izb, it, done, ichemo, icc, k, dCsum, msum, iters, ierr)
!do ichemo = 1,2 
do ic = 1,nchemo
	ichemo = chemomap(ic)
	if (chemo(ichemo)%constant) cycle
!	write(*,'(a,i2)') 'ichemo: ',ichemo
	ichemo_curr = ichemo
	icc = ichemo - 1
	allocate(rhs(nrow_b))
	allocate(x(nrow_b))
	allocate(a_b(MAX_CHEMO*nrow_b))
	Cave => Caverage(:,:,:,ichemo)
	Cprev => chemo(ichemo)%Cprev
	Fprev => chemo(ichemo)%Fprev
	Fcurr => Cflux(:,:,:,ichemo)
	Cave_b => chemo(ichemo)%Cave_b
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
!	write(*,*) 'Cave_b:'
!	write(*,'(5e15.6)') Cave_b(NXB/2,NYB/2,:)
		
	Fprev_b = Fcurr_b
	call makeF_b(Fcurr_b, Fcurr, dt)
!	if (ichemo == OXYGEN) write(*,'(a,3e12.3)') 'Cave_b,Fcurr_b(18,18,5): ',Cave_b(18,18,5),Fcurr_b(18,18,5)
	call make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs)		! coarse grid

	! Solve Cave_b(t+dt) on coarse grid
	!----------------------------------
	call itsol_create_matrix(icc,nrow_b,nnz_b,a_b,ja_b,ia_b,ierr)
	write(nflog,*) 'itsol_create_matrix: ierr: ',ierr
		
	if (ILUtype == 1) then
		call itsol_create_precond_ILUK(icc,nfill,ierr)
		write(nflog,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
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
!	write(nflog,*) 'call itsol_solve_fgmr_ILU'
	call itsol_solve_fgmr_ILU(icc, rhs, x, im_krylov, maxits, tol_b, iters, ierr)
	write(nflog,*) 'itsol_solve_fgmr_ILU: Cave_b: ierr, iters: ',ierr,iters
	call itsol_free_precond_ILU(icc, ierr)
!	write(nflog,*) 'did itsol_free_precond_ILU'
	call itsol_free_matrix(icc, ierr)
!	write(nflog,*) 'did itsol_free_matrix'

	Cprev_b = Cave_b
	msum = 0
	do izb = 1,NZB
		do iyb = 1,NYB
			do ixb = 1,NXB
				k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
				Cave_b(ixb,iyb,izb) = x(k)
				msum = msum + x(k)*dxb3		! this sums the mass of constituent in mumols
				if (x(k) < 0) then
					write(nflog,*) 'Cave_b < 0: ',ixb,iyb,izb,x(k)
				endif 
			enddo
		enddo
	enddo
!	write(*,'(a,e12.3)') 'mass msum: ',msum
!	write(nflog,*) 'Cave_b:'
!	write(nflog,'(5e15.6)') Cave_b(NXB/2,NYB/2,:)
	deallocate(a_b, x, rhs)
	
	if (.not.do_fine) then
		cycle
	else
	
		allocate(rhs(nrow))
		allocate(x(nrow))
		allocate(a(MAX_CHEMO*nrow))
		! interpolate Cave_b on fine grid boundary
		call interpolate_Cave(Cave, Cave_b)
		
		it = 0
		do
			it = it + 1
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
			!				if (x(k) < 0) then
			!					write(*,*) 'Cave < 0: ',ix,iy,iz,x(k)
			!					Cave(ix,iy,iz) = 0
			!					stop
			!				endif 
						enddo
					enddo
				enddo
			endif
			
!			write(nflog,*) 'call itsol_solve_fgmr_ILU'
			call itsol_solve_fgmr_ILU(icc,rhs, x, im_krylov, maxits, tol, iters, ierr)
			write(nflog,*) 'itsol_solve_fgmr_ILU: Cave: ierr, iters: ',ierr,iters
			call itsol_free_precond_ILU(icc,ierr)
			call itsol_free_matrix(icc,ierr)
			
			done = .false.
			! Convergence: the test for convergence needs to give equal or extra weight to the pts with low concentration.
			! This suggests using fractional change: (Cnew-Cprev)/Cprev
			dCsum = 0
			do iz = 1,NZ-1
				do iy = 2,NY-1
					do ix = 2,NX-1
						k = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
						Cave(ix,iy,iz) = (1-alfa)*x(k) + alfa*Cave(ix,iy,iz)
						dCsum = dCsum + abs((Cave(ix,iy,iz) - Cprev(ix,iy,iz))/Cave(ix,iy,iz))
					enddo
				enddo
			enddo
			dCsum = dCsum/(NX*NY*NZ)
			done = (abs(dCsum) < dCtol)
!			if (ichemo == OXYGEN) write(*,'(a,3e12.3)') 'Cave,Cflux(17,17,17): ',Cave(17,17,17),Cflux(17,17,17,1)
			
			Cprev = Cave
			if (use_const) then
				Cflux_prev(:,:,:,ichemo) = Fcurr
				call update_Cex_Cin_const(ichemo, dt)
				call getF_SS_const(ichemo, Fcurr)
			else
				Cflux_prev = Cflux
				call update_Cex_Cin(use_SS, dt)
				call getF_SS
			endif
			if (done) exit
		enddo
			
		deallocate(a, x, rhs)
	endif
	
enddo
!$omp end parallel do

! Need to update Cextra(:,:,:,:), cp%Cex(:) and cp%Cin(:) from Caverage(:,:,:,:), taking account of cell locations
!call update_Cex_Cin(use_SS, dt)

end subroutine

end module

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!program main
!use react_diff
!implicit none
!integer :: it, ncpu, nt = 100
!real(REAL_KIND) :: dt = 100
!real(8) :: t1, t2
!integer count_0, count_1, count_rate, count_max
!logical :: use_relax = .false.
!
!ncpu = 4
!call setup(ncpu)
!
!call system_clock(count_0, count_rate, count_max)
!t1 = count_0*1.0/count_rate
!istep = 0
!do it = 1,nt
!	call simulate_step(dt)
!enddo
!call system_clock(count_1, count_rate, count_max)
!t2 = count_1*1.0/count_rate
!write(*,*) 'secs: ',t2-t1
!stop
!
!
!do it = 1,nt
!	call gmres_solver(dt)
!!	write(*,'(i6,11f6.2)') it,Cave(1:11) 
!!	write(*,'(i6,11f6.2)') it,Cave(45:55,50,50)
!enddo
!call system_clock(count_1, count_rate, count_max)
!t2 = count_1*1.0/count_rate
!write(*,*) 'secs: ',t2-t1
!!do it =1,NX
!!	write(*,'(i4,e12.3)') it,Cave(it)
!!enddo
!end
