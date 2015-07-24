! To handle continuum modelling

module continuum

use global
!use fmotion
use nbr
!use reaction

implicit none


contains


!-----------------------------------------------------------------------------------------
! Grid-cell (i,j,k) contains all cells with centre in the range:
!    (ix-1)*DX <= x < ix*DX
!    (iy-1)*DX <= y < iy*DX
!    (iz-1)*DX <= z < iz*DX
!-----------------------------------------------------------------------------------------
subroutine setup_grid_cells
real(REAL_KIND) :: c(3)
integer :: kcell, ix, iy, iz
type(cell_type), pointer :: cp

grid(:,:,:)%nc = 0
do kcell = 1,nlist
!	write(*,*) 'setup_grid_cells: kcell: ',kcell
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	if (cp%nspheres == 1) then
		c = cp%centre(:,1)
	else
		c = 0.5*(cp%centre(:,1) + cp%centre(:,2))
	endif
	ix = c(1)/DELTA_X + 1
	iy = c(2)/DELTA_X + 1
	iz = c(3)/DELTA_X + 1
	grid(ix,iy,iz)%nc = grid(ix,iy,iz)%nc + 1
!	write(*,'(3e12.3,4i6)') c,ix,iy,iz,grid(ix,iy,iz)%nc
	grid(ix,iy,iz)%cell(grid(ix,iy,iz)%nc) = kcell
	cp%site = [ix,iy,iz]
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Estimates the extracellular concentrations at the location of the centre of a cell.
!-----------------------------------------------------------------------------------------
subroutine extra_concs(kcell,conc)
integer :: kcell
real(REAL_KIND) :: conc(:)
integer :: i, ix, iy, iz, ichemo
real(REAL_KIND) :: centre(3), alfa(3)
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: extra_concs: dead cell: ',kcell
	stop
endif
if (cp%nspheres == 1) then
	centre = cp%centre(:,1)
else
	centre = 0.5*(cp%centre(:,1) + cp%centre(:,2))
endif
ix = cp%site(1)
iy = cp%site(2)
iz = cp%site(3)
do i = 1,3
	alfa(i) = (centre(i) - (cp%site(i)-1)*DELTA_X)/DELTA_X
enddo
conc(:) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz,:)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz,:)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1,:)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1,:)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz,:)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz,:)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1,:)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1,:)
!do ichemo = 1,NCONST
!	if (conc(ichemo) < 0) then
!		write(*,'(a,4i4,4e12.3)') 'extra_concs: conc < 0: ',ichemo,ix,iy,iz,conc(ichemo),alfa
!		write(*,'(8e12.3)') Cextra(ix,iy,iz,ichemo),Cextra(ix,iy+1,iz,ichemo),Cextra(ix,iy,iz+1,ichemo),Cextra(ix,iy+1,iz+1,ichemo), &
!							Cextra(ix+1,iy,iz,ichemo),Cextra(ix+1,iy+1,iz,ichemo),Cextra(ix+1,iy,iz+1,ichemo),Cextra(ix+1,iy+1,iz+1,ichemo)
!		conc(ichemo) = 0
!	endif
!enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Flux contributions from a cell are accumulated at grid points given by cnr(:,:)
! This assumes that the cell flux rates have already been computed.
! Flux units: mumol/s
!-----------------------------------------------------------------------------------------
subroutine grid_flux_contribution(kcell,cnr)
integer :: kcell, cnr(3,8)
integer :: k
real(REAL_KIND) :: centre(3), gridpt(3), r(3), d(8), sum
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: grid_flux_contribution: dead cell: ',kcell
	stop
endif
if (cp%nspheres == 1) then
	centre = cp%centre(:,1)
else
	centre = 0.5*(cp%centre(:,1) + cp%centre(:,2))
endif
sum = 0
do k = 1,8
	gridpt(:) = (cnr(:,k)-1)*DELTA_X
	r = centre - gridpt
	d(k) = max(sqrt(dot_product(r,r)), small_d)
	sum = sum + 1/d(k)
enddo
! The grid flux weights are (1/d(k))/sum.  Note that dMdt > 0 for +ve flux into the cell, 
do k = 1,8
	Cflux(cnr(1,k),cnr(2,k),cnr(3,k),:) = Cflux(cnr(1,k),cnr(2,k),cnr(3,k),:) + cp%dMdt(:)*(1/d(k))/sum
enddo
end subroutine

#if 0
!-----------------------------------------------------------------------------------------
! useit(ix,1) is for ix-1
! useit(ix,2) is for ix+1
! ...
! Note: this dos not apply to O2 on one boundary!!
!-----------------------------------------------------------------------------------------
subroutine create_stencil
integer :: ix, iy, iz, i, j, n, del(3,3)
logical :: useit(3,2)

del = 0
do i = 1,3
	del(i,i) = 1
enddo
do ix = 1,NX
	do iy = 1,NY
		do iz = 1,NZ
			useit(1,1) = (ix /= 1)
			useit(1,2) = (ix /= NX)
			useit(2,1) = (iy /= 1)
			useit(2,2) = (iy /= NY)
			useit(3,1) = (iz /= 1)
			useit(3,2) = (iz /= NZ)
!			write(*,'(6L2)') useit
			n = 0
			do i = 1,3
				do j = 0,1
					if (useit(i,j+1)) then
						n = n+1
						stencil(ix,iy,iz)%indx(:,n) = [ix, iy, iz] + del(i,:)*(2*j-1)
					endif
				enddo
			enddo
			stencil(ix,iy,iz)%nindx = n
!			write(*,'(4i4)') ix,iy,iz,n
!			write(*,'(3i4)') (stencil(ix,iy,iz)%indx(:,j),j=1,n)
		enddo
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine solver1(dt)
real(REAL_KIND) :: dt
integer :: ic, ix, iy, iz, kcell, cnr(3,8)
type(cell_type), pointer :: cp
type(grid_type), pointer :: gp

! Can be parallelised
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	call extra_concs(kcell,	cp%Cex)
	call cell_reactor(cp, dt)
enddo

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
enddo

! Can be parallelised by striping in two passes
do iz = 1,NZ
	call diffuser(iz, dt)
enddo
	
end subroutine

!-----------------------------------------------------------------------------------------
! Need to keep units all consistent:
!   distance unit: um               converted to cm as required
!   time unit:     s
!   K diff unit:   cm^2/s			for consistency with spheroid-abm...
!   conc. unit:	   mM
!   flux unit:     mol/s			= 10^6 mM.cm^3/s
! OBSOLETE code
!-----------------------------------------------------------------------------------------
subroutine diffuser(iz, dt)
integer :: iz
real(REAL_KIND) :: dt
integer :: ix, iy, k, nc
real(REAL_KIND) :: vol_cm3, csum(NCONST), sumK(NCONST), dMdt(NCONST), Kdiff(NCONST,0:6), DX_cm
type(stencil_type) :: sten
type(grid_type), pointer :: gp

DX_cm = DELTA_X*1.0e4	! um -> cm
do iy = 1,NY
	do ix = 1,NX
		sten = stencil(ix,iy,iz)
		gp => grid(ix,iy,iz)
		nc = gp%nc
		Kdiff(:,0) = get_Kdiff(nc)
		csum = 0
		sumK = 0
		do k = 1,sten%nindx
			nc = grid(sten%indx(1,k),sten%indx(2,k),sten%indx(3,k))%nc
			Kdiff(:,k) = (get_Kdiff(nc) + Kdiff(:,0))/2
			sumK(:) = sumK(:) + Kdiff(:,k) 
			csum(:) = csum(:) + Kdiff(:,k)*Cextra(sten%indx(1,k),sten%indx(2,k),sten%indx(3,k),:)
		enddo
		if (iz == NZ) then		! O2 is a special case
			k = sten%nindx + 1
			nc = 0
			Kdiff(:,k) = (get_Kdiff(nc) + Kdiff(:,0))/2
			csum(OXYGEN) = csum(OXYGEN) + Kdiff(OXYGEN,k)*(C_O2_bdry - Cextra(ix,iy,iz,OXYGEN))
		endif
		csum(:) = csum(:) - sumK(:)*Cextra(ix,iy,iz,:)	! units here mM.cm^2/s, C is mM, Kdiff is in cm^2/s
		! Need to ensure correct scaling for unit consistency
		! flux = K.dA.dC/dx = K.dx.dC since dA/dx = dx.  Need dx in cm since K is cm^2/s
		dMdt = DX_cm*1.0e-6*csum + Cflux(ix,iy,iz,:)	! flux units mol/s (1 mM.cm^3 = 1.0e-6 mol)
		vol_cm3 = um3_cm3*get_volume(gp)
		Cextra(ix,iy,iz,:) = Cextra(ix,iy,iz,:) + dt*dMdt/vol_cm3
	enddo
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! The extracellular volume = (grid cell volume - sum of cell volumes), units um^3
!-----------------------------------------------------------------------------------------
function get_volume(gp) result(vol)
type(grid_type), pointer :: gp
real(REAL_KIND) :: vol
integer :: ic, kcell

vol = DELTA_X*DELTA_X*DELTA_X
do ic = 1,gp%nc
	kcell = gp%cell(ic)
	vol = vol - cell_list(kcell)%V
enddo
vol = max(vol,100.)	! maybe not needed
end function
#endif

end module

