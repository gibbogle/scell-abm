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
	grid(ix,iy,iz)%cell(grid(ix,iy,iz)%nc) = kcell
	cp%site = [ix,iy,iz]
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Estimates the extracellular concentrations at the location of the centre of a cell,
! from the field of Cextra at grid points
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
conc(:) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra_all(ix,iy,iz,:)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra_all(ix,iy+1,iz,:)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra_all(ix,iy+1,iz+1,:)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra_all(ix,iy,iz+1,:)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra_all(ix+1,iy,iz,:)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra_all(ix+1,iy+1,iz,:)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra_all(ix+1,iy+1,iz+1,:)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra_all(ix+1,iy,iz+1,:)
end subroutine

!-----------------------------------------------------------------------------------------
! Estimates the extracellular concentrations at the location of the centre of a cell,
! from the field of Cextra at grid points
!-----------------------------------------------------------------------------------------
subroutine extra_concs_const(kcell,Cextra,conc)
integer :: kcell
real(REAL_KIND) :: Cextra(:,:,:), conc
integer :: i, ix, iy, iz, ichemo
real(REAL_KIND) :: centre(3), alfa(3)
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: extra_concs_const: dead cell: ',kcell
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
conc = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine grid_interp(kcell,alfa)
integer :: kcell
real(REAL_KIND) :: alfa(3)
integer :: i, ix, iy, iz
real(REAL_KIND) :: centre(3)
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: extra_concs_const: dead cell: ',kcell
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
!conc = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
!        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
!        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
!        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
!        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
!        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
!        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)  &
!        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1)
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

!-----------------------------------------------------------------------------------------
! Flux contributions from a cell are accumulated at grid points given by cnr(:,:)
! This assumes that the cell flux rates have already been computed.
! Flux units: mumol/s
!-----------------------------------------------------------------------------------------
subroutine grid_flux_weights(kcell,cnr,wt)
integer :: kcell, cnr(3,8)
real(REAL_KIND) :: wt(8)
integer :: k
real(REAL_KIND) :: centre(3), gridpt(3), r(3), d(8), sum
type(cell_type), pointer :: cp

cp => cell_list(kcell)
if (cp%state == DEAD) then
	write(*,*) 'Error: grid_flux_weights: dead cell: ',kcell
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
!	Cflux(cnr(1,k),cnr(2,k),cnr(3,k),:) = Cflux(cnr(1,k),cnr(2,k),cnr(3,k),:) + cp%dMdt(:)*(1/d(k))/sum
	wt(k) = (1/d(k))/sum
enddo
end subroutine

end module

