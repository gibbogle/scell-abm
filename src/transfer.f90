! Transferring results to the GUI

module transfer

use global
use chemokine

implicit none

contains

!-----------------------------------------------------------------------------------------
! Note: the only requirement is that the array has the SAVE attribute
!-----------------------------------------------------------------------------------------
!subroutine test_array(array) bind(C)
subroutine test_array(narr1, narr2, cptr) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: test_array
use, intrinsic :: iso_c_binding
!TYPE, BIND(C) :: PASS
!   INTEGER (C_INT) :: lenf
!   TYPE (C_PTR)    :: f
!	real(c_double) :: f(*)
!END TYPE PASS
!TYPE (PASS), INTENT(INOUT) :: array
integer(c_int) :: narr1, narr2
TYPE (C_PTR)    :: cptr
!real(c_double), allocatable, save :: a(:,:)
real(c_double), save :: a(4,3)
integer :: i1, i2
integer :: n1=4, n2=3

!allocate(a(n1,n2))
do i1 = 1,n1
	do i2 = 1,n2
		a(i1,i2) = i1 + 10*i2
	enddo
enddo

narr1 = n1
narr2 = n2
cptr = c_loc(a)
write(*,'(4f6.0)') a

end subroutine

!-----------------------------------------------------------------------------------------
! Distances are still all cm
!-----------------------------------------------------------------------------------------
subroutine get_fielddata(axis, fraction, fdata, res) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fielddata
use, intrinsic :: iso_c_binding
integer(c_int) :: axis, res
real(c_double) :: fraction
type (fielddata_type)    :: fdata
type(celldata_type), save :: cdata(1000)
type(cell_type), pointer :: cp
real(REAL_KIND) :: x, y, z, c(3), r, rad
integer :: i, kcell, nc

fdata%NX = NX
fdata%NY = NY
fdata%NZ = NZ
fdata%NCONST = NCONST
fdata%DX = DELTA_X
fdata%conc_ptr = c_loc(Caverage)
fdata%cell_ptr = c_loc(cdata)
nc = 0
! start by fixing on central plane, ignoring fraction
if (axis == X_AXIS) then
	x = ((NX-1)/2)*DELTA_X		
elseif (axis == Y_AXIS) then
	y = ((NY-1)/2)*DELTA_X 
elseif (axis == Z_AXIS) then
	z = ((NZ-1)/2)*DELTA_X 
endif
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	do i = 1,cp%nspheres
		c = cp%centre(:,i)
		r = cp%radius(i)
		if (axis == X_AXIS) then
			if (c(1) + r < x .or. c(1) - r > x) cycle
			rad = sqrt(r**2-(c(1)-x)**2)
			nc = nc + 1
			cdata(nc)%radius = rad
			cdata(nc)%centre(1:2) = [c(2),c(3)]		! always use centre(1:2) to store the 2D coordinates
		elseif (axis == Y_AXIS) then
			if (c(2) + r < y .or. c(2) - r > y) cycle
			rad = sqrt(r**2-(c(2)-y)**2)
			nc = nc + 1
			cdata(nc)%radius = rad
			cdata(nc)%centre(1:2) = [c(1),c(3)]		! always use centre(1:2) to store the 2D coordinates
		elseif (axis == Z_AXIS) then
			if (c(3) + r < z .or. c(3) - r > z) cycle
			rad = sqrt(r**2-(c(3)-z)**2)
			nc = nc + 1
			cdata(nc)%radius = rad
			cdata(nc)%centre(1:2) = [c(1),c(2)]		! always use centre(1:2) to store the 2D coordinates
		endif
		if (cp%anoxia_tag) then
			cdata(nc)%status = 2	! tagged to die of anoxia
		elseif (cp%Cin(OXYGEN) < 4e-3) then
			cdata(nc)%status = 1	! radiobiological hypoxia
		elseif (cp%mitosis > 0) then
			cdata(nc)%status = 3	! in mitosis
		else
			cdata(nc)%status = 0
		endif
	enddo
enddo		
fdata%ncells = nc
res = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_constituents(nvars,cvar_index,nvarlen,name_array,narraylen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_constituents
use, intrinsic :: iso_c_binding
character(c_char) :: name_array(0:*)
integer(c_int) :: nvars, cvar_index(0:*), nvarlen, narraylen
integer :: ivar, k, ichemo
character*(24) :: name
character(c_char) :: c

write(nflog,*) 'get_cell_constituents'
nvarlen = 24
ivar = 0
k = ivar*nvarlen
cvar_index(ivar) = 0	! CFSE
name = 'CFSE'
call copyname(name,name_array(k),nvarlen)
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	ivar = ivar + 1
	k = ivar*nvarlen
	cvar_index(ivar) = ichemo
	name = chemo(ichemo)%name
	write(nflog,*) 'get_cell_constituents: ',ichemo,name
	call copyname(name,name_array(k),nvarlen)
enddo
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = GROWTH_RATE
name = 'Growth rate'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CELL_VOLUME
name = 'Cell volume'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = O2_BY_VOL
name = 'Cell O2xVol'
call copyname(name,name_array(k),nvarlen)
nvars = ivar + 1
write(nflog,*) 'did get_cell_constituents'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine copyname(name,name_array,n)
character*(*) :: name
character :: name_array(*)
integer :: n
integer :: k

do k = 1,n
	name_array(k) = name(k:k)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Get number of live cells
!-----------------------------------------------------------------------------------------
subroutine get_nFACS(n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nfacs
use, intrinsic :: iso_c_binding
integer(c_int) :: n
integer :: k, kcell

n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_FACS(facs_data) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_facs
use, intrinsic :: iso_c_binding
real(c_double) :: val, facs_data(*)
integer :: k, kcell, iextra, ichemo, ivar, nvars, var_index(32)
real(REAL_KIND) :: cfse_min

!write(nflog,*) 'get_FACS'
nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
do iextra = 1,N_EXTRA-1
	nvars = nvars + 1
	var_index(nvars) = MAX_CHEMO + iextra
enddo
cfse_min = 1.0e20
k = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
!	k = k+1
!	facs_data(k) = cell_list(kcell)%CFSE
!	k = k+1
!	facs_data(k) = cell_list(kcell)%dVdt
!	k = k+1
!	facs_data(k) = cell_list(kcell)%conc(OXYGEN)
!	if (cell_list(kcell)%conc(OXYGEN) <= 0.00001 .or. cell_list(kcell)%dVdt < 2.0e-6) then
!		write(nflog,'(2i6,2e12.3)') istep,kcell,cell_list(kcell)%dVdt,cell_list(kcell)%conc(OXYGEN)
!	endif
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
			cfse_min = min(val,cfse_min)
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%Cin(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = cell_list(kcell)%V
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%V*cell_list(kcell)%Cin(OXYGEN)
		endif
		k = k+1
		facs_data(k) = val
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! nhisto is the number of histogram boxes
! vmax(ivar) is the maximum value for variable ivar
! Probably need to adjust vmax() to a roundish value
!
! Compute 3 distributions: 1 = both cell types
!                          2 = type 1
!                          3 = type 2
! Stack three cases in vmax() and histo_data()
!-----------------------------------------------------------------------------------------
subroutine get_histo(nhisto, histo_data, vmin, vmax, histo_data_log, vmin_log, vmax_log) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_histo
use, intrinsic :: iso_c_binding
integer(c_int),value :: nhisto
real(c_double) :: vmin(*), vmax(*), histo_data(*)
real(c_double) :: vmin_log(*), vmax_log(*), histo_data_log(*)
real(REAL_KIND) :: val, val_log
integer :: n(3), i, ih, k, kcell, ict, ichemo, ivar, nvars, var_index(32)
integer,allocatable :: cnt(:,:,:)
real(REAL_KIND),allocatable :: dv(:,:), valmin(:,:), valmax(:,:)
integer,allocatable :: cnt_log(:,:,:)
real(REAL_KIND),allocatable :: dv_log(:,:), valmin_log(:,:), valmax_log(:,:)
!real(REAL_KIND) :: vmin_log(100), vmax_log(100)
!real(REAL_KIND),allocatable :: histo_data_log(:)

!write(nflog,*) 'get_histo'
nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
nvars = nvars + 1
var_index(nvars) = GROWTH_RATE
nvars = nvars + 1
var_index(nvars) = CELL_VOLUME
nvars = nvars + 1
var_index(nvars) = O2_BY_VOL

allocate(cnt(3,nvars,nhisto))
allocate(dv(3,nvars))
allocate(valmin(3,nvars))
allocate(valmax(3,nvars))
allocate(cnt_log(3,nvars,nhisto))
allocate(dv_log(3,nvars))
allocate(valmin_log(3,nvars))
allocate(valmax_log(3,nvars))
!allocate(histo_data_log(10000))
cnt = 0
valmin = 1.0e10
valmax = -1.0e10
cnt_log = 0
valmin_log = 1.0e10
valmax_log = -1.0e10
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ict = cell_list(kcell)%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%Cin(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cell_list(kcell)%V
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%Cin(OXYGEN)*Vcell_pL*cell_list(kcell)%V
		endif
		valmax(ict+1,ivar) = max(valmax(ict+1,ivar),val)	! cell type 1 or 2
		valmax(1,ivar) = max(valmax(1,ivar),val)			! both
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		valmin_log(ict+1,ivar) = min(valmin_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmin_log(1,ivar) = min(valmin_log(1,ivar),val_log)			! both
		valmax_log(ict+1,ivar) = max(valmax_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmax_log(1,ivar) = max(valmax_log(1,ivar),val_log)			! both
	enddo
	n(ict+1) = n(ict+1) + 1
	n(1) = n(1) + 1
enddo
do ivar = 1,nvars
	ichemo = var_index(ivar)
	if (ichemo == CELL_VOLUME) then
		valmin(:,ivar) = Vcell_pL*0.8
		valmin_log(:,ivar) = log10(Vcell_pL*0.8)
	else
		valmin(:,ivar) = 0
	endif
enddo

dv = (valmax - valmin)/nhisto
!write(nflog,*) 'dv'
!write(nflog,'(e12.3)') dv
dv_log = (valmax_log - valmin_log)/nhisto
!write(nflog,*) 'dv_log'
!write(nflog,'(e12.3)') dv_log
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ict = cell_list(kcell)%celltype
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%Cin(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL*cell_list(kcell)%V
		elseif (ichemo == O2_BY_VOL) then
			val = cell_list(kcell)%Cin(OXYGEN)*Vcell_pL*cell_list(kcell)%V
		endif
		k = (val-valmin(1,ivar))/dv(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(1,ivar,k) = cnt(1,ivar,k) + 1
		k = (val-valmin(ict+1,ivar))/dv(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(ict+1,ivar,k) = cnt(ict+1,ivar,k) + 1
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		k = (val_log-valmin_log(1,ivar))/dv_log(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(1,ivar,k) = cnt_log(1,ivar,k) + 1
		k = (val_log-valmin_log(ict+1,ivar))/dv_log(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(ict+1,ivar,k) = cnt_log(ict+1,ivar,k) + 1
	enddo
enddo

do i = 1,3
	if (n(i) == 0) then
		vmin((i-1)*nvars+1:i*nvars) = 0
		vmax((i-1)*nvars+1:i*nvars) = 0
		histo_data((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
		vmin_log((i-1)*nvars+1:i*nvars) = 0
		vmax_log((i-1)*nvars+1:i*nvars) = 0
		histo_data_log((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
	else
		do ivar = 1,nvars
			vmin((i-1)*nvars+ivar) = valmin(i,ivar)
			vmax((i-1)*nvars+ivar) = valmax(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data(k) = (100.*cnt(i,ivar,ih))/n(i)
			enddo
			vmin_log((i-1)*nvars+ivar) = valmin_log(i,ivar)
			vmax_log((i-1)*nvars+ivar) = valmax_log(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data_log(k) = (100.*cnt_log(i,ivar,ih))/n(i)
			enddo
		enddo
	endif
enddo
deallocate(cnt)
deallocate(dv)
deallocate(valmin)
deallocate(valmax)
deallocate(cnt_log)
deallocate(dv_log)
deallocate(valmin_log)
deallocate(valmax_log)
!deallocate(histo_data_log)
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of cell volume.
! nv is passed from the GUI
! Min divide volume = Vdivide0 - dVdivide
! therefore Vmin = (Vdivide0 - dVdivide)/2
! Max divide volume = Vmax = Vdivide0 + dVdivide
! dv = (Vmax - Vmin)/nv
! v0 = Vmin + dv/2
!--------------------------------------------------------------------------------
subroutine get_volprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_volprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax

!call logger('get_volprob')
Vmin = (Vdivide0 - dVdivide)/2
Vmax = Vdivide0 + dVdivide
dv = (Vmax - Vmin)/nv
v0 = Vmin + dv/2
prob(1:nv) = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%V
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of intracellular O2 level
! nv is passed from the GUI
!--------------------------------------------------------------------------------
subroutine get_oxyprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_oxyprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax, O2max

!call logger('get_oxyprob')
Vmin = 0
Vmax = chemo(OXYGEN)%bdry_conc
v0 = Vmin
dv = (Vmax - Vmin)/nv
prob(1:nv) = 0
n = 0
O2max = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%Cin(OXYGEN)
	O2max = max(O2max,v)
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getNecroticFraction(necrotic_fraction, vol_cm3)
real(REAL_KIND) :: necrotic_fraction, vol_cm3
real(REAL_KIND) :: necrotic_vol_cm3

! Until this is fixed
necrotic_vol_cm3 = 0
necrotic_fraction = 0
return
if (Ncells < 4000) then
	necrotic_vol_cm3 = 0
	necrotic_fraction = 0
else
	call getNecroticVolume(necrotic_vol_cm3)
	necrotic_fraction = necrotic_vol_cm3/vol_cm3
!	if (necrotic_fraction > 0) stop
endif
!write(*,'(a,4e12.3)') 'R,n_vol,vol,necrotic_fraction: ',Radius,necrotic_vol_cm3,vol_cm3,necrotic_fraction
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDiamVol(diam_cm,vol_cm3)
real(REAL_KIND) :: diam_cm, vol_cm3
real(REAL_KIND) :: radius, cntr(3), rng(3)
call getBlobCentreRange(cntr,rng,radius)
diam_cm = 2*radius
vol_cm3 = (4*PI/3)*radius**3
end subroutine

!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: diam_um, vol_mm3_1000, nhypoxic(3), ngrowth(3), hypoxic_percent_10, growth_percent_10, necrotic_percent_10, &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(2), npmm3
integer :: TNanoxia_dead, TNradiation_dead, TNdrug_dead(2),  &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_radiation, TNtagged_drug(2)
integer :: Tplate_eff_10   
real(REAL_KIND) :: diam_cm, vol_cm3, vol_mm3, hour, plate_eff(MAX_CELLTYPES), cmedium(MAX_CHEMO), necrotic_fraction

hour = istep*DELTA_T/3600.
call getDiamVol(diam_cm,vol_cm3)
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = diam_cm*10000
npmm3 = Ncells/vol_mm3

!Ntagged_anoxia(:) = Nanoxia_tag(:) - Nanoxia_dead(:)			! number currently tagged by anoxia
!Ntagged_radiation(:) = Nradiation_tag(:) - Nradiation_dead(:)	! number currently tagged by radiation
!Ntagged_drug(1,:) = Ndrug_tag(1,:) - Ndrug_dead(1,:)			! number currently tagged by drugA
!Ntagged_drug(2,:) = Ndrug_tag(2,:) - Ndrug_dead(2,:)			! number currently tagged by drugA
Ntagged_anoxia(:) = Nanoxia_tag(:) 			! number currently tagged by anoxia
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))

call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000*nhypoxic(i_hypoxia_cutoff))/Ncells
call getGrowthCount(ngrowth)
growth_percent_10 = (1000*ngrowth(i_growth_cutoff))/Ncells
call getNecroticFraction(necrotic_fraction,vol_cm3)
necrotic_percent_10 = 1000*necrotic_fraction
call getNviable(Nviable)
plate_eff = real(Nviable)/Ncells
plate_eff_10 = 1000*plate_eff
Tplate_eff_10 = sum(plate_eff_10(1:Ncelltypes))
call getMediumConc(cmedium)
medium_oxygen_1000 = cmedium(OXYGEN)*1000
medium_glucose_1000 = cmedium(GLUCOSE)*1000
medium_drug_1000(1) = cmedium(DRUG_A)*1000
medium_drug_1000(2) = cmedium(DRUG_B)*1000

summaryData(1:21) = [ istep, Ncells, TNanoxia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, growth_percent_10, necrotic_percent_10, Tplate_eff_10, &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(1), medium_drug_1000(2), npmm3 ]
write(nfres,'(2a12,i8,2e12.4,19i7,13e12.4)') gui_run_version, dll_run_version, istep, hour, vol_mm3, diam_um, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Ndrug_dead(1,1:2), Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_drug(1,1:2), Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), ngrowth(:)/real(Ncells), necrotic_fraction, plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A), cmedium(DRUG_B)

!write(nfres,'(a)') 'istep hour vol_mm3 diam_um Ncells(2) &
!Nanoxia_dead(2) NdrugA_dead(2) NdrugB_dead(2) Nradiation_dead(2) &
!Ntagged_anoxia(2) Ntagged_drugA(2) Ntagged_drugB(2) Ntagged_radiation(2) &
!f_hypox(3) f_growth(3) f_necrot plating_efficiency(2) &
!medium_oxygen medium_glucose medium_drugA medium_drugB'
		
call sum_dMdt(GLUCOSE)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_summary1(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Ntagged_anoxia(MAX_CELLTYPES), Ntagged_drug(2,MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES)
integer :: Nviable(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: diam_um, vol_mm3_1000, nhypoxic(3), ngrowth(3), hypoxic_percent_10, growth_percent_10, necrotic_percent_10, &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(2), npmm3
integer :: TNanoxia_dead, TNdrug_dead(2), TNradiation_dead, &
           TNtagged_anoxia, TNtagged_drug(2), TNtagged_radiation, Tplate_eff_10
real(REAL_KIND) :: diam_cm, vol_cm3, vol_mm3, hour, plate_eff(MAX_CELLTYPES), cmedium(MAX_CHEMO), necrotic_fraction
    
hour = istep*DELTA_T/3600.
call getDiamVol(diam_cm,vol_cm3)
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = diam_cm*10000
npmm3 = Ncells/vol_mm3
call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000*nhypoxic(i_hypoxia_cutoff))/Ncells
call getGrowthCount(ngrowth)
growth_percent_10 = (1000*ngrowth(i_growth_cutoff))/Ncells
call getNecroticFraction(necrotic_fraction, vol_cm3)
necrotic_percent_10 = 1000*necrotic_fraction	! need to estimate necrotic_fraction - how?
call getNviable(Nviable)
plate_eff = real(Nviable)/Ncells
plate_eff_10 = 1000*plate_eff

call getMediumConc(cmedium)		! need some sort of average concentrations
medium_oxygen_1000 = cmedium(OXYGEN)*1000
medium_glucose_1000 = cmedium(GLUCOSE)*1000
medium_drug_1000(1) = cmedium(DRUG_A)*1000
medium_drug_1000(2) = cmedium(DRUG_B)*1000

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNtagged_anoxia = sum(Nanoxia_tag(1:Ncelltypes))
TNtagged_drug(1) = sum(Ndrug_tag(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ndrug_tag(2,1:Ncelltypes))
TNtagged_radiation = sum(Nradiation_tag(1:Ncelltypes))
Tplate_eff_10 = sum(plate_eff_10(1:Ncelltypes))
summaryData(1:21) = [ istep, Ncells, TNanoxia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, growth_percent_10, necrotic_percent_10, Tplate_eff_10, &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(1), medium_drug_1000(2), npmm3 ]
write(nfres,'(2a12,i8,2e12.4,19i7,13e12.4)') gui_run_version, dll_run_version, istep, hour, vol_mm3, diam_um, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Ndrug_dead(1,1:2), Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ndrug_tag(1,1:2), Ndrug_tag(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), ngrowth(:)/real(Ncells), necrotic_fraction, plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A), cmedium(DRUG_B)
		
!	call sum_dMdt(GLUCOSE)
end subroutine



!--------------------------------------------------------------------------------
! Estimate average concentrations in the medium.
! In fact the average is over the coarse grid points that lie outside the blob.
!--------------------------------------------------------------------------------
subroutine getMediumConc(cmedium)
real(REAL_KIND) :: cmedium(:)
real(REAL_KIND) :: cntr(3), rng(3), radius, r2, d2
integer :: x, y, z
integer :: ixb, iyb, izb, nsum, ic, ichemo

call getBlobCentreRange(cntr,rng,radius)
r2 = radius**2
cmedium = 0
nsum = 0
do ixb = 1,NXB
	x = (ixb-1)*dxb
	do iyb = 1,NYB
		y = (iyb-1)*dxb
		do izb = 1,NZB
			z = (izb-1)*dxb
			d2 = (x - cntr(1))**2 + (y - cntr(2))**2 + (z - cntr(3))**2
			if (d2 > r2) then
				nsum = nsum + 1
				do ic = 1,nchemo
					ichemo = chemomap(ic)
					cmedium(ichemo) = cmedium(ichemo) + chemo(ichemo)%cave_b(ixb,iyb,izb)
				enddo
			endif
		enddo
	enddo
enddo
cmedium = cmedium/nsum
end subroutine

!--------------------------------------------------------------------------------
! Estimate average concentrations in the medium.
! In fact the average is over the coarse grid points, excluding the fine grid.
!--------------------------------------------------------------------------------
subroutine getMediumConc1(cmedium)
real(REAL_KIND) :: cmedium(:)
integer :: xb0, yb0, idxb, idyb, xb1, xb2, yb1, yb2, zb1, zb2
integer :: ixb, iyb, izb, nsum, ic, ichemo

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
cmedium = 0
nsum = 0
do ixb = 1,NXB
	do iyb = 1,NYB
		do izb = 1,NZB
			if (ixb >= xb1 .and. ixb <= xb2 .and. iyb >= yb1 .and. iyb <= yb2 .and. izb <= zb2) cycle
			nsum = nsum + 1
			do ic = 1,nchemo
				ichemo = chemomap(ic)
				cmedium(ichemo) = cmedium(ichemo) + chemo(ichemo)%cave_b(ixb,iyb,izb)
			enddo
		enddo
	enddo
enddo
cmedium = cmedium/nsum
end subroutine

!--------------------------------------------------------------------------------
! Compute total uptake rate for a constituent
!--------------------------------------------------------------------------------
subroutine sum_dMdt(ichemo)
integer :: ichemo
integer :: kcell, Nc
real(REAL_KIND) :: asum
type(cell_type), pointer :: cp

Nc = 0
asum = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	Nc = Nc + 1
	asum = asum + cp%dMdt(ichemo)
!	if (ichemo==GLUCOSE .and. cp%dMdt(ichemo) < 1.0e-15) then
!		write(*,*) 'sum_dMdt: glucose small dMdt: ',kcell,cp%dMdt(ichemo)
!		stop
!	endif
enddo
total_dMdt = total_dMdt + asum
!write(*,'(a,2i6,2e12.3)') 'sum_dMdt: ',ichemo,Nc,asum,total_dMdt*3600
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getNviable(Nviable)
integer :: Nviable(:)
integer :: kcell, ityp, idrug
logical :: tag

nviable = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%anoxia_tag .or. &
!	    cell_list(kcell)%drugA_tag .or. &
!	    cell_list(kcell)%drugB_tag .or. &
	    cell_list(kcell)%radiation_tag) cycle
	    tag = .false.
	    do idrug = 1,ndrugs_used
			if (cell_list(kcell)%drug_tag(idrug)) tag = .true.
		enddo
		if (tag) cycle
	    ityp = cell_list(kcell)%celltype
	Nviable(ityp) = Nviable(ityp) + 1
enddo	
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getHypoxicCount(nhypoxic)
integer :: nhypoxic(3)
integer :: kcell, i

nhypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%Cin(OXYGEN) < O2cutoff(i)) then
			nhypoxic(i) = nhypoxic(i) + 1
!			write(*,'(a,i6,2e12.3)') 'hypoxic: ',kcell,cell_list(kcell)%Cin(OXYGEN),O2cutoff(i)
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Need to compare growth rate with a fraction of average growth rate
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth)
integer :: ngrowth(3)
integer :: kcell, i, ityp
real(REAL_KIND) :: r_mean(2)

r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
ngrowth = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	do i = 1,3
		if (cell_list(kcell)%dVdt < growthcutoff(i)*r_mean(ityp)) ngrowth(i) = ngrowth(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
! The challenge is to find the appropriate projection/interpolation from the cells
! to the equi-spaced points on the line.
! Method:
! The (estimated) blob centre and range provides a line through the centre
! parallel to the X axis between the extremes of the blob in this direction.
! If the centre is c(3) and the range is rng(3), then the line is:
! x1 <= x <= x2  where x1 = c(1) - rng(1)/2, x2 = c(1) + rng(1)/2
! y = c(2)
! z = c(3)
! First the variable values need to be determined at all grid pts that enclose this line.
! In the case of extracellular concentrations, the grid pt values are used.
! For the non-concentration variables, i.e. those connected with cells, the
! grid pt values must be estimated by averaging over nearby cells.
! The relevant grid cells are those with "lower-left" corner pts given by:
! ixcnr = x1/dx + 1 ... x2/dx + 1
! iycnr = c(2)/dx + 1
! izcnr = c(3)/dx + 1
!--------------------------------------------------------------------------------
subroutine get_concdata(nvars, ns, dxc, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dxc, ex_conc(0:*)
real(REAL_KIND) ::  dx, x1, x2, y0, z0, x, y, z
real(REAL_KIND) :: cntr(3), rng(3), radius
real(REAL_KIND) :: alfax, alfay, alfaz
integer :: nxgpts, ix1, ix2, iy0, iz0, kx, ky, kz
integer :: ix, iy, iz, k, kp, ichemo, kcell, offset
integer :: gcell(100)
integer, allocatable ::  ngc(:)
real(REAL_KIND), allocatable :: ctemp(:,:,:,:)
type(cell_type), pointer :: cp

!call logger('get_concdata')
dx = DELTA_X

call getBlobCentreRange(cntr,rng,radius)

x1 = cntr(1) - rng(1)/2
x2 = cntr(1) + rng(1)/2
y0 = cntr(2)
z0 = cntr(3)
ix1 = x1/dx + 1
ix2 = x2/dx + 1
iy0 = y0/dx + 1
iz0 = z0/dx + 1
nxgpts = ix2 - ix1 + 2
nvars = 1 + MAX_CHEMO + N_EXTRA
allocate(ngc(NX))
allocate(ctemp(nxgpts,2,2,0:nvars-1))

! First calculate averages at the enclosing grid pts

ctemp = 0
ngc = 0
do kx = 1,nxgpts
	do ky = 1,2
		do kz = 1,2
			ix = ix1 + kx - 1
			iy = iy0 + ky - 1
			iz = iz0 + kz - 1
			call get_gridptcells(ix,iy,iz,ngc(ix),gcell)
			if (ngc(ix) > 0) then
				do k = 1,ngc(ix)
					kcell = gcell(k)
					cp => cell_list(kcell)
					ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
					ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
					ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
					ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
				enddo
				ctemp(kx,ky,kz,:) = ctemp(kx,ky,kz,:)/ngc(ix)
			endif
			do ichemo = 1,MAX_CHEMO
				ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
!				if (ichemo == OXYGEN .and. ctemp(kx,ky,kz,ichemo) > 0.18) then
!					write(*,*) 'bad Caverage: ',kx,ky,kz,ix,iy,iz,ctemp(kx,ky,kz,ichemo)
!					stop
!				endif
			enddo
		enddo
	enddo
enddo

ns = 20
dxc = (x2-x1)/(ns-1)

! Now need to interpolate values at ns points along the line from (x1,y0,z0) to (x2,y0,z0)
! alfax, alfay, alfaz each lies in (0,1)
y = y0
z = z0
iy = y/dx + 1
iz = z/dx + 1
alfay = (y - (iy-1)*dx)/dx
alfaz = (z - (iz-1)*dx)/dx
do kp = 1,ns
	x = x1 + (kp-1)*dxc
	ix = x/dx + 1
	alfax = (x - (ix-1)*dx)/dx
	kx = ix - ix1 + 1
!	if (kx >= nxgpts) then
!		write(*,*) 'bad kx: ',kx,nxgpts
!		write(*,'(a,i4,4f8.4,2i4)') 'kp,x,x1,x2,dxc,ix,kx: ',kp,x,x1,x2,dxc,ix,kx
!		stop
!	endif
	do ichemo = 0, nvars-1
		offset = ichemo*ns
		k = offset - 1 + kp
		ex_conc(k) = (1-alfax)*(1-alfay)*(1-alfaz)*ctemp(kx,1,1,ichemo) + &
					 (1-alfax)*alfay*(1-alfaz)*ctemp(kx,2,1,ichemo) + &
					 (1-alfax)*(1-alfay)*alfaz*ctemp(kx,1,2,ichemo) + &
					 (1-alfax)*alfay*alfaz*ctemp(kx,2,2,ichemo) + &
					 alfax*(1-alfay)*(1-alfaz)*ctemp(kx+1,1,1,ichemo) + &
					 alfax*alfay*(1-alfaz)*ctemp(kx+1,2,1,ichemo) + &
					 alfax*(1-alfay)*alfaz*ctemp(kx+1,1,2,ichemo) + &
					 alfax*alfay*alfaz*ctemp(kx+1,2,2,ichemo)
!		if (ichemo == OXYGEN) then
!			if (ex_conc(k) > 0.18) then
!				write(*,'(2i4,2f8.4)') kp,k,x,ex_conc(k)
!				 write(*,'(3i4,4f8.4)') kx,1,1,(1-alfax),(1-alfay),(1-alfaz),ctemp(kx,1,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx,2,1,(1-alfax),alfay,(1-alfaz),ctemp(kx,2,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx,1,2,(1-alfax),(1-alfay),alfaz,ctemp(kx,1,2,ichemo)
!				 write(*,'(3i4,4f8.4)') kx,2,2,(1-alfax),alfay,alfaz,ctemp(kx,2,2,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,1,1,alfax,(1-alfay),(1-alfaz),ctemp(kx+1,1,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,2,1,alfax,alfay,(1-alfaz),ctemp(kx+1,2,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,1,2,alfax,(1-alfay),alfaz,ctemp(kx+1,1,2,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,2,2,alfax,alfay,alfaz,ctemp(kx+1,2,2,ichemo)
!				stop
!			endif
!		endif
	enddo
enddo

deallocate(ctemp)
deallocate(ngc)
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
! The challenge is to find the appropriate projection/interpolation from the cells
! to the equi-spaced points on the line.
!--------------------------------------------------------------------------------
subroutine get_concdata1(nvars, ns, dx, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dx, ex_conc(0:*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
integer :: ix, iy, iz, i, ic, k, ichemo, kcell, x, y, z, x1, x2, offset
integer :: gcell(100)
integer, allocatable ::  ngc(:)
real(REAL_KIND), allocatable :: ctemp(:,:)
type(cell_type), pointer :: cp

!call logger('get_concdata')
nvars = 1 + MAX_CHEMO + N_EXTRA
allocate(ngc(NX))
allocate(ctemp(NX,0:nvars-1))

dx = DELTA_X

iy = NY/2+1
iz = NZ/2+1
!write(*,*) 'get_concdata: iy,iz: ',iy,iz

! First calculate averages at the grid pts

x1 = 0
x2 = 0
ctemp = 0
ngc = 0
do ix = 1,NX
	call get_gridptcells(ix,iy,iz,ngc(ix),gcell)
	if (ngc(ix) == 0) cycle
	if (x1 == 0) x1 = ix
	do k = 1,ngc(ix)
		kcell = gcell(k)
		cp => cell_list(kcell)
		ctemp(ix,0) = ctemp(ix,0) + cp%CFSE
		ctemp(ix,1:MAX_CHEMO) = ctemp(ix,1:MAX_CHEMO) + cp%Cex(1:MAX_CHEMO)
		ctemp(ix,GROWTH_RATE) = ctemp(ix,GROWTH_RATE) + cp%dVdt
		ctemp(ix,CELL_VOLUME) = ctemp(ix,CELL_VOLUME) + cp%V
		ctemp(ix,O2_BY_VOL) = ctemp(ix,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
	enddo
	ctemp(ix,:) = ctemp(ix,:)/ngc(ix)
	x2 = ix
enddo
do ix = x1,x2
	if (ngc(ix) == 0) then
		ctemp(ix,:) = ctemp(ix-1,:)
	endif
enddo
!write(*,'(10i6)') ngc(:)
!write(*,'(10f7.3)') ctemp(:,GLUCOSE)
ns = x2 - x1 + 1
!write(*,*) 'x1, x2, ns: ',x1,x2,ns

do ichemo = 0,nvars-1
	offset = ichemo*ns
	k = offset - 1
	do x = x1, x2
		k = k + 1
		ex_conc(k) = ctemp(x,ichemo)
!		if (ichemo == 1) write(*,*) x,k,ex_conc(k)
	enddo
enddo

deallocate(ctemp)
deallocate(ngc)

#if 0
rng(:,1) = Centre(:) - (adrop*Radius + 2)
rng(:,2) = Centre(:) + (adrop*Radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = Centre(2) + 0.5
z = Centre(3) + 0.5

! First need to establish the range of x that is inside the blob: (x1,x2)
x1 = 0
do x = rng(1,1),rng(1,2)
	kcell = occupancy(x,y,z)%indx(1)
	if (kcell <= OUTSIDE_TAG) then
		if (x1 == 0) then
			cycle
		else
			exit
		endif
	elseif (x1 == 0) then
		x1 = x
	endif
	x2 = x
enddo

ns = x2 - x1 + 1 
do ichemo = 0,nvars-1
	offset = ichemo*ns
	k = offset - 1
	do x = x1, x2
		k = k + 1
		kcell = occupancy(x,y,z)%indx(1)
		if (kcell <= OUTSIDE_TAG) then
			ex_conc(k) = 0
			cycle
		endif
		i = ODEdiff%ivar(x,y,z)
        if (ichemo == 0) then	! CFSE
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%CFSE
			else
				ex_conc(k) = 0
			endif
       elseif (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					ex_conc(k) = allstate(i,ichemo)	
				else
					ex_conc(k) = 0
				endif
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == GROWTH_RATE) then	
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%dVdt
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == CELL_VOLUME) then	
			if (kcell > 0) then
				ex_conc(k) = Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == O2_BY_VOL) then	
			if (kcell > 0) then
				ex_conc(k) = allstate(i,OXYGEN)*Vcell_pL*cell_list(kcell)%volume
			else
				ex_conc(k) = 0
			endif
		endif
    enddo
enddo
#endif
! Add concentrations at the two boundaries 
! At ns=1, at at ns=ns+1
!ns = ns+1
!do ic = 1,nvars
!	ichemo = ic - 1
!	if (ichemo == 0) then	! CFSE
!		k = ic	
!		ex_conc(k) = 0
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0	
!	elseif (ichemo <= MAX_CHEMO) then
!		k = ic
!		if (chemo(ichemo)%used) then
!			ex_conc(k) = BdryConc(ichemo,t_simulation)
!		else
!			ex_conc(k) = 0
!		endif      
!		k = (ns-1)*nvars + ic
!		if (chemo(ichemo)%used) then
!			ex_conc(k) = BdryConc(ichemo,t_simulation)
!		else
!			ex_conc(k) = 0
!		endif      
!    elseif (ichemo == MAX_CHEMO+1) then	! growth rate
!		k = ic
!		ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
!    elseif (ichemo == MAX_CHEMO+2) then	! cell volume
!		k = ic
!		ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
!    endif
!enddo
end subroutine

!--------------------------------------------------------------------------------
! Generates a list of cells with centres within the "cube of influence" of the
! grid pt (ix,iy,iz).
! The boundaries of the grid pt (ix,iy,iz) "cube of influence" (dotted lines) are:
!	x in (ix-1)*dxf - dxf/2, (ix-1)*dxf + dxf/2
!	y in (iy-1)*dxf - dyf/2, (iy-1)*dyf + dyf/2
!	z in (iz-1)*dxf - dzf/2, (iz-1)*dzf + dzf/2
! (respecting grid boundaries)
!--------------------------------------------------------------------------------
subroutine get_gridptcells(ix,iy,iz,ngc,gcell)
integer :: ix,iy,iz,ngc,gcell(:)
integer :: idx, idy, idz, ixx, iyy, izz, nc, k, kcell
real(REAL_KIND) :: xmax, ymax, zmax, x1, x2, y1, y2, z1, z2, c(3)
type(cell_type), pointer :: cp

xmax = (NX-1)*dxf
ymax = (NY-1)*dxf
zmax = (NZ-1)*dxf
x1 = max(0.0,(ix-1.5)*dxf)
x2 = min(xmax,x1 + dxf)
y1 = max(0.0,(iy-1.5)*dxf)
y2 = min(ymax,y1 + dxf)
z1 = max(0.0,(iz-1.5)*dxf)
z2 = min(zmax,z1 + dxf)
ngc = 0
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
				if ((c(1) < x1 .or. c(1) > x2) .or. (c(2) < y1 .or. c(2) > y2) .or. (c(3) < z1 .or. c(3) > z2)) cycle
				ngc = ngc + 1
				gcell(ngc) = kcell
			enddo
		enddo
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Estimate necrotic volume by determining the range in the three axis directions
!-----------------------------------------------------------------------------------------
subroutine getNecroticVolume(vol)
real(REAL_KIND) :: vol
real(REAL_KIND) :: rmin(3), rmax(3), rng(3), cblob(3), c(3), r, smin, smax, R3, Rn(3), Rdiv
integer :: kcell, k, js, ix, iy, iz, ixx, iyy, izz
logical :: overlap, hit1, hit2
type(cell_type), pointer :: cp
real(REAL_KIND) :: fac = 1.5

! First estimate blob centre location
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
cblob(:) = (rmin(:) + rmax(:))/2	! approximate blob centre
ix = cblob(1)/DELTA_X + 1
iy = cblob(2)/DELTA_X + 1
iz = cblob(3)/DELTA_X + 1
!write(*,*) 'getNecroticVolume:'
!write(*,'(a,3f8.4)') 'centre: ',cblob
!write(*,'(a,3i4)') 'centre gridsite: ',ix,iy,iz
! Check if any cells in this gridcell overlap with C(:)
! A cell is treated here as a cube with side = cell diameter
overlap = .false.
do k = 1,grid(ix,iy,iz)%nc
	kcell = grid(ix,iy,iz)%cell(k)
	cp => cell_list(kcell)
	do js = 1,cp%nspheres
		r = fac*cp%radius(js)
		c = cp%centre(:,js)
!		write(*,'(a,i6,3(2f8.4,2x))') 'kcell,bnds: ',kcell,c(1)-r,c(1)+r,c(2)-r,c(2)+r,c(3)-r,c(3)+r
		if ((cblob(1) >= c(1)-r .and. cblob(1) <= c(1)+r) .and. &
			(cblob(2) >= c(2)-r .and. cblob(2) <= c(2)+r) .and. &
			(cblob(3) >= c(3)-r .and. cblob(3) <= c(3)+r)) then
			overlap = .true.
			exit
		endif
	enddo
	if (overlap) exit
enddo
if (overlap) then
	vol = 0
	return
endif
! There is a gap in the centre, find the gap extent in each axis direction
! X direction.  
!write(*,*) 'seeking rmin, rmax in X direction'
smin = 1.0e10
smax = 0
hit1 = .false.
hit2 = .false.
do ixx = 1,NX-1
	do k = 1,grid(ixx,iy,iz)%nc
		kcell = grid(ixx,iy,iz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
			if ((cblob(2) > c(2)-r .and. cblob(2) < c(2)+r) .and. &
				(cblob(3) > c(3)-r .and. cblob(3) < c(3)+r)) then
!				write(*,'(3(2f8.4,2x))') c(1)-r,c(1)+r,c(2)-r,c(2)+r,c(3)-r,c(3)+r
				if (smax < c(1)+r .and. cblob(1) > c(1)+r) then
					smax = c(1)+r
!					write(*,'(a,f8.4)') 'smax: ',smax
					hit1 = .true.
				endif
				if (smin > c(1)-r .and. cblob(1) < c(1)-r) then
					smin = c(1)-r
!					write(*,'(a,f8.4)') 'smin: ',smin
					hit2 = .true.
				endif
			endif
		enddo
	enddo
enddo
if (hit1) then
	rmin(1) = smax
else
	write(*,*) 'Error: not hit1'
	stop
endif
if (hit2) then
	rmax(1) = smin
else
	write(*,*) 'Error: not hit2'
	stop
endif
! Y direction.  
!write(*,*) 'seeking rmin, rmax in Y direction'
smin = 1.0e10
smax = 0
hit1 = .false.
hit2 = .false.
do iyy = 1,NY-1
	do k = 1,grid(ix,iyy,iz)%nc
		kcell = grid(ix,iyy,iz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
			if ((cblob(1) > c(1)-r .and. cblob(1) < c(1)+r) .and. &
				(cblob(3) > c(3)-r .and. cblob(3) < c(3)+r)) then
!				write(*,'(3(2f8.4,2x))') c(1)-r,c(1)+r,c(2)-r,c(2)+r,c(3)-r,c(3)+r
				if (smax < c(2)+r .and. cblob(2) > c(2)+r) then
					smax = c(2)+r
!					write(*,'(a,f8.4)') 'smax: ',smax
					hit1 = .true.
				endif
				if (smin > c(2)-r .and. cblob(2) < c(2)-r) then
					smin = c(2)-r
!					write(*,'(a,f8.4)') 'smin: ',smin
					hit2 = .true.
				endif
			endif
		enddo
	enddo
enddo
if (hit1) then
	rmin(2) = smax
else
	write(*,*) 'Error: not hit1'
	stop
endif
if (hit2) then
	rmax(2) = smin
else
	write(*,*) 'Error: not hit2'
	stop
endif
!write(*,*) 'seeking rmin, rmax in Z direction'
smin = 1.0e10
smax = 0
hit1 = .false.
hit2 = .false.
do izz = 1,NZ-1
	do k = 1,grid(ix,iy,izz)%nc
		kcell = grid(ix,iy,izz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
			if ((cblob(2) > c(2)-r .and. cblob(2) < c(2)+r) .and. &
				(cblob(1) > c(1)-r .and. cblob(1) < c(1)+r)) then
!				write(*,'(3(2f8.4,2x))') c(1)-r,c(1)+r,c(2)-r,c(2)+r,c(3)-r,c(3)+r
				if (smax < c(3)+r .and. cblob(3) > c(3)+r) then
					smax = c(3)+r
!					write(*,'(a,f8.4)') 'smax: ',smax
					hit1 = .true.
				endif
				if (smin > c(3)-r .and. cblob(3) < c(3)-r) then
					smin = c(3)-r
!					write(*,'(a,f8.4)') 'smin: ',smin
					hit2 = .true.
				endif
			endif
		enddo
	enddo
enddo
if (hit1) then
	rmin(3) = smax
else
	write(*,*) 'Error: not hit1'
	stop
endif
if (hit2) then
	rmax(3) = smin
else
	write(*,*) 'Error: not hit2'
	stop
endif

#if 0
write(*,*) 'seeking rmax in X direction'
smin = 1.0e10
hit = .false.
do ixx = 1,NX-1
	do k = 1,grid(ixx,iy,iz)%nc
		kcell = grid(ixx,iy,iz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
			if ((cblob(2) > c(2)-r .and. cblob(2) < c(2)+r) .and. &
				(cblob(3) > c(3)-r .and. cblob(3) < c(3)+r)) then
				write(*,'(3(2f8.4,2x))') c(1)-r,c(1)+r,c(2)-r,c(2)+r,c(3)-r,c(3)+r
				if (smin > c(1)-r .and. cblob(1) < c(1)-r) then
					smin = c(1)-r
					write(*,'(a,f8.4)') 'smin: ',smin
					hit = .true.
				endif
			endif
		enddo
	enddo
enddo
if (hit) then
	rmax(1) = smin
else
	write(*,*) 'Error: not hit'
	stop
endif

! Y direction.  
write(*,*) 'Y direction'
do iyy = iy,1,-1
	hit = .false.
	smax = 0
	do k = 1,grid(ix,iyy,iz)%nc
		kcell = grid(ix,iyy,iz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
!			write(*,'(i2,f8.4)') kcell,c(2)+r
			if ((smax < c(2)+r .and. cblob(2) > c(2)+r) .and. &
				(cblob(1) > c(1)-r .and. cblob(1) < c(1)+r) .and. &
				(cblob(3) > c(3)-r .and. cblob(3) < c(3)+r)) then
				smax = c(2)+r
				hit = .true.
			endif
		enddo
	enddo
	if (hit) then
		rmin(2) = smax
		write(*,'(a,f8.4)') 'rmin(2): ',rmin(2)
		exit
	endif
enddo
if (.not.hit) then
	write(*,*) 'Error: not hit'
	stop
endif
do iyy = iy,NY-1
	hit = .false.
	smin = 1.0e10
	do k = 1,grid(ix,iyy,iz)%nc
		kcell = grid(ix,iyy,iz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
!			write(*,'(i2,f8.4)') kcell,c(2)-r
			if ((smin > c(2)-r .and. cblob(2) < c(2)-r) .and. &
				(cblob(1) > c(1)-r .and. cblob(1) < c(1)+r) .and. &
				(cblob(3) > c(3)-r .and. cblob(3) < c(3)+r)) then
				smin = c(2)-r
				hit = .true.
			endif
		enddo
	enddo
	if (hit) then
		rmax(2) = smin
		write(*,'(a,f8.4)') 'rmax(2): ',rmax(2)
		exit
	endif
enddo
if (.not.hit) then
	write(*,*) 'Error: not hit'
	stop
endif
! Z direction.  
write(*,*) 'Z direction'
do izz = iz,1,-1
	hit = .false.
	smax = 0
	do k = 1,grid(ix,iy,izz)%nc
		kcell = grid(ix,iy,izz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
			if ((smax < c(3)+r .and. cblob(3) > c(3)+r) .and. &
				(cblob(2) > c(2)-r .and. cblob(2) < c(2)+r) .and. &
				(cblob(1) > c(1)-r .and. cblob(1) < c(1)+r)) then
				smax = c(3)+r
				hit = .true.
			endif
		enddo
	enddo
	if (hit) then
		rmin(3) = smax
		write(*,'(a,f8.4)') 'rmin(3): ',rmin(3)
		exit
	endif
enddo
if (.not.hit) then
	write(*,*) 'Error: not hit'
	stop
endif
do izz = iz,NZ-1
	hit = .false.
	smin = 1.0e10
	do k = 1,grid(ix,iy,izz)%nc
		kcell = grid(ix,iy,izz)%cell(k)
		cp => cell_list(kcell)
		do js = 1,cp%nspheres
			r = fac*cp%radius(js)
			c = cp%centre(:,js)
			if ((cblob(2) > c(2)-r .and. cblob(2) < c(2)+r) .and. &
				(cblob(1) > c(1)-r .and. cblob(1) < c(1)+r)) then
				write(*,*) cblob(3),c(3)-r
				if ((smin > c(3)-r .and. cblob(3) < c(3)-r)) then
					smin = c(3)-r
					hit = .true.
				endif
			endif
		enddo
	enddo
	if (hit) then
		rmax(3) = smin
		write(*,'(a,f8.4)') 'rmax(3): ',rmax(3)
		exit
	endif
enddo
if (.not.hit) then
	write(*,*) 'Error: not hit'
	stop
endif
Rn = (rmax - rmin)/2
#endif

Rn = (rmax - rmin)/2
R3 = (1./3.)*(Rn(1)**3 + Rn(2)**3 + Rn(3)**3)
Rdiv = (3*Vdivide0/(4*PI))**(1./3.)
!write(*,'(a,e12.3)') 'Rdiv: ',Rdiv
!write(*,'(a,3e12.3)') 'Rn: ',Rn
if (Rn(1) < 0 .or. Rn(2) < 0 .or. Rn(3) < 0) then
	write(*,*) 'Rn < 0'
	stop
endif
if (R3 > 0) then
	R = R3**(1./3)
	if (R < 4*Rdiv) then
		R3 = 0
		vol = 0
	endif
elseif (R3 == 0) then
	R = 0
else
	write(*,*) 'Error: getNecroticVolume: R3 < 0: ',R3
	stop
endif
vol = (4./3.)*PI*R3 
!write(*,'(a,2e12.3)') 'getNecroticVolume: R, vol: ',R, vol
end subroutine

subroutine write_bdryconcs
integer :: bcells(3,2)
real(REAL_KIND) :: c(3,2)
integer :: ichemo, iaxis, iend, kcell

call getBdryCells(bcells)
do ichemo = 1,2
	do iaxis = 1,3
		do iend = 1,2
			kcell = bcells(iaxis,iend)
			c(iaxis,iend) = cell_list(kcell)%Cex(ichemo)
		enddo
	enddo
	write(nflog,'(a,6f8.4,a,f8.4)') chemo(ichemo)%name,c,'  ave: ',sum(c)/6
enddo
	
end subroutine

end module