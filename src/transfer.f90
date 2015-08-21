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
	enddo
enddo		
fdata%ncells = nc
res = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_cell_constituents(nvars,cvar_index,nvarlen,name_array,narraylen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_cell_constituents
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
do iextra = 1,N_EXTRA
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
subroutine get_summary(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Ntagged_anoxia(MAX_CELLTYPES), Ntagged_drugA(MAX_CELLTYPES), Ntagged_drugB(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES)
integer :: Nviable(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: diam_um, vol_mm3_1000, nhypoxic(3), ngrowth(3), hypoxic_percent_10, growth_percent_10, necrotic_percent_10, &
    medium_oxygen_100, medium_glucose_100, medium_TPZ_drug_1000, medium_DNB_drug_1000
integer :: TNanoxia_dead, TNdrugA_dead, TNdrugB_dead, TNradiation_dead, &
           TNtagged_anoxia, TNtagged_drugA, TNtagged_drugB, TNtagged_radiation, Tplate_eff_10
real(REAL_KIND) :: vol_cm3, vol_mm3, hour, plate_eff(MAX_CELLTYPES)
! temporary
integer :: Nsites
real(REAL_KIND) :: Radius
    
Radius = getRadius()	! cm
Nsites = Ncells		! just temporary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ncells_type(1) = Ncells
!Ncells_type(2) = 0
hour = istep*DELTA_T/3600.
!vol_cm3 = Vsite_cm3*Nsites			! total volume in cm^3
vol_cm3 = (4*PI/3)*Radius**3
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_1000 = vol_mm3*1000			! 1000 * volume in mm^3
diam_um = 2*Radius*10000
Ntagged_anoxia = Nanoxia_tag - Nanoxia_dead				! number currently tagged by anoxia
!Ntagged_drugA = NdrugA_tag - NdrugA_dead				! number currently tagged by drugA
!Ntagged_drugB = NdrugB_tag - NdrugB_dead				! number currently tagged by drugB
Ntagged_drugA = 0
Ntagged_drugB = 0
Ntagged_radiation = Nradiation_tag - Nradiation_dead	! number currently tagged by radiation
call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000*nhypoxic(i_hypoxia_cutoff))/Ncells
call getGrowthCount(ngrowth)
growth_percent_10 = (1000*ngrowth(i_growth_cutoff))/Ncells
necrotic_percent_10 = (1000*(Nsites-Ncells))/Nsites
call getNviable(Nviable)
plate_eff = real(Nviable)/Ncells
plate_eff_10 = 1000*plate_eff
medium_oxygen_100 = 0		! need some sort of average concentrations
medium_glucose_100 = 0
medium_TPZ_drug_1000 = 0
medium_DNB_drug_1000 = 0
TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
!TNdrugA_dead = sum(NdrugA_dead(1:Ncelltypes))
!TNdrugB_dead = sum(NdrugB_dead(1:Ncelltypes))
TNdrugA_dead = 0
TNdrugB_dead = 0
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_drugA = sum(Ntagged_drugA(1:Ncelltypes))
TNtagged_drugB = sum(Ntagged_drugB(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
Tplate_eff_10 = sum(plate_eff_10(1:Ncelltypes))
summaryData(1:20) = [ istep, Ncells, TNanoxia_dead, TNdrugA_dead, TNdrugB_dead, TNradiation_dead, &
    TNtagged_anoxia, TNtagged_drugA, TNtagged_drugB, TNtagged_radiation, &
	diam_um, vol_mm3_1000, hypoxic_percent_10, growth_percent_10, necrotic_percent_10, Tplate_eff_10, &
	medium_oxygen_100, medium_glucose_100, medium_TPZ_drug_1000, medium_DNB_drug_1000 ]
write(nfres,'(2a12,i8,2e12.4,19i7,13e12.4)') gui_run_version, dll_run_version, istep, hour, vol_mm3, diam_um, Ncells_type(1:2), &
!    Nanoxia_dead(1:2), NdrugA_dead(1:2), NdrugB_dead(1:2), Nradiation_dead(1:2), &
    Nanoxia_dead(1:2), Ndrug_dead(1,1:2), Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_drugA(1:2), Ntagged_drugB(1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), ngrowth(:)/real(Ncells), (Nsites-Ncells)/real(Nsites), plate_eff(1:2), &
	chemo(OXYGEN)%medium_Cext, chemo(GLUCOSE)%medium_Cext, chemo(TPZ_DRUG)%medium_Cext, chemo(DNB_DRUG)%medium_Cext
		
!	call sum_dMdt(GLUCOSE)
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
	if (ichemo==GLUCOSE .and. cp%dMdt(ichemo) < 1.0e-15) then
		write(*,*) 'sum_dMdt: glucose small dMdt: ',kcell,cp%dMdt(ichemo)
		stop
	endif
enddo
total_dMdt = total_dMdt + asum
write(*,'(a,2i6,2e12.3)') 'sum_dMdt: ',ichemo,Nc,asum,total_dMdt*3600
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
!--------------------------------------------------------------------------------
subroutine get_concdata(nvars, ns, dx, ex_conc) BIND(C)
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
write(*,'(10i6)') ngc(:)
write(*,'(10f7.3)') ctemp(:,GLUCOSE)
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

end module