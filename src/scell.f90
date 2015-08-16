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
use transfer
use cellstate

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
if (allocated(Cextra_all)) deallocate(Cextra_all)
if (allocated(Caverage)) deallocate(Caverage)
if (allocated(Cflux)) deallocate(Cflux)
if (allocated(Cflux_prev)) deallocate(Cflux_prev)
if (allocated(stencil)) deallocate(stencil)
call logger('did deallocation')
allocate(cell_list(MAX_NLIST))
allocate(grid(NX,NY,NZ))
allocate(perm_index(MAX_NLIST))
allocate(Cextra_all(NX,NY,NZ,NCONST))
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
integer :: i, idrug, nmetab, im, ichemo
integer :: Nmm3, itestcase, ictype, ishow_progeny
integer :: iuse_oxygen, iuse_glucose, iuse_tracer, iuse_drug, iuse_metab, iV_depend, iV_random
integer :: iuse_extra, iuse_relax, iuse_par_relax
real(REAL_KIND) :: days, percent, fluid_fraction, d_layer, sigma(MAX_CELLTYPES), Vsite_cm3, bdry_conc, spcrad_value
integer :: iuse_drop, iconstant, isaveprofiledata
logical :: use_metabolites
character*(12) :: drug_name
character*(1) :: numstr

ok = .true.

open(nfcell,file=inputfile,status='old')

read(nfcell,*) gui_run_version				! program run version number
read(nfcell,*) dll_run_version				! DLL run version number
read(nfcell,*) NX							! size of fine grid
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median(1)
read(nfcell,*) divide_time_shape(1)
read(nfcell,*) divide_time_median(2)
read(nfcell,*) divide_time_shape(2)
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) NXB							! size of coarse grid
read(nfcell,*) DELTA_X						! grid size (um)
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

do i = 1,2			! currently allowing for just two different drugs: 1 = TPZ-type, 2 = DNB-type
	read(nfcell,*) iuse_drug
	read(nfcell,'(a12)') drug_name
	read(nfcell,*) bdry_conc
	read(nfcell,*) iconstant
	read(nfcell,*) iuse_metab
	call getIndices(i, idrug, nmetab)
	if (idrug < 0 .and. iuse_drug /= 0) then
		write(logmsg,*) 'Unrecognized drug name: ',drug_name
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (idrug == 0) cycle
	chemo(idrug)%name = drug_name
	chemo(idrug)%used = (iuse_drug == 1)
	chemo(idrug)%bdry_conc = bdry_conc
	chemo(idrug)%constant = (iconstant == 1)
!	chemo(idrug)%decay = (idrug_decay == 1)
	if (chemo(idrug)%used) then
		use_metabolites = (iuse_metab == 1)
		write(nflog,*) 'drug: ',idrug,'  name: ',chemo(idrug)%name,' use_metabolites: ',use_metabolites
	else
		use_metabolites = .false.
	endif
!	chemo(imetab)%decay = (imetab_decay == 1)
	if (idrug == TPZ_DRUG) then
		TPZ%name = drug_name
		TPZ%nmetabolites = nmetab
		TPZ%use_metabolites = use_metabolites
		do im = 0,nmetab
			read(nfcell,*) TPZ%diff_coef(im)
			read(nfcell,*) TPZ%medium_diff_coef(im)
			read(nfcell,*) TPZ%membrane_diff_in(im)
			TPZ%membrane_diff_in(im) = TPZ%membrane_diff_in(im)*Vsite_cm3/60		! /min -> /sec
			read(nfcell,*) TPZ%membrane_diff_out(im)
			TPZ%membrane_diff_out(im) = TPZ%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
			read(nfcell,*) TPZ%halflife(im)
			ichemo = idrug + im
			if (im > 0) then
				write(numstr,'(i1)') im
				chemo(ichemo)%name = 'TPZ_metab_'//numstr
			endif
			chemo(ichemo)%diff_coef = TPZ%diff_coef(im)
			chemo(ichemo)%medium_diff_coef = TPZ%medium_diff_coef(im)
			chemo(ichemo)%membrane_diff_in = TPZ%membrane_diff_in(im)
			chemo(ichemo)%membrane_diff_out = TPZ%membrane_diff_out(im)
			chemo(ichemo)%halflife = TPZ%halflife(im)
		enddo
		do ictype = 1,Ncelltypes
			do im = 0,nmetab
				read(nfcell,*) TPZ%Kmet0(ictype,im)
				read(nfcell,*) TPZ%C2(ictype,im)
				read(nfcell,*) TPZ%KO2(ictype,im)
				read(nfcell,*) TPZ%Vmax(ictype,im)
				read(nfcell,*) TPZ%Km(ictype,im)
				read(nfcell,*) TPZ%Klesion(ictype,im)
				if (im == 0) then
					read(nfcell,*) TPZ%kill_model(ictype)
					read(nfcell,*) TPZ%kill_O2(ictype)
					read(nfcell,*) TPZ%kill_drug(ictype)
					read(nfcell,*) TPZ%kill_duration(ictype)
					read(nfcell,*) TPZ%kill_fraction(ictype)
				endif
			enddo
			TPZ%Kmet0(ictype,:) = TPZ%Kmet0(ictype,:)/60				! /min -> /sec
			TPZ%KO2(ictype,:) = 1.0e-3*TPZ%KO2(ictype,:)                ! um -> mM
			TPZ%kill_duration(ictype) = 60*TPZ%kill_duration(ictype)    ! minutes -> seconds
		enddo

	elseif (idrug == DNB_DRUG) then
		DNB%name = drug_name
		DNB%nmetabolites = nmetab
		DNB%use_metabolites = use_metabolites
		do im = 0,nmetab
			read(nfcell,*) DNB%diff_coef(im)
			read(nfcell,*) DNB%medium_diff_coef(im)
			read(nfcell,*) DNB%membrane_diff_in(im)
			DNB%membrane_diff_in(im) = DNB%membrane_diff_in(im)*Vsite_cm3/60		! /min -> /sec
			read(nfcell,*) DNB%membrane_diff_out(im)
			DNB%membrane_diff_out(im) = DNB%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
			read(nfcell,*) DNB%halflife(im)
			ichemo = idrug + im
			if (im > 0) then
				write(numstr,'(i1)') im
				chemo(ichemo)%name = 'DNB_metab_'//numstr
			endif
			chemo(ichemo)%diff_coef = DNB%diff_coef(im)
			chemo(ichemo)%medium_diff_coef = DNB%medium_diff_coef(im)
			chemo(ichemo)%membrane_diff_in = DNB%membrane_diff_in(im)
			chemo(ichemo)%membrane_diff_out = DNB%membrane_diff_out(im)
			chemo(ichemo)%halflife = DNB%halflife(im)	
		enddo
		do ictype = 1,Ncelltypes
			do im = 0,nmetab
				read(nfcell,*) DNB%Kmet0(ictype,im)
				read(nfcell,*) DNB%C2(ictype,im)
				read(nfcell,*) DNB%KO2(ictype,im)
				read(nfcell,*) DNB%Vmax(ictype,im)
				read(nfcell,*) DNB%Km(ictype,im)
				read(nfcell,*) DNB%Klesion(ictype,im)
				read(nfcell,*) DNB%kill_model(ictype,im)
				read(nfcell,*) DNB%kill_O2(ictype,im)
				read(nfcell,*) DNB%kill_drug(ictype,im)
				read(nfcell,*) DNB%kill_duration(ictype,im)
				read(nfcell,*) DNB%kill_fraction(ictype,im)
			enddo
			DNB%Kmet0(ictype,:) = DNB%Kmet0(ictype,:)/60					! /min -> /sec
			DNB%KO2(ictype,:) = 1.0e-3*DNB%KO2(ictype,:)					! um -> mM
			DNB%kill_duration(ictype,:) = 60*DNB%kill_duration(ictype,:)    ! minutes -> seconds
		enddo
	endif
	do im = 0,nmetab
		ichemo = idrug + im
		chemo(ichemo)%decay = (chemo(ichemo)%halflife > 0)
		if (chemo(ichemo)%used .and. chemo(ichemo)%decay) then
			chemo(ichemo)%decay_rate = DecayRate(chemo(ichemo)%halflife)
		else
			chemo(ichemo)%decay_rate = 0
		endif
!		if (chemo(imetab)%used .and. chemo(imetab)%decay) then
!			chemo(imetab)%decay_rate = DecayRate(chemo(imetab)%halflife)
!		else
!			chemo(imetab)%decay_rate = 0
!		endif
	enddo
enddo

read(nfcell,*) LQ(1)%alpha_H
read(nfcell,*) LQ(1)%beta_H
read(nfcell,*) LQ(1)%OER_am
read(nfcell,*) LQ(1)%OER_bm
read(nfcell,*) LQ(1)%K_ms
read(nfcell,*) LQ(1)%death_prob
read(nfcell,*) LQ(2)%alpha_H
read(nfcell,*) LQ(2)%beta_H
read(nfcell,*) LQ(2)%OER_am
read(nfcell,*) LQ(2)%OER_bm
read(nfcell,*) LQ(2)%K_ms
read(nfcell,*) LQ(2)%death_prob
read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) growthcutoff(1)
read(nfcell,*) growthcutoff(2)
read(nfcell,*) growthcutoff(3)
read(nfcell,*) spcrad_value
read(nfcell,*) iuse_extra
read(nfcell,*) iuse_relax
read(nfcell,*) iuse_par_relax
read(nfcell,*) iuse_drop
read(nfcell,*) Ndrop
read(nfcell,*) alpha_shape
read(nfcell,*) beta_shape
read(nfcell,*) isaveprofiledata
read(nfcell,*) profiledatafilebase
read(nfcell,*) dt_saveprofiledata
read(nfcell,*) nt_saveprofiledata

if (use_events) then
	call ReadProtocol(nfcell)
	use_treatment = .false.
endif

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
O2cutoff = O2cutoff/1000							! uM -> mM
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

divide_dist(1:2)%class = LOGNORMAL_DIST
divide_time_median(1:2) = 60*60*divide_time_median(1:2)		! hours -> seconds
sigma(1:2) = log(divide_time_shape(1:2))
!divide_dist%p1 = log(divide_time_mean/exp(sigma*sigma/2))	
divide_dist(1:2)%p1 = log(divide_time_median(1:2))	
divide_dist(1:2)%p2 = sigma(1:2)
divide_time_mean(1:2) = exp(divide_dist(1:2)%p1 + 0.5*divide_dist(1:2)%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,24e12.4)') 'shape, sigma: ',divide_time_shape(1:2),sigma(1:2)
call logger(logmsg)
write(logmsg,'(a,4e12.4)') 'Median, mean divide time: ',divide_time_median(1:2)/3600,divide_time_mean(1:2)/3600
call logger(logmsg)

use_V_dependence = (iV_depend == 1)
saveprofiledata = (isaveprofiledata == 1)
dt_saveprofiledata = 60*dt_saveprofiledata			! mins -> seconds
use_dropper = (iuse_drop == 1)

randomise_initial_volume = (iV_random == 1)
use_extracellular_O2 = (iuse_extra == 1)

open(nfout,file=outputfile,status='replace')
write(nfout,'(a,a)') 'GUI version: ',gui_run_version
write(nfout,'(a,a)') 'DLL version: ',dll_run_version
write(nfout,*)

write(nflog,*)
write(nflog,'(a,a)') 'GUI version: ',gui_run_version
write(nflog,'(a,a)') 'DLL version: ',dll_run_version
write(nflog,*)

open(nfres,file='scell_ts.out',status='replace')
write(nfres,'(a,a)') 'GUI version: ',gui_run_version
write(nfres,'(a,a)') 'DLL version: ',dll_run_version
write(nfres,*)
write(nfres,'(a)') 'istep hour vol_mm3 diam_um Ncells(2) &
Nanoxia_dead(2) NdrugA_dead(2) NdrugB_dead(2) Nradiation_dead(2) &
Ntagged_anoxia(2) Ntagged_drugA(2) Ntagged_drugB(2) Ntagged_radiation(2) &
f_hypox_1 f_hypox_2 f_hypox_3 f_growth_1 f_growth_2 f_growth_3 f_necrot plating_efficiency(2) &
medium_oxygen medium_glucose medium_TPZ medium_DNB'

write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)
call DetermineKd
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getIndices(indx, idrug, nmetab)
integer :: indx, idrug, nmetab

if (indx == 1) then
	idrug = TPZ_DRUG
	nmetab = 2
elseif (indx == 2) then
	idrug = DNB_DRUG
	nmetab = 2
else
	idrug = 0
	nmetab = 0
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Skip lines until the 'PROTOCOL' line
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
integer :: itime, ntimes, kevent, ichemo, im
character*(64) :: line
real(REAL_KIND) :: t, dt, vol, conc, dose
type(event_type) :: E

chemo(TRACER+1:)%used = .false.
do
	read(nf,'(a)') line
	if (trim(line) == 'PROTOCOL') exit
enddo
read(nf,*) ntimes
if (allocated(event)) deallocate(event)
allocate(event(2*ntimes))
kevent = 0
do itime = 1,ntimes
	read(nf,'(a)') line
	if (trim(line) == 'DRUG') then
		kevent = kevent + 1
		event(kevent)%etype = DRUG_EVENT
		read(nf,'(a)') line
		if (trim(line) == 'SN30000' .or. trim(line) == 'TPZ') then
			ichemo = TPZ_DRUG
		elseif (trim(line) == 'PR104A') then
			ichemo = DNB_DRUG
		endif
		read(nf,*) t
		read(nf,*) dt
		read(nf,*) vol
		read(nf,*) conc
		event(kevent)%time = t
		event(kevent)%ichemo = ichemo
		event(kevent)%volume = vol
		event(kevent)%conc = conc
		event(kevent)%dose = 0
		chemo(ichemo)%used = .true.
		if (ichemo == TPZ_DRUG .and. TPZ%use_metabolites) then
			do im = 1,TPZ%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif
		if (ichemo == DNB_DRUG .and. DNB%use_metabolites) then
			do im = 1,DNB%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		event(kevent)%time = t + dt
		event(kevent)%ichemo = 0
		event(kevent)%volume = medium_volume0*1.0e3		! -> uL for consistency
		event(kevent)%conc = 0
		event(kevent)%dose = 0
	elseif (trim(line) == 'MEDIUM') then
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		read(nf,*) t
		read(nf,*) vol
		event(kevent)%time = t
		event(kevent)%volume = vol	
		event(kevent)%ichemo = 0
		event(kevent)%conc = 0
		event(kevent)%dose = 0
	elseif (trim(line) == 'RADIATION') then
		kevent = kevent + 1
		event(kevent)%etype = RADIATION_EVENT
		read(nf,*) t
		read(nf,*) dose
		event(kevent)%time = t
		event(kevent)%dose = dose	
		event(kevent)%ichemo = 0
		event(kevent)%volume = 0
		event(kevent)%conc = 0
	endif
enddo
Nevents = kevent
! Set events not done
! convert time from hours to seconds
! convert volume from uL to cm^3
do kevent = 1,Nevents
	event(kevent)%done = .false.
	event(kevent)%time = event(kevent)%time*60*60
	event(kevent)%volume = event(kevent)%volume*1.0e-3
	E = event(kevent)
!	write(*,'(a,i3,f8.0,2i3,3f8.4)') 'event: ',kevent,E%time,E%etype,E%ichemo,E%volume,E%conc,E%dose
enddo
! Check that events are sequential
do kevent = 1,Nevents-1
	if (event(kevent)%time >= event(kevent+1)%time) then
		write(logmsg,*) 'Error: non-sequential event: ',kevent,event(kevent)%time
		call logger(logmsg)
		stop
	endif
enddo
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
	test_growthrate = Vdivide0/(2*(divide_time_median(1) - mitosis_duration))	! um3/sec
else
	mitosis_duration = 0
	test_growthrate = Vdivide0/(2*divide_time_median(1))	! um3/sec
endif
write(logmsg,'(a,2e12.3)') 'Raverage,Vdivide0: ',Raverage,Vdivide0
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
total_dMdt = 0

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
! SN30000 CT model
! ----------------
! The rate of cell killing, which becomes a cell death probability rate, is inferred from
! the cell kill experiment.
! The basic assumption is that the rate of killing depends on the drug metabolism rate.
! There are five models:
! kill_model = 1:
!   killing rate = c = Kd.dM/dt
! kill_model = 2:
!   killing rate = c = Kd.Ci.dM/dt
! kill_model = 3:
!   killing rate = c = Kd.(dM/dt)^2
! kill_model = 4:
!   killing rate = c = Kd.Ci
! kill_model = 5:
!   killing rate = c = Kd.Ci^2
! where dM/dt = F(O2).kmet0.Ci
! In the kill experiment both O2 and Ci are held constant:
! O2 = CkillO2, Ci = Ckill
! In this case c is constant and the cell population N(t) is given by:
! N(t) = N(0).exp(-ct), i.e.
! c = -log(N(T)/N(0))/T where T is the duration of the experiment
! N(T)/N(0) = 1 - f, where f = kill fraction
! kill_model = 1:
!   c = Kd.F(CkillO2).kmet0.Ckill => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill)
! kill_model = 2:
!   c = Kd.F(CkillO2).kmet0.Ckill^2 => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill^2)
! kill_model = 3:
!   c = Kd.(F(CkillO2).kmet0.Ckill)^2 => Kd = -log(1-f)/(T.(F(CkillO2).kmet0.Ckill)^2)
! kill_model = 4:
!   c = Kd.Ckill => Kd = -log(1-f)/(T.Ckill)
! kill_model = 5:
!   c = Kd.Ckill^2 => Kd = -log(1-f)/(T.Ckill^2)
!-----------------------------------------------------------------------------------------
subroutine DetermineKd
real(REAL_KIND) :: C2, KO2, Kmet0, kmet 
real(REAL_KIND) :: f, T, Ckill, Ckill_O2, Kd
integer :: idrug, i, im, kill_model

do idrug = 1,2
	if (idrug == 1 .and. .not.chemo(TPZ_DRUG)%used) cycle
	if (idrug == 2 .and. .not.chemo(DNB_DRUG)%used) cycle
	do i = 1,Ncelltypes
		do im = 0,2
			if (idrug == 1) then		! TPZ
				if (im /= 0) cycle
				C2 = TPZ%C2(i,im)
				KO2 = TPZ%KO2(i,im)
				Kmet0 = TPZ%Kmet0(i,im)
				kill_model = TPZ%kill_model(i)
				Ckill_O2 = TPZ%kill_O2(i)
				f = TPZ%kill_fraction(i)
				T = TPZ%kill_duration(i)
				Ckill = TPZ%kill_drug(i)
			elseif (idrug == 2) then	! DNB
				C2 = DNB%C2(i,im)
				KO2 = DNB%KO2(i,im)
				Kmet0 = TPZ%Kmet0(i,im)
				kill_model = DNB%kill_model(i,im)
				Ckill_O2 = DNB%kill_O2(i,im)
				f = DNB%kill_fraction(i,im)
				T = DNB%kill_duration(i,im)
				Ckill = DNB%kill_drug(i,im)
			endif
			kmet = (1 - C2 + C2*KO2/(KO2 + Ckill_O2))*Kmet0
			if (kill_model == 1) then
				Kd = -log(1-f)/(T*kmet*Ckill)
			elseif (kill_model == 2) then
				Kd = -log(1-f)/(T*kmet*Ckill**2)
			elseif (kill_model == 3) then
				Kd = -log(1-f)/(T*(kmet*Ckill)**2)
			elseif (kill_model == 4) then
				Kd = -log(1-f)/(T*Ckill)
			elseif (kill_model == 5) then
				Kd = -log(1-f)/(T*Ckill**2)
			endif
			if (idrug == 1) then
				TPZ%Kd(i) = Kd
			elseif (idrug == 2) then
				DNB%Kd(i,im) = Kd
			endif
		enddo
	enddo
enddo

!if (chemo(TPZ_DRUG)%used) then
!	do i = 1,Ncelltypes
!		kmet = (1 - TPZ%C2(i,0) + TPZ%C2(i,0)*TPZ%KO2(i,0)/(TPZ%KO2(i,0) + TPZ%kill_O2(i)))*TPZ%Kmet0(i,0)
!		kill_model = TPZ%kill_model(i)
!		if (kill_model == 1) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*kmet*TPZ%kill_drug(i))
!		elseif (kill_model == 2) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*kmet*TPZ%kill_drug(i)**2)
!		elseif (kill_model == 3) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*(kmet*TPZ%kill_drug(i))**2)
!		elseif (kill_model == 4) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*TPZ%kill_drug(i))
!		elseif (kill_model == 5) then
!			Kd = -log(1-TPZ%kill_fraction(i))/(TPZ%kill_duration(i)*TPZ%kill_drug(i)**2)
!		endif
!		TPZ%Kd(i) = Kd
!	enddo
!endif
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
d = 2.2*Raverage
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
cp%ID = kcell
cp%state = ALIVE
cp%celltype = 1		! temporary!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!cp%growthrate = test_growthrate
!cp2%V_divide = get_divide_volume()
cp%d_divide = (3*cp%V_divide/PI)**(1./3.)
cp%mitosis = 0
cp%drugA_tag = .false.
cp%drugB_tag = .false.
cp%anoxia_tag = .false.
cp%t_hypoxic = 0
call get_random_vector3(r)	! set initial axis direction
cp%d = 0.1*small_d
c = cp%centre(:,1)
cp%centre(:,1) = c + (cp%d/2)*r
cp%centre(:,2) = c - (cp%d/2)*r
cp%nbrs = 0
cp%Cex = Caverage(1,1,1,:)	! initially the concentrations are the same everywhere
cp%Cin = cp%Cex
cp%CFSE = generate_CFSE(1.d0)
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
real(REAL_KIND) :: radiation_dose, dt
integer :: i, k, nit, nt_diff, it_diff, ncells0, nhypoxic(3)
integer :: nshow = 100
integer :: Nhop
integer :: nvars, ns
real(REAL_KIND) :: dxc, ex_conc(35*O2_BY_VOL+1)		! just for testing
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
t_simulation = (istep-1)*DELTA_T	! seconds
radiation_dose = 0
if (use_events) then
	call ProcessEvent(radiation_dose)
endif
if (radiation_dose > 0) then
	write(logmsg,'(a,f6.1)') 'Radiation dose: ',radiation_dose
	call logger(logmsg)
endif
call SetupChemomap
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
!	call grower(dt,changed,ok)
	call GrowCells(radiation_dose,dt,changed,ok)
	if (.not.ok) then
		call logger('grower error')
		res = 2
		return
	endif
	radiation_dose = 0
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

call make_grid_flux_weights

! Reaction-diffusion system
! Assuming DELTA_T = 600 ...
!if (ncells < 2000) then
!	nt_diff = 1
!elseif (ncells < 3000) then
!	nt_diff = 2
!elseif (ncells < 4000) then
!	nt_diff = 3
!elseif (ncells < 5000) then
!	nt_diff = 4
!elseif (ncells < 7000) then
!	nt_diff = 5
!elseif (ncells < 10000) then
!	nt_diff = 6
!else
!	nt_diff = 7
!endif
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
if (.not.use_TCP .and. (mod(istep,6) == 0)) then
	call get_concdata(nvars, ns, dxc, ex_conc)
!	write(*,'(a,3f8.4)') 'cell #1: ',cell_list(1)%Cex(1),cell_list(1)%Cin(1),cell_list(1)%Cex(1)-cell_list(1)%Cin(1)
endif
res = 0
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine ProcessEvent(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: kevent, ichemo, im, nmetab
real(REAL_KIND) :: V, C(MAX_CHEMO)
type(event_type) :: E

do kevent = 1,Nevents
	E = event(kevent)
	if (t_simulation >= E%time .and. .not.E%done) then
!	write(*,'(i3,2f8.0,i3,2f10.4)') E%etype,t_simulation,E%time,E%ichemo,E%volume,E%conc
		if (E%etype == RADIATION_EVENT) then
			radiation_dose = E%dose
!			write(*,*) 'radiation_dose: ',radiation_dose
		elseif (E%etype == MEDIUM_EVENT) then
			C = 0
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			V = E%volume
			call MediumChange(V,C)
		elseif (E%etype == DRUG_EVENT) then
			C = 0
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			ichemo = E%ichemo
			C(ichemo) = E%conc
			V = E%volume
			! set %present
			chemo(ichemo)%present = .true.
			chemo(ichemo)%bdry_conc = 0
			if (ichemo == TPZ_DRUG) then
				nmetab = TPZ%nmetabolites
			elseif (ichemo == DNB_DRUG) then
				nmetab = DNB%nmetabolites
			endif			
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					chemo(ichemo + im)%bdry_conc = 0
				endif
			enddo
			call MediumChange(V,C)
		endif
		event(kevent)%done = .true.
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! If the volume removed is Vr, the fraction of constituent mass that is retained
! in the medium is (Vm - Vr)/Vm.  The calculation does not apply to oxygen.
! Usually Vr = Ve.
! With grid, need to exclude gridcells that contain cells and those interior gridcells.
! The total mass in the external gridcells is determined, and a fraction of it may be retained
!-----------------------------------------------------------------------------------------
subroutine MediumChange(Ve,Ce)
real(REAL_KIND) :: Ve, Ce(:)
real(REAL_KIND) :: R, Vm_old, Vm_new, Vr, Vkeep, Vblob, fkeep
integer :: kcell, site(3), siteb(3), ixb, iyb, izb, izb_1, izb_2, ichemo
integer, allocatable :: zinrng(:,:,:)
real(REAL_KIND), allocatable :: exmass(:), exconc(:)

allocate(ngcells(NXB,NYB,NZB))
allocate(zinrng(2,NXB,NYB))
allocate(exmass(MAX_CHEMO))
allocate(exconc(MAX_CHEMO))
exmass = 0
ngcells = 0
! need to identify gridcells containing cells, to exclude from the calculation
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	site = cell_list(kcell)%site	! this gives the location in the fine grid
	call getSiteb(site,siteb)
	ngcells(siteb(1),siteb(2),siteb(3)) = ngcells(siteb(1),siteb(2),siteb(3)) + 1
enddo
do ixb = 1,NXB
	do iyb = 1,NYB
		! Need to include all internal empty gridcells (necrotic)
		! First find first and last non-empty gridcells
		izb_1 = 0
		izb_2 = 0
		do izb = 1,NZB
			if (ngcells(ixb,iyb,izb) > 0) then
				if (izb_1 == 0) izb_1 = izb
				izb_2 = izb
			endif
		enddo
		if (izb_1 /= 0) then
			Vblob = Vblob + izb_2 - izb_1 + 1
		endif
		do izb = 1,NZB
			if (izb_1 /= 0 .and. (izb >= izb_1 .and. izb <= izb_2)) cycle	! blob gridcells
			do ichemo = 1,MAX_CHEMO
				exmass(ichemo) = exmass(ichemo) + dxb3*chemo(ichemo)%Cave_b(ixb,iyb,izb)
			enddo
		enddo
		zinrng(:,ixb,iyb) = [izb_1,izb_2]
	enddo
enddo
Vblob = dxb3*Vblob

R = getRadius()
Vblob = (4./3.)*PI*R**3			! cm3
Vm_old = total_volume - Vblob	! this is the initial external (medium) volume
Vr = min(Vm_old,Ve)				! this is the amount of medium volume that is exchanged
Vkeep = Vm_old - Vr
fkeep = Vkeep/Vm_old			! fraction of initial medium volume that is kept. i.e. fraction of exmass(:)
Vm_new = Ve + Vkeep				! new medium volume
total_volume = Vm_new + Vblob	! new total volume
! Concentrations in the gridcells external to the blob (those with no cells) are set
! to the values of the mixture of old values and added medium.
do ichemo = 1,MAX_CHEMO
	exconc(ichemo) = (Ve*Ce(ichemo) + fkeep*exmass(ichemo))/Vm_new
enddo
do ixb = 1,NXB
	do iyb = 1,NYB
		if (zinrng(1,ixb,iyb) == 0) then
			do ichemo = 1,MAX_CHEMO
				chemo(ichemo)%Cave_b(ixb,iyb,1:NZB) = exconc(ichemo)
			enddo
		else
			izb_1 = zinrng(1,ixb,iyb)
			izb_2 = zinrng(2,ixb,iyb)
			do izb = 1,NZB
				if ((izb >= izb_1 .and. izb <= izb_2)) cycle	! blob gridcells
				do ichemo = 1,MAX_CHEMO
					chemo(ichemo)%Cave_b(ixb,iyb,izb) = exconc(ichemo)
				enddo
			enddo
		endif
	enddo
enddo
!chemo(OXYGEN+1:)%medium_M = ((Vm - Vr)/Vm)*chemo(OXYGEN+1:)%medium_M + Vr*Ce(OXYGEN+1:)
!chemo(OXYGEN+1:)%medium_Cext = chemo(OXYGEN+1:)%medium_M/(total_volume - Vblob)
!chemo(OXYGEN)%medium_Cext = chemo(OXYGEN)%bdry_conc

deallocate(ngcells)
deallocate(zinrng)
deallocate(exmass)
deallocate(exconc)
end subroutine

!-----------------------------------------------------------------------------------------
! Determine which gridcell in the big grid (siteb(:)) the fine gridcell site(:) falls in.
! Note that %site(:) and %centre(:,:) always refer to the fine grid, because it is
! assumed that cells never move outside the fine grid.
!-----------------------------------------------------------------------------------------
subroutine getSiteb(site,siteb)
integer :: site(3), siteb(3)
integer :: ixb1, iyb1, ixb, iyb, izb

ixb1 = (NXB+1)/2 - (NX-1)/(2*NRF)
ixb = ixb1 + (site(1)-1)/NRF
iyb1 = (NYB+1)/2 - (NY-1)/(2*NRF)
iyb = iyb1 + (site(2)-1)/NRF
izb = 1 + (site(3)-1)/NRF
siteb = [ixb, iyb, izb]
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