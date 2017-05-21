! Cancer cell state development

module cellstate
use global
use chemokine
use nbr

implicit none

integer :: kcell_dividing = 0

contains

!-----------------------------------------------------------------------------------------
! Need to initialize site and cell concentrations when a cell divides and when there is
! cell death.
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dt,changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok

ok = .true.
!if (dose > 0) then
!	call Irradiation(dose, ok)
!	if (.not.ok) return
!endif
changed = .false.
!call CellGrowth(dt,ok)
call grower(dt,changed,OK)
if (.not.ok) return
if (use_death) then
	call CellDeath(dt,changed, ok)
	if (.not.ok) return
endif
if (use_migration) then
	call CellMigration(ok)
	if (.not.ok) return
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The O2 concentration to use with cell kcell is either the intracellular concentration,
! or is use_extracellular_O2, the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getO2conc(kcell, C_O2)
integer :: kcell
real(REAL_KIND) :: C_O2
integer :: iv, site(3)
real(REAL_KIND) :: tnow

if (use_extracellular_O2 .and. istep > 1) then		! fix 30/04/2015
	C_O2 = cell_list(kcell)%Cex(OXYGEN)
else
	C_O2 = cell_list(kcell)%Cin(OXYGEN)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The glucose concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2 (!), the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getGlucoseconc(kcell, C_glucose)
integer :: kcell
real(REAL_KIND) :: C_glucose
integer :: iv, site(3)
real(REAL_KIND) :: tnow

if (use_extracellular_O2 .and. istep > 1) then
	C_glucose = cell_list(kcell)%Cex(GLUCOSE)
else
	C_glucose = cell_list(kcell)%Cin(GLUCOSE)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose.
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, ityp, idrug, im, ichemo, kpar=0
real(REAL_KIND) :: C_O2, SER, p_death, p_recovery, R, kill_prob, tnow
real(REAL_KIND) :: Cs							! concentration of radiosensitising drug
real(REAL_KIND) :: SER_max0, SER_Km, SER_KO2	! SER parameters of the drug
real(REAL_KIND) :: SERmax						! max sensitisation at the drug concentration

ok = .true.
call logger('Irradiation')
tnow = istep*DELTA_T	! seconds
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%radiation_tag) cycle	! we do not tag twice (yet)
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	! Compute sensitisation SER
	SER = 1.0
	do idrug = 1,Ndrugs_used
		ichemo = 4 + 3*(idrug-1)
		if (.not.chemo(ichemo)%present) cycle
		do im = 0,2
			ichemo = 4 + 3*(idrug-1) + im
			if (drug(idrug)%sensitises(ityp,im)) then
				Cs = cell_list(kcell)%Cin(ichemo)	! concentration of drug/metabolite in the cell
				SER_max0 = drug(idrug)%SER_max(ityp,im)
				SER_Km = drug(idrug)%SER_Km(ityp,im)
				SER_KO2 = drug(idrug)%SER_KO2(ityp,im)
				SERmax = (Cs*SER_max0 + SER_Km)/(Cs + SER_Km)
				SER = SER*(C_O2 + SER_KO2*SERmax)/(C_O2 + SER_KO2)
			endif
		enddo
	enddo
	call get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
	kill_prob = 1 - p_recovery
	R = par_uni(kpar)
	if (R < kill_prob) then
		cell_list(kcell)%radiation_tag = .true.
		Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
		cell_list(kcell)%p_rad_death = p_death
		if (LQ(ityp)%growth_delay_N > 0 .and. cell_list(kcell)%Iphase) then
			cell_list(kcell)%growth_delay = .true.
			cell_list(kcell)%dt_delay = LQ(ityp)%growth_delay_factor*dose
			cell_list(kcell)%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
		else
			cell_list(kcell)%growth_delay = .false.
		endif
	elseif (use_radiation_growth_delay_all .and. LQ(ityp)%growth_delay_N > 0) then
		cell_list(kcell)%growth_delay = .true.
		cell_list(kcell)%dt_delay = LQ(ityp)%growth_delay_factor*dose
		cell_list(kcell)%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
	else
		cell_list(kcell)%growth_delay = .false.
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A cell that receives a dose of radiation either recovers completely before reaching 
! mitosis or retains damage that has a probability of causing cell death during mitosis.
! A damaged cell that does not die at this point passes the damage on to the progeny cells.
! The probability of complete recovery = p_recovery = p_r
! The probability of death for a famaged cell at mitosis = p_death = p_d
! To ensure that the short-term death probability is consistent with the previous
! LQ formulation, we require p_d(1-p_r) = kill_prob as previously calculated.
! If p_d is determined (currently it is fixed), then 1-p_r = kill_prob/p_d,
! therefore p_r = 1 - kill_prob/p_d
!-----------------------------------------------------------------------------------------
subroutine get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_recovery, p_death
real(REAL_KIND) :: SER						! sensitisation
real(REAL_KIND) :: OER_alpha_d, OER_beta_d, expon, kill_prob_orig

OER_alpha_d = dose*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
OER_beta_d = dose*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
!expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_alpha_d**2		! 07/08/2015

OER_alpha_d = OER_alpha_d*SER
OER_beta_d = OER_beta_d*SER

expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_beta_d**2
p_recovery = exp(-expon)	! = SF
p_death = LQ(ityp)%death_prob

!kill_prob_orig = 1 - exp(-expon)
!call get_pdeath(ityp,dose,C_O2,p_death)
!p_recovery = 1 - kill_prob_orig/p_death
end subroutine

!-----------------------------------------------------------------------------------------
! This is the probability of death at time of division of cell that received a dose of 
! radiation and did not recover.
! In general one would expect this to depend of damage, i.e. on dose and C_O2, but
! for now it is a constant for a cell type
!-----------------------------------------------------------------------------------------
subroutine get_pdeath(ityp,dose,C_O2,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_death

p_death = LQ(ityp)%death_prob
end subroutine

!-----------------------------------------------------------------------------------------
! Cells move to preferable nearby sites.
! For now this is turned off - need to formulate a sensible interpretation of "preferable"
!-----------------------------------------------------------------------------------------
subroutine CellMigration(ok)
logical :: ok
integer :: kcell, j, indx, site0(3), site(3), jmax
real(REAL_KIND) :: C0(MAX_CHEMO), C(MAX_CHEMO), v0, v, vmax, d0, d

call logger('CellMigration is not yet implemented')
ok = .false.
return

end subroutine


!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia or glucosia, or they can be tagged 
! for death at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,changed,ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, killmodel, kpar=0 
real(REAL_KIND) :: C_O2, C_glucose, kmet, Kd, dMdt, Cdrug, n_O2, kill_prob, dkill_prob, death_prob, tnow
real(REAL_KIND) :: SF, survival_prob
logical :: anoxia_death, aglucosia_death
logical :: flag
real(REAL_KIND) :: p_tag = 1.0		! 0.3 ?? 
type(drug_type), pointer :: dp

!call logger('CellDeath')
ok = .true.
tnow = istep*DELTA_T	! seconds
anoxia_death = chemo(OXYGEN)%controls_death
aglucosia_death = chemo(GLUCOSE)%controls_death
nlist0 = nlist
nlt_threshold = 0
flag = .false.
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	if (cell_list(kcell)%anoxia_tag) then
		if (tnow >= cell_list(kcell)%t_anoxia_die) then	! .and. par_uni(kpar) < p_tag*dt/DELTA_T) then	???
			call CellDies(kcell)
			changed = .true.
			Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
!			write(logmsg,*) 'death from anoxia: Nanoxia_tag, Nanoxia_dead: ',Nanoxia_tag(ityp),Nanoxia_dead(ityp)
!			call logger(logmsg)
			cycle
		endif
	else
		if (anoxia_death .and. C_O2 < anoxia_threshold) then
			nlt_threshold = nlt_threshold + 1
			cell_list(kcell)%t_anoxia = cell_list(kcell)%t_anoxia + dt
			if (cell_list(kcell)%t_anoxia > t_anoxia_limit) then
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
				cell_list(kcell)%t_anoxia_tag = tnow
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
!				write(logmsg,'(a,2i6,2f8.1)') 'Tagged for anoxia: kcell, Nanoxia_tag: ',kcell,Nanoxia_tag(ityp),tnow,cell_list(kcell)%t_anoxia_die
!				call logger(logmsg)
			endif
		else
			cell_list(kcell)%t_anoxia = 0
		endif
	endif
	call getGlucoseconc(kcell,C_glucose)
	if (cell_list(kcell)%aglucosia_tag) then
		if (tnow >= cell_list(kcell)%t_aglucosia_die) then
			call CellDies(kcell)
			changed = .true.
			Naglucosia_dead(ityp) = Naglucosia_dead(ityp) + 1
			cycle
		endif
	else
		if (aglucosia_death .and. C_O2 < aglucosia_threshold) then
			cell_list(kcell)%t_aglucosia = cell_list(kcell)%t_aglucosia + dt
			if (cell_list(kcell)%t_aglucosia > t_aglucosia_limit) then
				cell_list(kcell)%aglucosia_tag = .true.						    ! tagged to die later
				cell_list(kcell)%t_aglucosia_die = tnow + aglucosia_death_delay	! time that the cell will die
				Naglucosia_tag(ityp) = Naglucosia_tag(ityp) + 1
			endif
		else
			cell_list(kcell)%t_aglucosia = 0
		endif
	endif
	do idrug = 1,ndrugs_used	
		ichemo = DRUG_A + 3*(idrug-1)	
		if (.not.chemo(ichemo)%present) cycle
		if (cell_list(kcell)%drug_tag(idrug)) cycle	! don't tag more than once
		dp => drug(idrug)
		kill_prob = 0
		death_prob = 0
		survival_prob = 1
		do im = 0,2
			if (.not.dp%kills(ityp,im)) cycle
			Cdrug = cell_list(kcell)%Cin(ichemo + im)
			if (Cdrug == 0) cycle
!			if (.not.flag) write(nflog,*) 'im: ',im
			killmodel = dp%kill_model(ityp,im)		! could use %drugclass to separate kill modes
			Kd = dp%Kd(ityp,im)
			n_O2 = dp%n_O2(ityp,im)
			kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)**n_O2/(dp%KO2(ityp,im)**n_O2 + C_O2**n_O2))*dp%Kmet0(ityp,im)
			if (.not.flag .and. C_O2 < 0.01) write(nflog,'(a,6e12.3)') 'C_O2,Kmet0,kmet,Kd,Cdrug,dt: ',C_O2,dp%Kmet0(ityp,im),kmet,Kd,Cdrug,dt
			dMdt = kmet*Cdrug
			call getDrugSF(killmodel,Kd,dMdt,Cdrug,dt,SF)
			survival_prob = SF*survival_prob
			death_prob = max(death_prob,dp%death_prob(ityp,im))
		enddo
		kill_prob = 1 - survival_prob
		if (kill_prob == 0) cycle
		if (.not.flag .and. im==0 .and. idrug==1) then
			total_dose_time = total_dose_time + dt
			write(nflog,'(a,2f7.4,2f8.1)') 'kill_prob, SF: ',kill_prob,1-kill_prob,dt,total_dose_time
		endif
	    if (.not.cell_list(kcell)%drug_tag(idrug) .and. par_uni(kpar) < kill_prob) then	! don't tag more than once
			cell_list(kcell)%p_drug_death(idrug) = death_prob
			cell_list(kcell)%drug_tag(idrug) = .true.
            Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
		endif
		flag = .true.
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDrugSF(kill_model,Kd,dMdt,Cdrug,dt,SF)
integer :: kill_model
real(REAL_KIND) :: Kd, dMdt, Cdrug, dt, SF
real(REAL_KIND) :: c

if (kill_model == 1) then
	c = Kd*dMdt
elseif (kill_model == 2) then
	c = Kd*dMdt*Cdrug
elseif (kill_model == 3) then
	c = Kd*dMdt**2
elseif (kill_model == 4) then
	c = Kd*Cdrug
elseif (kill_model == 5) then
	c = Kd*Cdrug**2
endif
SF = exp(-c*dt)
end subroutine

!-----------------------------------------------------------------------------------------
! If the dying cell site is less than a specified fraction f_migrate of the blob radius,
! the site migrates towards the blob centre.
! %indx -> 0
! If the site is on the boundary, it is removed from the boundary list, and %indx -> OUTSIDE_TAG
! The cell contents should be released into the site.
!-----------------------------------------------------------------------------------------
subroutine CellDies(kcell)
integer :: kcell
integer :: ityp, idrug

cell_list(kcell)%state = DEAD
ityp = cell_list(kcell)%celltype
Ncells = Ncells - 1
Ncells_type(ityp) = Ncells_type(ityp) - 1
if (cell_list(kcell)%anoxia_tag) then
	Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
endif
if (cell_list(kcell)%aglucosia_tag) then
	Naglucosia_tag(ityp) = Naglucosia_tag(ityp) - 1
endif
do idrug = 1,ndrugs_used
	if (cell_list(kcell)%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
	endif
enddo
if (cell_list(kcell)%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
endif
ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(logmsg,'(a,i6,i6)') 'CellDies: ngaps > max_ngaps: ',ngaps,max_ngaps
    call logger(logmsg)
    stop
endif
gaplist(ngaps) = kcell

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
integer :: k, kcell, ityp, idrug, kpar=0
type(cell_type), pointer :: cp
real(REAL_KIND) :: rr(3), c(3), rad, tnow, d_desired, tgrowth(MAX_CELLTYPES), c_rate(MAX_CELLTYPES), r_mean(MAX_CELLTYPES)
real(REAL_KIND) :: R
integer :: ndivide, divide_list(1000)
logical :: drugkilled
logical :: divide

ok = .true.
changed = .false.
! Cancel CONTINUOUS_MITOSIS option 13/02/2016
!if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then
!	tgrowth = divide_time_mean
!else
tgrowth = divide_time_mean - mitosis_duration
c_rate(1:2) = log(2.0)/tgrowth(1:2)		! Note: to randomise divide time need to use random number, not mean!
r_mean(1:2) = Vdivide0/(2*tgrowth(1:2))

ndivide = 0
tnow = istep*DELTA_T + t_fmover
do k = 1,ncells
	kcell = perm_index(k)
	kcell_now = kcell
	cp => cell_list(kcell)
	ityp = cell_list(kcell)%celltype
	divide = .false.
!	if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then	! always Mphase, growing and separating
!		! need to set initial mitosis axis and d_divide at birth
!		call growcell(cp,dt,c_rate(ityp),r_mean(ityp))
!		cp%mitosis = cp%V/cp%divide_volume
!		cp%d = cp%mitosis*cp%d_divide
!		rr = cp%centre(:,2) - cp%centre(:,1)
!		rr = rr/sqrt(dot_product(r,r))	! axis direction unit vector
!		c = (cp%centre(:,1) + cp%centre(:,2))/2
!		cp%centre(:,1) = c - (cp%d/2)*rr
!		cp%centre(:,2) = c + (cp%d/2)*rr
!		call cubic_solver(cp%d,cp%V,rad)
!		cp%radius = rad
!		if (cp%d > 2*rad) then			! time for completion of cell division
!			divide = .true.
!		endif
!		write(nflog,'(a,i4,3f8.3)') 'V,d,rad: ',kcell,cp%V,cp%d,rad
!	elseif (MITOSIS_MODE == TERMINAL_MITOSIS) then
	if (cp%Iphase) then
!		call growcell(cp,dt,c_rate(ityp),r_mean(ityp))
		call growcell(cp,dt)
		cp%radius(1) = (3*cp%V/(4*PI))**(1./3.)
!		write(*,*) 'grower: V: ',cp%V
		if (cp%V > cp%divide_volume) then	! time to divide
!			if (cp%radiation_tag .and..not.cp%G2_M) then
!				if (par_uni(kpar) < cp%p_rad_death) then
!					call CellDies(kcell)
!					changed = .true.
!					Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
!					cycle
!				endif
!			endif
			drugkilled = .false.
			do idrug = 1,ndrugs_used
				if (cp%drug_tag(idrug)) then
					call CellDies(kcell)
					changed = .true.
					Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
					drugkilled = .true.
					exit
				endif
			enddo
			if (drugkilled) cycle
			
			if (cell_list(kcell)%growth_delay) then
				if (cell_list(kcell)%G2_M) then
					if (tnow > cell_list(kcell)%t_growth_delay_end) then
						cell_list(kcell)%G2_M = .false.
					else
						cycle
					endif
				else
					cell_list(kcell)%t_growth_delay_end = tnow + cell_list(kcell)%dt_delay
					cell_list(kcell)%G2_M = .true.
					cycle
				endif
			endif
			! try moving death prob test to here
			if (cell_list(kcell)%radiation_tag) then
				R = par_uni(kpar)
				if (R < cp%p_rad_death) then
					call CellDies(kcell)
					changed = .true.
					Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
					cycle
				endif
			endif
		
!			if (cp%radiation_tag) then
!				if (cp%growth_delay) then
!					if (.not.cp%G2_M) then	! this is the first time to divide
!						cp%t_growth_delay_end = tnow + cp%dt_delay
!						cp%G2_M = .true.
!					endif
!					if (tnow > cp%t_growth_delay_end) then
!						cp%G2_M = .false.	! delay ended, time to enter mitosis
!					else
!						cycle				! continue delaying mitosis
!					endif
!				endif
!			endif
		
			cp%Iphase = .false.
            cp%nspheres = 2
			ncells_mphase = ncells_mphase + 1
			call get_random_vector3(rr)	! set initial axis direction
			cp%d = 0.1*small_d
			c = cp%centre(:,1)
			cp%centre(:,1) = c + (cp%d/2)*rr
			cp%centre(:,2) = c - (cp%d/2)*rr
			cp%mitosis = 0
			cp%t_start_mitosis = tnow
			cp%d_divide = 2.0**(2./3)*cp%radius(1)
		endif
	else
		cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
		d_desired = max(cp%mitosis*cp%d_divide,small_d)
		rr = cp%centre(:,2) - cp%centre(:,1)
		cp%d = sqrt(dot_product(rr,rr))
		rr = rr/cp%d	! axis direction
		c = (cp%centre(:,1) + cp%centre(:,2))/2
		cp%site = c/DELTA_X + 1
!			cp%centre(1,:) = c - (cp%d/2)*r		! For fmover we do not need to set centre positions
!			cp%centre(2,:) = c + (cp%d/2)*r
		call cubic_solver(d_desired,cp%V,rad)
		cp%radius = rad
		if (cp%d > 2*rad) then			! time for completion of cell division
			divide = .true.
		endif
	endif
	if (divide) then
!		if (cp%radiation_tag) then
!			if (par_uni(kpar) < cp%p_rad_death) then
!				call CellDies(kcell)
!				changed = .true.
!				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
!			endif
!			cycle
!		endif
!		drugkilled = .false.
!		do idrug = 1,ndrugs_used
!			if (cp%drug_tag(idrug)) then
!				call CellDies(kcell)
!				changed = .true.
!				Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
!				drugkilled = .true.
!				exit
!			endif
!		enddo
!		if (drugkilled) cycle
		ndivide = ndivide + 1
		divide_list(ndivide) = kcell
	endif
enddo
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	call divider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine old_growcell(cp, dt, c_rate, r_mean)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt, c_rate, r_mean
real(REAL_KIND) :: Cin_0(NCONST), Cex_0(NCONST)		! at some point NCONST -> MAX_CHEMO
real(REAL_KIND) :: dVdt,  Vin_0, dV, metab_O2, metab_glucose, metab, dVdt_new
logical :: glucose_growth
integer :: C_option = 1	! we must use this

glucose_growth = chemo(GLUCOSE)%controls_growth
Cin_0 = cp%Cin
metab_O2 = O2_metab(Cin_0(OXYGEN))	! Note that currently growth depends only on O2
metab_glucose = glucose_metab(Cin_0(GLUCOSE))
if (glucose_growth) then
	metab = metab_O2*metab_glucose
else
	metab = metab_O2
endif
!if (use_V_dependence) then
!	dVdt = c_rate*metab*cp%V/(Vdivide0/2)
!else
!	dVdt = r_mean*metab
!!	write(*,'(a,2e12.3)') 'Vdivide0,tgrowth: ',Vdivide0,divide_time_mean(1) - mitosis_duration
!!	write(*,'(a,3e12.3)') 'r_mean, metab, dVdt: ',r_mean, metab, dVdt
!endif
dVdt = get_dVdt(cp,metab)
if (suppress_growth) then	! for checking solvers
	dVdt = 0
endif
Cex_0 = cp%Cex
cp%dVdt = dVdt
!cp%growthrate = dVdt
Vin_0 = cp%V
dV = dVdt*dt
cp%V = Vin_0 + dV
if (C_option == 1) then
	! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
	cp%Cin = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
elseif (C_option == 2) then
	! Calculation based on change in volumes without mass transfer of constituents
	cp%Cin = Vin_0*Cin_0/(Vin_0 + dV)
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine growcell(cp, dt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt
real(REAL_KIND) :: Cin_0(NCONST), Cex_0(NCONST)		! at some point NCONST -> MAX_CHEMO
real(REAL_KIND) :: dVdt,  Vin_0, dV, metab_O2, metab_glucose, metab, dVdt_new
logical :: oxygen_growth, glucose_growth, tagged
integer :: C_option = 1	! we must use this

tagged = cp%anoxia_tag .or. cp%aglucosia_tag
if (tagged) then
	cp%dVdt = 0
	return		! cells tagged to die of anoxia or aglucosia do not metabolise
endif
oxygen_growth = chemo(OXYGEN)%controls_growth
glucose_growth = chemo(GLUCOSE)%controls_growth
Cin_0 = cp%Cin
metab_O2 = O2_metab(Cin_0(OXYGEN))	! Note that currently growth depends only on O2 
metab_glucose = glucose_metab(Cin_0(GLUCOSE))
if (oxygen_growth .and. glucose_growth) then
	metab = metab_O2*metab_glucose
elseif (oxygen_growth) then
	metab = metab_O2
elseif (glucose_growth) then
	metab = metab_glucose
endif
dVdt = get_dVdt(cp,metab)
if (suppress_growth) then	! for checking solvers
	dVdt = 0
endif
Cex_0 = cp%Cex
cp%dVdt = dVdt
Vin_0 = cp%V
dV = dVdt*dt
cp%V = Vin_0 + dV
if (C_option == 1) then
	! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
	cp%Cin = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
elseif (C_option == 2) then
	! Calculation based on change in volumes without mass transfer of constituents
	cp%Cin = Vin_0*Cin_0/(Vin_0 + dV)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Need to account for time spent in mitosis.  Because there is no growth during mitosis,
! the effective divide_time must be reduced by mitosis_duration.
! Note that TERMINAL_MITOSIS is the only option.
!-----------------------------------------------------------------------------------------
function get_dVdt(cp, metab) result(dVdt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: metab, dVdt
integer :: ityp
real(REAL_KIND) :: r_mean, c_rate

if (use_V_dependence) then
	if (use_constant_divide_volume) then
		dVdt = metab*log(2.0)*cp%V/(cp%divide_time - mitosis_duration)
	else
		ityp = cp%celltype
		c_rate = log(2.0)/(divide_time_mean(ityp) - mitosis_duration)
		dVdt = c_rate*cp%V*metab
	endif
else
	if (use_constant_divide_volume) then
		dVdt = 0.5*metab*Vdivide0/(cp%divide_time  - mitosis_duration)
	else
		ityp = cp%celltype
		r_mean = 0.5*Vdivide0/(divide_time_mean(ityp) - mitosis_duration)
		dVdt = r_mean*metab
	endif
endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetInitialGrowthRate
integer :: kcell, ityp
real(REAL_KIND) :: C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt
logical :: oxygen_growth, glucose_growth
type(cell_type), pointer :: cp

oxygen_growth = chemo(OXYGEN)%controls_growth
glucose_growth = chemo(GLUCOSE)%controls_growth
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	C_O2 = chemo(OXYGEN)%bdry_conc
	C_glucose = cp%Cin(GLUCOSE)
	if (oxygen_growth .and. glucose_growth) then
	    metab_O2 = O2_metab(C_O2)
		metab_glucose = glucose_metab(C_glucose)
		metab = metab_O2*metab_glucose
	elseif (oxygen_growth) then
	    metab_O2 = O2_metab(C_O2)
		metab = metab_O2
	elseif (glucose_growth) then
		metab_glucose = glucose_metab(C_glucose)
		metab = metab_glucose
	endif
	dVdt = get_dVdt(cp,metab)
	if (suppress_growth) then	! for checking solvers
		dVdt = 0
	endif
	cp%dVdt = dVdt
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A single cell is replaced by two.
!-----------------------------------------------------------------------------------------
subroutine divider(kcell1, ok)
integer :: kcell1
logical :: ok
integer :: kcell2, ityp, nbrs0
real(REAL_KIND) :: r(3), c(3), cfse0, cfse1, V0, Tdiv
type(cell_type), pointer :: cp1, cp2

!write(*,*) 'divider:'
!write(logmsg,*) 'divider: ',kcell1 
!call logger(logmsg)
ok = .true.
tnow = istep*DELTA_T
cp1 => cell_list(kcell1)
if (ngaps > 0) then
    kcell2 = gaplist(ngaps)
    ngaps = ngaps - 1
else
	nlist = nlist + 1
	if (nlist > MAX_NLIST) then
		ok = .false.
		call logger('Error: Maximum number of cells MAX_NLIST has been exceeded.  Increase and rebuild.')
		return
	endif
	kcell2 = nlist
endif
ncells = ncells + 1
ityp = cp1%celltype
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
cp2 => cell_list(kcell2)
cp1%state = ALIVE
cp1%generation = cp1%generation + 1
V0 = cp1%V/2
cp1%V = V0
cp1%site = cp1%centre(:,1)/DELTA_X + 1
cp1%d = 0
cp1%birthtime = tnow
!cp1%divide_volume = get_divide_volume1()
cp1%divide_volume = get_divide_volume(ityp,V0,Tdiv)
cp1%divide_time = Tdiv
cp1%d_divide = (3*cp1%divide_volume/PI)**(1./3.)
cp1%mitosis = 0
cp1%t_divide_last = tnow
cfse0 = cp1%CFSE
cp1%CFSE = generate_CFSE(cfse0/2)
cfse1 = cfse0 - cp1%CFSE

cp1%drug_tag = .false.
cp1%anoxia_tag = .false.
cp1%t_anoxia = 0
cp1%aglucosia_tag = .false.
cp1%t_aglucosia = 0

if (cp1%growth_delay) then
	cp1%N_delayed_cycles_left = cp1%N_delayed_cycles_left - 1
	cp1%growth_delay = (cp1%N_delayed_cycles_left > 0)
!	write(*,*) 'growth_delay cell divides: ',kcell1,kcell2,cp1%N_delayed_cycles_left
endif
cp1%G2_M = .false.

nbrs0 = cp1%nbrs
cp1%nbrs = nbrs0 + 1
cp1%nbrlist(cp1%nbrs)%indx = kcell2
cp1%nbrlist(cp1%nbrs)%contact = .false.
cp1%nbrlist(cp1%nbrs)%contact(1,1) = .true.

cp2%ID = cp1%ID
cp2%celltype = cp1%celltype
cp2%state = ALIVE
cp2%generation = cp1%generation
cp2%V = V0
cp2%radius(1) = cp1%radius(2)
cp2%centre(:,1) = cp1%centre(:,2)
cp2%site = cp2%centre(:,1)/DELTA_X + 1
cp2%d = 0
cp2%birthtime = tnow
!cp2%divide_volume = get_divide_volume1()
cp2%divide_volume = get_divide_volume(ityp,V0,Tdiv)
cp2%divide_time = Tdiv
cp2%d_divide = (3*cp2%divide_volume/PI)**(1./3.)
cp2%mitosis = 0
cp2%t_divide_last = tnow
cp2%CFSE = cfse1

cp2%ID = cp1%ID
cp2%p_rad_death = cp1%p_rad_death
cp2%radiation_tag = cp1%radiation_tag
if (cp2%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
cp2%drug_tag = .false.
cp2%anoxia_tag = .false.
cp2%t_anoxia = 0
cp2%aglucosia_tag = .false.
cp2%t_aglucosia = 0

cp2%growth_delay = cp1%growth_delay
if (cp2%growth_delay) then
	cp2%dt_delay = cp1%dt_delay
	cp2%N_delayed_cycles_left = cp1%N_delayed_cycles_left
endif
cp2%G2_M = .false.

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
! Cancel CONTINUOUS_MITOSIS option 13/02/2016
!if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then
!	cp1%Iphase = .false.
!	cp2%Iphase = .false.
!	call get_random_vector3(r)	! set initial axis direction
!	cp1%d = 0.1*small_d
!	c = cp1%centre(:,1)
!	cp1%centre(:,1) = c + (cp1%d/2)*r
!	cp1%centre(:,2) = c - (cp1%d/2)*r
!	cp2%d = 0.1*small_d
!	c = cp2%centre(:,1)
!	cp2%centre(:,1) = c + (cp2%d/2)*r
!	cp2%centre(:,2) = c - (cp2%d/2)*r
!else
cp1%Iphase = .true.
cp1%nspheres = 1
cp2%Iphase = .true.
cp2%nspheres = 1
call update_nbrlist(kcell1)
call update_nbrlist(kcell2)
! Note: any cell that has kcell1 in its nbrlist now needs to add kcell2 to the nbrlist.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function get_divide_volume1() result(vol)
real(REAL_KIND) :: vol
integer :: kpar = 0
real(REAL_KIND) :: U

U = par_uni(kpar)
vol = Vdivide0 + dVdivide*(2*U-1)
end function

!-----------------------------------------------------------------------------------------
! Generate a random value for CFSE from a distribution with mean = average
! In the simplest case we can allow a uniform distribution about the average.
! Multiplying factor in the range (1-a, 1+a)
! Better to make it a Gaussian distribution: 
!  = average*(1+s*R)
! where R = N(0,1), s = std deviation
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function generate_CFSE(average)
real(REAL_KIND) :: average, std
integer :: kpar = 0
real(REAL_KIND) :: R

! Uniform distribution
!R = par_uni(kpar)
!generate_CFSE = (1 - a + 2*a*R)*average
! Gaussian distribution
R = par_rnor(kpar)	! N(0,1)
generate_CFSE = (1 + CFSE_std*R)*average
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

end module
