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
subroutine GrowCells(dose,dt,changed, ok)
real(REAL_KIND) :: dose, dt
logical :: changed, ok

!call logger('GrowCells')
ok = .true.
if (dose > 0) then
	call Irradiation(dose, ok)
	if (.not.ok) return
endif
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
! Irradiate cells with dose.
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, ityp, kpar=0
real(REAL_KIND) :: C_O2, p_death, p_recovery, R, kill_prob
!real(REAL_KIND) :: OER_alpha_d, OER_beta_d, expon

ok = .true.
call logger('Irradiation')
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%radiation_tag) cycle	! we do not tag twice (yet)
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
!	OER_alpha_d = dose*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
!	OER_beta_d = dose*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
!	expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_alpha_d**2
!	kill_prob = 1 - exp(-expon)
	call get_kill_probs(ityp,dose,C_O2,p_recovery,p_death)
	kill_prob = 1 - p_recovery
	R = par_uni(kpar)
	if (R < kill_prob) then
		cell_list(kcell)%radiation_tag = .true.
		Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
		cell_list(kcell)%p_death = p_death
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
subroutine get_kill_probs(ityp,dose,C_O2,p_recovery,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_recovery, p_death
real(REAL_KIND) :: OER_alpha_d, OER_beta_d, expon, kill_prob_orig

real(REAL_KIND) :: Cs						! concentration of radiosensitising drug
real(REAL_KIND) :: SERmax0, KmSER, KO2SER	! SER parameters of the drug
real(REAL_KIND) :: SERmax					! max sensitisation at the drug concentration
real(REAL_KIND) :: SER						! sensitisation

OER_alpha_d = dose*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
OER_beta_d = dose*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
!expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_alpha_d**2		! 07/08/2015

!SERmax = (Cs*SERmax0 + KmSER)/(Cs + KmSER)
!SER = (C_O2 + KO2SER*SERmax)/(C_O2 + KO2SER)
!OER_alpha_d = OER_alpha_d*SER
!OER_beta_d = OER_beta_d*SER

expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_beta_d**2
kill_prob_orig = 1 - exp(-expon)
call get_pdeath(ityp,dose,C_O2,p_death)
p_recovery = 1 - kill_prob_orig/p_death
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

#if 0
!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia, or they can be tagged for death 
! at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine CellDeath1(dt,changed,ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: kcell, ict, nlist0, site(3), i, im, ityp, kpar=0 
real(REAL_KIND) :: C_O2, kmet, Kd, dMdt, kill_prob, tnow
logical :: use_TPZ_DRUG, use_DNB_DRUG

!call logger('CellDeath')
ok = .true.
use_TPZ_DRUG = chemo(TPZ_DRUG)%used
use_DNB_DRUG = chemo(DNB_DRUG)%used
tnow = istep*DELTA_T	! seconds
nlist0 = nlist
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	if (cell_list(kcell)%anoxia_tag) then
!		write(logmsg,*) 'anoxia_tag: ',kcell,cell_list(kcell)%state,tnow,cell_list(kcell)%t_anoxia_die
!		call logger(logmsg)
		if (tnow >= cell_list(kcell)%t_anoxia_die) then
!			call logger('cell dies')
			call CellDies(kcell)
			changed = .true.
			Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
			if (cell_list(kcell)%drugA_tag) then
				NdrugA_tag(ityp) = NdrugA_tag(ityp) - 1
			endif
			if (cell_list(kcell)%drugB_tag) then
				NdrugB_tag(ityp) = NdrugB_tag(ityp) - 1
			endif
			if (cell_list(kcell)%radiation_tag) then
				Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
			endif
			cycle
		endif
	else
		if (C_O2 < ANOXIA_THRESHOLD) then
			cell_list(kcell)%t_hypoxic = cell_list(kcell)%t_hypoxic + dt
			if (cell_list(kcell)%t_hypoxic > t_anoxic_limit) then
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
			endif
		else
			cell_list(kcell)%t_hypoxic = 0
		endif
	endif
	if (use_TPZ_DRUG .and. .not.cell_list(kcell)%drugA_tag) then	
		ict = cell_list(kcell)%celltype
		Kd = TPZ%Kd(ict)
	    kmet = (1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + C_O2))*TPZ%Kmet0(ict,0)
	    dMdt = kmet*cell_list(kcell)%Cin(TPZ_DRUG)
	    if (TPZ%kill_model(ict) == 1) then
		    kill_prob = Kd*dMdt*dt
	    elseif (TPZ%kill_model(ict) == 2) then
		    kill_prob = Kd*dMdt*cell_list(kcell)%Cin(TPZ_DRUG)*dt
	    elseif (TPZ%kill_model(ict) == 3) then
		    kill_prob = Kd*dMdt**2*dt
	    elseif (TPZ%kill_model(ict) == 4) then
		    kill_prob = Kd*cell_list(kcell)%Cin(TPZ_DRUG)*dt
	    elseif (TPZ%kill_model(ict) == 5) then
		    kill_prob = Kd*cell_list(kcell)%Cin(TPZ_DRUG)**2*dt
		endif
	    if (par_uni(kpar) < kill_prob) then
            cell_list(kcell)%drugA_tag = .true.
            NdrugA_tag(ityp) = NdrugA_tag(ityp) + 1
!            write(nflog,'(a,2i6)') 'TPZ tagged: ',kcell,ict
		endif
	endif
	if (use_DNB_DRUG .and. .not.cell_list(kcell)%drugB_tag) then
		ict = cell_list(kcell)%celltype
!	    kmet = (1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + C_O2))*TPZ%Kmet0(ict,0)
!	    dMdt = kmet*cell_list(kcell)%Cin(TPZ_DRUG)
	    if (DNB%kill_model(ict,1) < 4 .or. DNB%kill_model(ict,2) < 4) then
			write(logmsg,*) 'Error: CellDeath: DNB kill model must be 4 or 5, not: ',DNB%kill_model(ict,1),DNB%kill_model(ict,2)
			call logger(logmsg)
			ok = .false.
			return
		endif
		kill_prob = 0
		do im = 1,2
			if (DNB%kill_model(ict,im) == 4) then
				kill_prob = kill_prob + DNB%Kd(ict,im)*cell_list(kcell)%Cin(DNB_DRUG + im)*dt
			elseif (DNB%kill_model(ict,im) == 5) then
				kill_prob = kill_prob + DNB%Kd(ict,im)*(cell_list(kcell)%Cin(DNB_DRUG + im)**2)*dt
			endif
		enddo
	    if (par_uni(kpar) < kill_prob) then
            cell_list(kcell)%drugB_tag = .true.
            NdrugB_tag(ityp) = NdrugB_tag(ityp) + 1
            write(nflog,'(a,2i6)') 'DNB tagged: ',kcell,ict
		endif
	endif
enddo
end subroutine
#endif

!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia, or they can be tagged for death 
! at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,changed,ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, kpar=0 
real(REAL_KIND) :: C_O2, kmet, Kd, dMdt, killmodel, kill_prob, tnow
!logical :: use_TPZ_DRUG, use_DNB_DRUG
type(drug_type), pointer :: dp

!call logger('CellDeath')
ok = .true.
!use_TPZ_DRUG = chemo(TPZ_DRUG)%used
!use_DNB_DRUG = chemo(DNB_DRUG)%used
tnow = istep*DELTA_T	! seconds
nlist0 = nlist
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	if (cell_list(kcell)%anoxia_tag) then
!		write(logmsg,*) 'anoxia_tag: ',kcell,cell_list(kcell)%state,tnow,cell_list(kcell)%t_anoxia_die
!		call logger(logmsg)
		if (tnow >= cell_list(kcell)%t_anoxia_die) then
!			call logger('cell dies')
			call CellDies(kcell)
			changed = .true.
			Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
!			if (cell_list(kcell)%drugA_tag) then
!				NdrugA_tag(ityp) = NdrugA_tag(ityp) - 1
!			endif
!			if (cell_list(kcell)%drugB_tag) then
!				NdrugB_tag(ityp) = NdrugB_tag(ityp) - 1
!			endif
			do idrug = 1,ndrugs_used
				if (cell_list(kcell)%drug_tag(idrug)) then
					Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
				endif
			enddo
			if (cell_list(kcell)%radiation_tag) then
				Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
			endif
			cycle
		endif
	else
		if (C_O2 < ANOXIA_THRESHOLD) then
			cell_list(kcell)%t_hypoxic = cell_list(kcell)%t_hypoxic + dt
			if (cell_list(kcell)%t_hypoxic > t_anoxic_limit) then
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
			endif
		else
			cell_list(kcell)%t_hypoxic = 0
		endif
	endif
	
	do idrug = 1,ndrugs_used	
		dp => drug(idrug)
		ichemo = TRACER + 1 + 3*(idrug-1)	
		kill_prob = 0
		do im = 0,2
			if (.not.dp%kills(ityp,im)) cycle
			killmodel = dp%kill_model(ityp,im)		! could use %drugclass to separate kill modes
			Kd = dp%Kd(ityp,im)
			kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)/(dp%KO2(ityp,im) + C_O2))*dp%Kmet0(ityp,im)
!			dMdt = kmet*cell_list(kcell)%conc(ichemo + im)
			dMdt = kmet*cell_list(kcell)%Cin(ichemo + im)
			if (killmodel == 1) then
				kill_prob = kill_prob + Kd*dMdt*dt
			elseif (killmodel == 2) then
!				kill_prob = kill_prob + Kd*dMdt*cell_list(kcell)%conc(ichemo + im)*dt
				kill_prob = kill_prob + Kd*dMdt*cell_list(kcell)%Cin(ichemo + im)*dt
			elseif (killmodel == 3) then
				kill_prob = kill_prob + Kd*dMdt**2*dt
			elseif (killmodel == 4) then
!				kill_prob = kill_prob + Kd*cell_list(kcell)%conc(ichemo + im)*dt
				kill_prob = kill_prob + Kd*cell_list(kcell)%Cin(ichemo + im)*dt
			elseif (killmodel == 5) then
!				kill_prob = kill_prob + Kd*(cell_list(kcell)%conc(ichemo + im)**2)*dt
				kill_prob = kill_prob + Kd*(cell_list(kcell)%Cin(ichemo + im)**2)*dt
			endif
		enddo
	    if (par_uni(kpar) < kill_prob) then
!            cell_list(kcell)%drugB_tag = .true.			! actually either drugA_tag or drugB_tag
!            NdrugB_tag(ityp) = NdrugB_tag(ityp) + 1
			cell_list(kcell)%drug_tag(idrug) = .true.
            Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
		endif
	enddo
enddo
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
integer :: ityp

cell_list(kcell)%state = DEAD
ityp = cell_list(kcell)%celltype
Ncells = Ncells - 1
Ncells_type(ityp) = Ncells_type(ityp) - 1

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
real(REAL_KIND) :: r(3), c(3), rad, tnow, d_desired, tgrowth(MAX_CELLTYPES), c_rate(MAX_CELLTYPES), r_mean(MAX_CELLTYPES)
integer :: ndivide, divide_list(1000)
logical :: drugkilled
logical :: divide

ok = .true.
changed = .false.
if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then
	tgrowth = divide_time_mean
else
	tgrowth = divide_time_mean - mitosis_duration
endif
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
!	write(*,*) 'test_growthrate, dVdt: ',test_growthrate*dt, cp%dVdt
!	if (cp%state == DEAD) cycle
	if (MITOSIS_MODE == CONTINUOUS_MITOSIS) then	! always Mphase, growing and separating
		! need to set initial mitosis axis and d_divide at birth
!		cp%V = cp%V + cp%growthrate*DELTA_T
		call growcell(cp,dt,c_rate(ityp),r_mean(ityp))
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
			divide = .true.
!			ndivide = ndivide + 1
!			divide_list(ndivide) = kcell
		endif
		write(nflog,'(a,i4,3f8.3)') 'V,d,rad: ',kcell,cp%V,cp%d,rad
	elseif (MITOSIS_MODE == TERMINAL_MITOSIS) then
		if (cp%Iphase) then
!			cp%V = cp%V + cp%growthrate*dt
			call growcell(cp,dt,c_rate(ityp),r_mean(ityp))
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
			cp%radius = rad
			if (cp%d > 2*rad) then			! time for completion of cell division
				divide = .true.
!				ndivide = ndivide + 1
!				divide_list(ndivide) = kcell
			endif
		endif
	endif
	if (divide) then
		if (cp%radiation_tag) then
			if (par_uni(kpar) < cp%p_death) then
				call CellDies(kcell)
				changed = .true.
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
			endif
!			if (cp%drugA_tag) then
!				NdrugA_tag(ityp) = NdrugA_tag(ityp) - 1
!			endif
!			if (cp%drugB_tag) then
!				NdrugB_tag(ityp) = NdrugB_tag(ityp) - 1
!			endif
			do idrug = 1,ndrugs_used
				if (cp%drug_tag(idrug)) then
					Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
				endif
			enddo
			if (cp%anoxia_tag) then
				Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
			endif
			cycle
		endif
!		if (cp%drugA_tag) then
!			call CellDies(kcell)
!			NdrugA_dead(ityp) = NdrugA_dead(ityp) + 1
!			if (cp%anoxia_tag) then
!				Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
!			endif
!			if (cp%radiation_tag) then
!				Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
!			endif
!			cycle
!		endif
!		if (cp%drugB_tag) then
!			call CellDies(kcell)
!			NdrugB_dead(ityp) = NdrugB_dead(ityp) + 1
!			if (cp%anoxia_tag) then
!				Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
!			endif
!			if (cp%radiation_tag) then
!				Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
!			endif
!			cycle
!		endif
		drugkilled = .false.
		do idrug = 1,ndrugs_used
			if (cp%drug_tag(idrug)) then
				call CellDies(kcell)
				changed = .true.
				Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
				if (cp%anoxia_tag) then
					Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
				endif
				if (cp%radiation_tag) then
					Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
				endif
				drugkilled = .true.
				exit
			endif
		enddo
		if (drugkilled) cycle
		ndivide = ndivide + 1
		divide_list(ndivide) = kcell
	endif
enddo
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	call divider(kcell)
enddo
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
!	write(*,'(a,2e12.3)') 'Vdivide0,tgrowth: ',Vdivide0,divide_time_mean(1) - mitosis_duration
!	write(*,'(a,3e12.3)') 'r_mean, metab, dVdt: ',r_mean, metab, dVdt
endif
!if (suppress_growth) then	! for checking solvers
!	dVdt = 0
!endif
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
! A single cell is replaced by two.
!-----------------------------------------------------------------------------------------
subroutine divider(kcell1)
integer :: kcell1
integer :: kcell2, ityp, nbrs0
real(REAL_KIND) :: r(3), c(3), cfse0, cfse1
type(cell_type), pointer :: cp1, cp2

!write(*,*) 'divider:'
!write(logmsg,*) 'divider: ',kcell1 
!call logger(logmsg)
tnow = istep*DELTA_T
cp1 => cell_list(kcell1)
nlist = nlist + 1
ncells = ncells + 1
ityp = cp1%celltype
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
kcell2 = nlist
write(nflog,*) 'divider: ',kcell1,kcell2,cp1%ID
cp2 => cell_list(kcell2)
cp1%state = ALIVE
cp1%V = cp1%V/2
cp1%site = cp1%centre(:,1)/DELTA_X + 1
cp1%d = 0
cp1%birthtime = tnow
cp1%V_divide = get_divide_volume()
cp1%d_divide = (3*cp1%V_divide/PI)**(1./3.)
cp1%mitosis = 0
cfse0 = cp1%CFSE
cp1%CFSE = generate_CFSE(cfse0/2)
cfse1 = cfse0 - cp1%CFSE

!cp1%drugA_tag = .false.
!cp1%drugB_tag = .false.
cp1%drug_tag = .false.
cp1%anoxia_tag = .false.
cp1%t_hypoxic = 0

nbrs0 = cp1%nbrs
cp1%nbrs = nbrs0 + 1
cp1%nbrlist(cp1%nbrs)%indx = kcell2
cp1%nbrlist(cp1%nbrs)%contact = .false.
cp1%nbrlist(cp1%nbrs)%contact(1,1) = .true.

cp2%ID = cp1%ID
cp2%celltype = cp1%celltype
cp2%state = ALIVE
cp2%V = cp1%V
cp2%radius(1) = cp1%radius(2)
cp2%centre(:,1) = cp1%centre(:,2)
cp2%site = cp2%centre(:,1)/DELTA_X + 1
cp2%d = 0
cp2%birthtime = tnow
cp2%V_divide = get_divide_volume()
cp2%d_divide = (3*cp2%V_divide/PI)**(1./3.)
cp2%mitosis = 0
cp2%CFSE = cfse1

cp2%ID = cp1%ID
cp2%p_death = cp1%p_death
cp2%radiation_tag = cp1%radiation_tag
!cp2%drugA_tag = .false.
!cp2%drugB_tag = .false.
cp2%drug_tag = .false.
cp2%anoxia_tag = .false.
cp2%t_hypoxic = 0

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
