! To determine distribution of colony size
! Note: this was developed originally for the on-lattice model, for which %volume, %divide_volume are normalised!
! This is the off-lattice version, volumes are NOT normalised

module colony

use global
implicit none

integer, parameter :: n_colony_days=10
integer :: nmax
type(cell_type), target, allocatable :: ccell_list(:)

contains

!---------------------------------------------------------------------------------------------------
! Simulate fate of cells grown with no nutrient constraints.
! The only determinants of the colony size for a cell are (considering radiation only):
! volume
! divide_time_mean(ityp)
! radiation_tag
! p_rad_death
! growth_delay
! G2_M
! dt_delay
! t_growth_delay_end
! N_delayed_cycles_left
! For simplicity, assume that each cell has the mean growth rate, and that for a cell that is not
! at the G2/M checkpoint, the fraction of time through the growth period is V0/divide_volume, where
! V0 = volume
! and
! divide_volume = Vdivide0 + dVdivide*(2*R-1)

!---------------------------------------------------------------------------------------------------
subroutine make_colony_distribution(tnow)
real(REAL_KIND) :: tnow
real(REAL_KIND) :: V0, dVdt, dt, t, tend
integer, parameter :: ndist=40
real(REAL_KIND) :: ddist, dist(ndist)
integer :: kcell, ityp, n, idist, ncycmax, ntot
type (cell_type), pointer :: cp

write(nflog,*) 'make_colony_distribution: nlist: ',nlist
ncycmax = 24*3600*n_colony_days/divide_time_mean(1) + 3
nmax = 2**ncycmax
allocate(ccell_list(nmax))
ddist = 50
dist = 0
ntot = 0
tend = tnow + 10*24*3600
do kcell = 1, nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	if (.not.cp%G2_M) then
		! approximate time to reach divide volume
		if (use_divide_time_distribution) then
			cp%divide_volume = Vdivide0
		endif
		V0 = cp%V	! actual volume cm3
		dt = ((cp%divide_volume - V0)/cp%divide_volume)*divide_time_mean(ityp)
	else
		dt = 0
	endif
	t = tnow + dt
!	cp%V = cp%divide_volume/Vcell_mm3
	cp%V = cp%divide_volume	! actual volume cm3
	! Now simulate colony growth
	call make_colony(kcell,t,tend,n)
	write(nflog,*) 'cell: n: ',kcell,n
	ntot = ntot + n
	idist = n/ddist + 1
	dist(idist) = dist(idist) + 1
enddo
dist = dist/sum(dist)
write(nflog,*) 'Colony size distribution: ', nlist,ntot
write(nfout,*) 'Colony size distribution: ', nlist,ntot
do idist = 1,ndist
	write(nfout,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
enddo
deallocate(ccell_list)
end subroutine

!---------------------------------------------------------------------------------------------------
! The cell is at the point of division - possibly G2_M (arrested at G2/M checkpoint)
! For now only radiation tagging is handled
! Growth rate dVdt (mean) is used only to estimate the time of next division
!---------------------------------------------------------------------------------------------------
subroutine make_colony(kcell,t,tend,n)
integer :: kcell, n
real(REAL_KIND) :: t, tend, dt 
integer :: icell, ityp, nl, nl0, kpar=0
real(REAL_KIND) :: V0, Tdiv0, r_mean, c_rate, dVdt, Tmean, R
logical :: ok
type (cell_type), pointer :: cp

!write(*,'(a,i6,2f8.0)') 'make_colony: ',kcell,t,tend
ccell_list(1) = cell_list(kcell)
ccell_list(1)%anoxia_tag = .false.
ccell_list(1)%drug_tag = .false.
dt = 3600	! use big time steps, 1h
ityp = cell_list(kcell)%celltype
Tdiv0 = divide_time_mean(ityp)
r_mean = Vdivide0/Tdiv0
c_rate = log(2.0)/Tdiv0
ccell_list(1)%V = cell_list(kcell)%divide_volume
ccell_list(1)%t_divide_next = t
nl = 1
do while (t < tend)
	t = t + dt
	nl0 = nl
	do icell = 1,nl0
		cp => ccell_list(icell)
		if (cp%state == DEAD) cycle
		if (t >= cp%t_divide_next) then
!			if (kcell == nlist) write(*,*) 't_divide_next: ',icell,cp%t_divide_next
			! Time to divide, unless delayed
			if (cp%growth_delay) then
				if (cp%G2_M) then
					if (t > cp%t_growth_delay_end) then
						cp%G2_M = .false.
					else
						cycle
					endif
				else
					cp%t_growth_delay_end = t + cp%dt_delay
					cp%G2_M = .true.
					cycle
				endif
			endif
			if (cp%radiation_tag) then
				R = par_uni(kpar)
				if (R < cp%p_rad_death) then
					! Cell dies
					cp%state = DEAD
					cycle
				endif
			endif
			! Can divide
!			write(*,*) 'divides: ',icell
			nl = nl+1
			if (nl > nmax) then
				write(nflog,*) 'nmax exceeded: ',nmax
				stop
			endif
			call CloneColonyCell(icell,nl,t,ok)
		endif
		if (use_V_dependence) then
			dVdt = c_rate*cp%V
		else
			dVdt = r_mean
		endif
		cp%V = cp%V + dt*dVdt
	enddo
enddo
n = 0
do icell = 1,nl
	if (ccell_list(icell)%state /= DEAD) n = n+1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CloneColonyCell(kcell0,kcell1,tnow,ok)
integer :: kcell0, kcell1, ityp, idrug
logical :: ok
integer :: kpar = 0
real(REAL_KIND) :: tnow, R, Tdiv, Tdiv0

ok = .true.
ityp = ccell_list(kcell0)%celltype
Tdiv0 = divide_time_mean(ityp)
ccell_list(kcell0)%V = ccell_list(kcell0)%V/2
ccell_list(kcell0)%generation = ccell_list(kcell0)%generation + 1
if (ccell_list(kcell0)%growth_delay) then
	ccell_list(kcell0)%N_delayed_cycles_left = ccell_list(kcell0)%N_delayed_cycles_left - 1
	ccell_list(kcell0)%growth_delay = (ccell_list(kcell0)%N_delayed_cycles_left > 0)
endif
ccell_list(kcell0)%G2_M = .false.
ccell_list(kcell0)%t_divide_last = tnow
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	ccell_list(kcell0)%t_divide_next = tnow + Tdiv
else
	R = par_uni(kpar)
	ccell_list(kcell0)%t_divide_next = tnow + Tdiv0*(1 + (2*dVdivide/Vdivide0)*(2*R-1))
endif
!write(*,*) 'divide time: ',Tdiv0*(1 + (2*dVdivide/Vdivide0)*(2*R-1))/3600
ccell_list(kcell1)%celltype = ccell_list(kcell0)%celltype
ccell_list(kcell1)%state = ccell_list(kcell0)%state
ccell_list(kcell1)%V = ccell_list(kcell0)%V
ccell_list(kcell1)%generation = ccell_list(kcell0)%generation
ccell_list(kcell1)%ID = ccell_list(kcell0)%ID
ccell_list(kcell1)%p_rad_death = ccell_list(kcell0)%p_rad_death
ccell_list(kcell1)%p_drug_death = ccell_list(kcell0)%p_drug_death
ccell_list(kcell1)%radiation_tag = ccell_list(kcell0)%radiation_tag
ccell_list(kcell1)%anoxia_tag = .false.
!ccell_list(kcell1)%exists = .true.
!ccell_list(kcell1)%active = .true.
ccell_list(kcell1)%growth_delay = ccell_list(kcell0)%growth_delay
if (ccell_list(kcell1)%growth_delay) then
	ccell_list(kcell1)%dt_delay = ccell_list(kcell0)%dt_delay
	ccell_list(kcell1)%N_delayed_cycles_left = ccell_list(kcell0)%N_delayed_cycles_left
endif
ccell_list(kcell1)%G2_M = .false.
ccell_list(kcell1)%t_divide_last = tnow
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	ccell_list(kcell1)%t_divide_next = tnow + Tdiv
else
	R = par_uni(kpar)
	ccell_list(kcell1)%t_divide_next = tnow + Tdiv0*(1 + (2*dVdivide/Vdivide0)*(2*R-1))
endif
ccell_list(kcell1)%t_hypoxic = 0
!ccell_list(kcell1)%conc = ccell_list(kcell0)%conc
ccell_list(kcell1)%Cex = ccell_list(kcell0)%Cex
ccell_list(kcell1)%dCdt = ccell_list(kcell0)%dCdt
ccell_list(kcell1)%dMdt = ccell_list(kcell0)%dMdt
!ccell_list(kcell1)%M = ccell_list(kcell0)%M
end subroutine


end module
