! Chemokine data

module chemokine

use global

implicit none

type chemokine_type
	character(24) :: name
	logical :: used
	logical :: present
	logical :: constant
	real(REAL_KIND) :: bdry_rate
	real(REAL_KIND) :: bdry_conc
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: membrane_diff_in
	real(REAL_KIND) :: membrane_diff_out
	logical :: decay
	real(REAL_KIND) :: halflife			! h
	real(REAL_KIND) :: decay_rate		! /sec
	real(REAL_KIND) :: max_cell_rate	! Vmax
	real(REAL_KIND) :: MM_C0			! Km
	real(REAL_KIND) :: Hill_N
	real(REAL_KIND) :: medium_diff_coef	! diffusion coefficient in the medium
	real(REAL_KIND) :: diff_reduction_factor
	real(REAL_KIND) :: medium_Cext		! far-field concentration
!	Fine grid
	real(REAL_KIND), allocatable :: Cprev(:,:,:)
	real(REAL_KIND), allocatable :: Fprev(:,:,:)
!	Coarse grid
	real(REAL_KIND), allocatable :: Cave_b(:,:,:)
	real(REAL_KIND), allocatable :: Cprev_b(:,:,:)
	real(REAL_KIND), allocatable :: Fprev_b(:,:,:), Fcurr_b(:,:,:)
end type

type(chemokine_type), target :: chemo(MAX_CHEMO)

integer :: nchemo, chemomap(MAX_CHEMO)

contains

!----------------------------------------------------------------------------------------
! Convert halflife in hours to a decay rate /sec
!----------------------------------------------------------------------------------------
real(REAL_KIND) function DecayRate(halflife)
real(REAL_KIND) :: halflife

if (halflife == 0) then		! No decay
	DecayRate = 0
else
	DecayRate = log(2.0)/(halflife*60*60)
endif
end function

!----------------------------------------------------------------------------------------
! Note units:
! distance		cm
! volume		cm^3	
! time			s
! diff coeff	cm^2.s^-1
! mass			mol
! concentration	mM where 1 mM = 1.0e-6 mol.cm^-3 = 1 mumol.cm^-3
! consumption	mol.cell^-1.s^-1
! production	mol.cell^-1.s^-1
! flux			mumol.s^-1
!----------------------------------------------------------------------------------------
subroutine SetupChemo
integer :: ichemo

chemo(OXYGEN)%name = 'Oxygen'
chemo(GLUCOSE)%name = 'Glucose'
chemo(TRACER)%name = 'Tracer'
chemo(OXYGEN)%decay_rate = 0
chemo(GLUCOSE)%decay_rate = 0
chemo(TRACER)%decay_rate = 0

!chemo(DRUG_A)%name = 'Drug_A'
!chemo(DRUG_B)%name = 'Drug_B'
do ichemo = 1,MAX_CHEMO
	chemo(ichemo)%present = .false.
	if (chemo(ichemo)%used) then
		if (ichemo == OXYGEN .or. ichemo == GLUCOSE .or. ichemo == TRACER) then
			chemo(ichemo)%present = .true.
		endif
	endif
enddo
! constituent fields, fine and coarse grids
do ichemo = 1,2
	if (allocated(chemo(ichemo)%Cprev)) deallocate(chemo(ichemo)%Cprev)
	if (allocated(chemo(ichemo)%Fprev)) deallocate(chemo(ichemo)%Fprev)
	if (allocated(chemo(ichemo)%Cave_b)) deallocate(chemo(ichemo)%Cave_b)
	if (allocated(chemo(ichemo)%Cprev_b)) deallocate(chemo(ichemo)%Cprev_b)
	if (allocated(chemo(ichemo)%Fprev_b)) deallocate(chemo(ichemo)%Fprev_b)
	if (allocated(chemo(ichemo)%Fcurr_b)) deallocate(chemo(ichemo)%Fcurr_b)
	allocate(chemo(ichemo)%Cprev(NX,NY,NZ))
	allocate(chemo(ichemo)%Fprev(NX,NY,NZ))
	allocate(chemo(ichemo)%Cave_b(NXB,NYB,NZB))
	allocate(chemo(ichemo)%Cprev_b(NXB,NYB,NZB))
	allocate(chemo(ichemo)%Fprev_b(NXB,NYB,NZB))
	allocate(chemo(ichemo)%Fcurr_b(NXB,NYB,NZB))
enddo
chemo(OXYGEN)%diff_reduction_factor = 0.3
chemo(GLUCOSE)%diff_reduction_factor = 0.5

call SetupChemomap
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetupChemomap
integer :: ichemo

nchemo = 0
do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%present) then
		nchemo = nchemo + 1
		chemomap(nchemo) = ichemo
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
! Consumption rate = Rmax*C/(MM_C0 + C)  (Michaelis-Menten)
! where Km = MM_C0 >= Rmax*T*10^6/Vextra_cm3 is the O2 concentration at which cell uptake is halved,
! and T is the maximum expected time step
! This is to ensure that all the O2 in the site is not depleted in a time step, making
! O2 conc go negative.  For now it seems reasonable to assume that glucose uptake
! varies in proportion to O2 uptake, i.e. we only need Oxygen M-M
!----------------------------------------------------------------------------------------
subroutine SetMMParameters

chemo(OXYGEN)%MM_C0 = 0.00133		! 1 mmHg = 1.33 uM (Kevin suggests 1 - 2 uM)
chemo(GLUCOSE)%MM_C0 = chemo(OXYGEN)%MM_C0*chemo(GLUCOSE)%max_cell_rate/chemo(OXYGEN)%max_cell_rate
write(logmsg,'(a,e12.4)') 'Oxygen MM_C0: ',chemo(OXYGEN)%MM_C0
call logger(logmsg)
write(logmsg,'(a,e12.4)') 'Glucose MM_C0: ',chemo(GLUCOSE)%MM_C0
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
! Use the "soft landing" option for Hill_N = 1 if MM_threshold = 0
!----------------------------------------------------------------------------------
real(REAL_KIND) function O2_metab(C)
!real(REAL_KIND) function metabolic_rate(ichemo,C)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: metab
real(REAL_KIND) :: deltaC_soft = 0.0, C1_soft = 0.0, k_soft = 0.0

ichemo = OXYGEN
if (ichemo == OXYGEN) then
	if (chemo(ichemo)%Hill_N == 2) then
		if (C > 0) then
			metab = C*C/(chemo(ichemo)%MM_C0*chemo(ichemo)%MM_C0 + C*C)
		else
			metab = 0
		endif
	else
		if (MM_THRESHOLD > 0) then
			if (C > C1_soft) then
				metab = (C-deltaC_soft)/(chemo(ichemo)%MM_C0 + C - deltaC_soft)
			elseif (C > 0) then
				metab = k_soft*C*C
			else
				metab = 0
			endif
		else
			if (C > 0) then
				metab = C/(chemo(ichemo)%MM_C0 + C)
			else
				metab = 0
			endif
		endif
	endif
endif
O2_metab = metab
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!subroutine AllocateConcArrays
!integer :: ic
!
!write(logmsg,*) 'AllocateConcArrays'
!call logger(logmsg)
!do ic = 1,MAX_CHEMO
!	if (chemo(ic)%used) then
!		if (allocated(chemo(ic)%Cin)) then
!			call logger("chemo(ic)%Cin already allocated")
!			deallocate(chemo(ic)%Cin)
!!			stop
!		endif
!		allocate(chemo(ic)%Cin(NX,NY,NZ))
!		if (allocated(chemo(ic)%grad)) then
!			call logger("chemo(ic)%grad already allocated")
!			deallocate(chemo(ic)%grad)
!!			stop
!		endif
!		allocate(chemo(ic)%grad(3,NX,NY,NZ))
!		chemo(ic)%Cin = chemo(ic)%bdry_conc	
!	endif
!enddo
!write(logmsg,*) 'did AllocateConcArrays'
!call logger(logmsg)
!end subroutine

!----------------------------------------------------------------------------------------
! Set the concentrations at a site when it is on the boundary
!----------------------------------------------------------------------------------------
!subroutine SetBdryConcs(site)
!integer :: site(3)
!integer :: ichemo
!
!do ichemo = 1,MAX_CHEMO
!	if (chemo(ichemo)%used) then
!		chemo(ichemo)%Cin(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
!	endif
!enddo
!end subroutine



end module