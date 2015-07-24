module reaction

use global
use chemokine

implicit none

contains

#if 0
!----------------------------------------------------------------------------------
subroutine f_rkc(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: i, k, ie, ki, kv, nextra, nintra, ichemo, site(3), kcell, ict, ith, Ng
real(REAL_KIND) :: dCsum, dCdiff, dCreact,  DX2, DX3, vol_cm3, val, Cin(MAX_CHEMO), Cex
real(REAL_KIND) :: decay_rate, dc1, dc6, cbnd, yy, C, membrane_kin, membrane_kout, membrane_flux, area_factor
logical :: bnd, dbug
logical :: TPZ_metabolised(MAX_CELLTYPES,0:2), DNB_metabolised(MAX_CELLTYPES,0:2)
real(REAL_KIND) :: metab, dMdt, KmetC
logical :: intracellular, cell_exists
logical :: use_actual_cell_volume = .false.

ichemo = icase
if (ichemo == GLUCOSE) then
	Ng = chemo(GLUCOSE)%Hill_N
endif

DX2 = DELTA_X*DELTA_X
decay_rate = chemo(ichemo)%decay_rate
dc1 = chemo(ichemo)%diff_coef/DX2
dc6 = 6*dc1 + decay_rate
cbnd = BdryConc(ichemo,t_simulation)
membrane_kin = chemo(ichemo)%membrane_diff_in
membrane_kout = chemo(ichemo)%membrane_diff_out
TPZ_metabolised(:,:) = (TPZ%Kmet0(:,:) > 0)	
DNB_metabolised(:,:) = (DNB%Kmet0(:,:) > 0)	

do i = 1,neqn
	yy = y(i)
	if (isnan(yy)) then
		write(nflog,*) 'f_rkc: isnan: ',i,ichemo,yy
		stop
	endif
    if (ODEdiff%vartype(i) == EXTRA) then
        intracellular = .false.
		vol_cm3 = Vsite_cm3
        Cex = yy
        cell_exists = .false.
        if (i < neqn) then
            if (ODEdiff%vartype(i+1) == INTRA) then
                cell_exists = .true.
				kcell = ODEdiff%cell_index(i+1)		! for access to cell-specific parameters (note that the intra variable follows the extra variable)
				if (use_actual_cell_volume) then
					vol_cm3 = Vsite_cm3 - Vcell_cm3*cell_list(kcell)%volume	! accounting for cell volume change
				else
	                vol_cm3 = Vextra_cm3									! for now, ignoring cell volume change!!!!!
	            endif
				area_factor = (cell_list(kcell)%volume)**(2./3.)
	            Cin = allstate(i+1,:)
	            Cin(ichemo) = y(i+1)
	        endif
	    endif
	else
        intracellular = .true.
	    kcell = ODEdiff%cell_index(i)					! for access to cell-specific parameters
	    if (use_actual_cell_volume) then
			vol_cm3 = Vcell_cm3*cell_list(kcell)%volume	! accounting for cell volume change
		else
			vol_cm3 = Vsite_cm3 - Vextra_cm3			! for now, ignoring cell volume change!!!!!
		endif
		area_factor = (cell_list(kcell)%volume)**(2./3.)
        Cex = y(i-1)
	    Cin = allstate(i,:)
	    Cin(ichemo) = yy
	    ict = cell_list(kcell)%celltype
	endif
	if (.not.intracellular) then
		! Need to check diffusion eqtn. when Vextra_cm3 < Vsite_cm3 = DX^3 !!!!!!!!!!!!!!!!!!!!!!
	    dCsum = 0
	    do k = 1,7
		    kv = ODEdiff%icoef(i,k)
		    if (k == 1) then
			    dCdiff = -dc6
		    else
			    dCdiff = dc1
		    endif
		    if (kv == OUTSIDE_TAG) then
			    val = cbnd
		    elseif (kv == UNREACHABLE_TAG) then
				val = y(ODEdiff%icoef(i,1))		! reflect concentration --> no flux
		    else
			    val = y(kv)
		    endif
		    dCsum = dCsum + dCdiff*val
	    enddo
	    if (cell_exists) then
!			membrane_flux = chemo(ichemo)%membrane_diff*(Cex - Cin(ichemo))*Vsite_cm3	! just a scaling.  We should account for change in surface area
			membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*Cin(ichemo))*Vsite_cm3	! just a scaling.
		else
		    membrane_flux = 0
		endif
    	dydt(i) = dCsum - membrane_flux/vol_cm3
	else
		C = Cin(ichemo)
!		membrane_flux = chemo(ichemo)%membrane_diff*(Cex - C)*Vsite_cm3		! just a scaling.  We should account for change in surface area
		membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)*Vsite_cm3	! just a scaling.
		dCreact = 0
		select case (ichemo)
		case (OXYGEN)
			metab = O2_metab(C)
			dCreact = (-metab*chemo(ichemo)%max_cell_rate*1.0e6 + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
		case (GLUCOSE)
			metab = C**Ng/(chemo(ichemo)%MM_C0**Ng + C**Ng)
			dCreact = (-metab*chemo(ichemo)%max_cell_rate*1.0e6 + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
		case (TRACER)
			dCreact = membrane_flux/vol_cm3
		case (TPZ_DRUG)
		    if (TPZ_metabolised(ict,0) .and. C > 0) then
				KmetC = TPZ%Kmet0(ict,0)*C
				if (TPZ%Vmax(ict,0) > 0) then
					KmetC = KmetC + TPZ%Vmax(ict,0)*C/(TPZ%Km(ict,0) + C)
				endif
				dCreact = -(1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + Cin(OXYGEN)))*KmetC
			endif
			dCreact = dCreact + membrane_flux/vol_cm3
		case (TPZ_DRUG_METAB_1)
			if (TPZ_metabolised(ict,0) .and. Cin(TPZ_DRUG) > 0) then
				dCreact = (1 - TPZ%C2(ict,0) + TPZ%C2(ict,0)*TPZ%KO2(ict,0)/(TPZ%KO2(ict,0) + Cin(OXYGEN)))*TPZ%Kmet0(ict,0)*Cin(TPZ_DRUG)
			endif
			if (TPZ_metabolised(ict,1) .and. C > 0) then
				dCreact = dCreact - (1 - TPZ%C2(ict,1) + TPZ%C2(ict,1)*TPZ%KO2(ict,1)/(TPZ%KO2(ict,1) + Cin(OXYGEN)))*TPZ%Kmet0(ict,1)*C
			endif
			dCreact = dCreact + membrane_flux/vol_cm3
		case (TPZ_DRUG_METAB_2)
			if (TPZ_metabolised(ict,1) .and. Cin(TPZ_DRUG_METAB_1) > 0) then
				dCreact = (1 - TPZ%C2(ict,1) + TPZ%C2(ict,1)*TPZ%KO2(ict,1)/(TPZ%KO2(ict,1) + Cin(OXYGEN)))*TPZ%Kmet0(ict,1)*Cin(TPZ_DRUG_METAB_1)
			endif
			if (TPZ_metabolised(ict,2) .and. C > 0) then
				dCreact = dCreact - (1 - TPZ%C2(ict,2) + TPZ%C2(ict,2)*TPZ%KO2(ict,2)/(TPZ%KO2(ict,2) + Cin(OXYGEN)))*TPZ%Kmet0(ict,2)*C
			endif
			dCreact = dCreact + membrane_flux/vol_cm3
		case (DNB_DRUG)
		    if (DNB_metabolised(ict,0) .and. C > 0) then
				KmetC = DNB%Kmet0(ict,0)*C
				if (DNB%Vmax(ict,0) > 0) then
					KmetC = KmetC + DNB%Vmax(ict,0)*C/(DNB%Km(ict,0) + C)
				endif
				dCreact = -(1 - DNB%C2(ict,0) + DNB%C2(ict,0)*DNB%KO2(ict,0)/(DNB%KO2(ict,0) + Cin(OXYGEN)))*KmetC
			endif
			dCreact = dCreact + membrane_flux/vol_cm3
		case (DNB_DRUG_METAB_1)
			if (DNB_metabolised(ict,0) .and. Cin(DNB_DRUG) > 0) then
				dCreact = (1 - DNB%C2(ict,0) + DNB%C2(ict,0)*DNB%KO2(ict,0)/(DNB%KO2(ict,0) + Cin(OXYGEN)))*DNB%Kmet0(ict,0)*Cin(DNB_DRUG)
			endif
			if (DNB_metabolised(ict,1) .and. C > 0) then
				dCreact = dCreact - (1 - DNB%C2(ict,1) + DNB%C2(ict,1)*DNB%KO2(ict,1)/(DNB%KO2(ict,1) + Cin(OXYGEN)))*DNB%Kmet0(ict,1)*C
			endif
			dCreact = dCreact + membrane_flux/vol_cm3
		case (DNB_DRUG_METAB_2)
			if (DNB_metabolised(ict,1) .and. Cin(DNB_DRUG_METAB_1) > 0) then
				dCreact = (1 - DNB%C2(ict,1) + DNB%C2(ict,1)*DNB%KO2(ict,1)/(DNB%KO2(ict,1) + Cin(OXYGEN)))*DNB%Kmet0(ict,1)*Cin(DNB_DRUG_METAB_1)
			endif
			if (DNB_metabolised(ict,2) .and. C > 0) then
				dCreact = dCreact - (1 - DNB%C2(ict,2) + DNB%C2(ict,2)*DNB%KO2(ict,2)/(DNB%KO2(ict,2) + Cin(OXYGEN)))*DNB%Kmet0(ict,2)*C
			endif
			dCreact = dCreact + membrane_flux/vol_cm3
		end select
	    dydt(i) = dCreact - yy*decay_rate
	endif
enddo
end subroutine
#endif



!-----------------------------------------------------------------------------------------
! Mass flux dMdt has units mol/s, therefore %membrane_diff_* has units of flow: cm^3/s
! Note that we will separate Kin and Kout.
! chemo(:)%max_cell_rate has units mol/cell/s
! CHANGED!!!!!!!!!!
! Now mass flux has unit mumol/s
! This code needs to be fixed or removed
!-----------------------------------------------------------------------------------------
subroutine cell_reactor(cp, dt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt
real(REAL_KIND) :: dMdt(NCONST), C, metab, vol_cm3, dCdt(NCONST)
integer :: ichemo, Ng

!dMdt(:) = Kmemb(:)*(cp%Cave(:) - cp%Cin(:))	! flux into the cell (mol/s)
dMdt(:) = chemo(1:NCONST)%membrane_diff_in*cp%Cex(:) - chemo(1:NCONST)%membrane_diff_out*cp%Cin(:)	! flux into the cell (mol/s)
write(*,'(a,2e12.3)') 'Cin:  ',cp%Cin(:)
write(*,'(a,2e12.3)') 'Cex:  ',cp%Cex(:)
write(*,'(a,2e12.3)') 'Kin:  ',chemo(1:NCONST)%membrane_diff_in
write(*,'(a,2e12.3)') 'Kout: ',chemo(1:NCONST)%membrane_diff_out
write(*,'(a,2e12.3)') 'dMdt: ',dMdt(:)
cp%dMdt = dMdt
vol_cm3 = um3_cm3*cp%V	! um^3 -> cm^3
do ichemo = 1,NCONST
	C = cp%Cin(ichemo)
	select case (ichemo)
		case (OXYGEN)
			metab = O2_metab(C)
			dMdt(ichemo) = dMdt(ichemo) - metab*chemo(ichemo)%max_cell_rate
		case (GLUCOSE)
			Ng = chemo(ichemo)%Hill_N
			metab = C**Ng/(chemo(ichemo)%MM_C0**Ng + C**Ng)
			dMdt(ichemo) = dMdt(ichemo) - metab*chemo(ichemo)%max_cell_rate
	end select
!	write(*,'(a,i2,5e12.3)') 'ichemo: ',ichemo,C,metab,chemo(ichemo)%max_cell_rate,metab*chemo(ichemo)%max_cell_rate,vol_cm3
enddo
dCdt = dMdt(:)*1.0e6/vol_cm3	! mol/cm^3/s -> mM/s
!write(*,'(a,4e12.3)') 'Cin: ',cp%Cin(:),dCdt(:)
!cp%Cin(:) = cp%Cin(:) + dt*dCdt		
!stop
end subroutine

end module