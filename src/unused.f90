!-----------------------------------------------------------------------------------------
! Actually place a single cell
!-----------------------------------------------------------------------------------------
subroutine PlaceCells1(ok)
logical :: ok
real(REAL_KIND) :: r(3), c(3)
integer, parameter :: ncx=10
integer :: kcell, icx, icy, icz
type(cell_type), pointer :: cp

nlist = 0
ncells_mphase = 0
do icx = 1,ncx
	do icy = 1,ncx
		do icz = 1,ncx
			nlist = nlist + 1
			cp=>cell_list(nlist)
			cp%state = ALIVE
			if (MITOSIS_MODE == TERMINAL_MITOSIS) then
				cp%Iphase = .true.
				cp%nspheres = 1
			else
				cp%Iphase = .false.
				cp%nspheres = 2
				ncells_mphase = ncells_mphase + 1
			endif
			cp%radius(1) = Raverage
			cp%V = (4.*PI/3.)*Raverage**3
			cp%site = [NX/2-ncx/2+icx,NY/2-ncx/2+icy,NZ/2-ncx/2+icz]
			cp%centre(:,1) = [(cp%site(1)-0.5)*DELTA_X,(cp%site(2)-0.5)*DELTA_X,(cp%site(3)-0.5)*DELTA_X]
			!cp%centre(:,1) = [0,0,0]
			cp%d = 0
			cp%birthtime = 0
			cp%growthrate = test_growthrate
			cp%V_divide = Vdivide0
			!cp2%V_divide = get_divide_volume()
			cp%d_divide = (3*cp%V_divide/PI)**(1./3.)
			cp%mitosis = 0
			call get_random_vector3(r)	! set initial axis direction
			cp%d = 0.1*small_d
			c = cp%centre(:,1)
			cp%centre(:,1) = c + (cp%d/2)*r
			cp%centre(:,2) = c - (cp%d/2)*r
			cp%nbrs = 0
			cp%Cex = Cextra(1,1,1,:)	! initially the concentrations are the same everywhere
			cp%Cin = cp%Cex
		enddo
	enddo
enddo
ncells = nlist
end subroutine

!-------------------------------------------------------------------------------------------
! Flux into the cell.  The membrane flux is the equilibrium value, i.e. 
! rate of consumption is balanced by the rate of mass transport across the membrane.  
! The rate is positive when it adds to the intracellular concentration.
! Should dVindt be included?? I think not.
!-------------------------------------------------------------------------------------------
subroutine getF_SS1(ichemo,W,F)
integer :: ichemo
real(REAL_KIND) :: W(:,:,:), F(:,:,:)
real(REAL_KIND) :: Km, Cin, Cex, Vin, Fsum
real(REAL_KIND) :: Kin, Kout
integer :: ix, iy, iz

Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
!Fsum = 0
!!$omp parallel do private(iy, ix, Cex, Cin)
do iz = 1,NZ
do iy = 1,NY
do ix = 1,NX
	if (ngcells(ix,iy,iz) > 0) then
		Cex = W(ix,iy,iz)
		Cin = getCin_SS(ichemo,Vin,Cex)
		F(ix,iy,iz) = (Kin*Cex - Kout*Cin)	! + Cex*dVindt(ix,iy,iz) - CexVex**decayrate_ex
!		if (ix==NX/2 .and. iy==NY/2 .and. iz== NZ/2) write(*,*) 'dC: ',Cin-Cex,F(ix,iy,iz)
!		write(*,'(3i6,4e12.3)') ix,iy,iz,Cex,Cin,Vex,F(ix,iy,iz)
	else
		F(ix,iy,iz) = 0		! - Cex*Vex*decayrate_ex
	endif
!	Fsum = Fsum + F(ix,iy,iz)
enddo
enddo
enddo
!!$omp end parallel do
!write(*,*) 'Fsum: ',Fsum
end subroutine

!-------------------------------------------------------------------------------------------
! This valid for steady-state with any constituent.
! Includes decay and recognises that the PDE actually solves for Cave:
! Cave = (Vex.Cex + Vin.Cin)/(Vex+Vin)
! The input Cave is the average grid concentration, i.e. what was previously called Cex elsewhere.
!-------------------------------------------------------------------------------------------
function getCin_SS1(ichemo, Cave, Vin) result(Cin)
integer :: ichemo
real(REAL_KIND) :: Cave, Vin, Cin
real(REAL_KIND) :: Vex, Vg, Kin, Kout, Kmax, Kdecay, K1, K2, C0, a, b, c, D, r(3)
integer :: i, n

Vg = dx3
Vex = Vg - Vin
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_in
Kmax = chemo(ichemo)%max_cell_rate
Kdecay = chemo(ichemo)%decay_rate
C0 = chemo(ichemo)%MM_C0
K1 = (Kout + Kdecay + Kin*Vin/Vex)/Kmax
K2 = (Kin*Vg*Cave)/(Kmax*Vex)
if (chemo(ichemo)%Hill_N == 2) then
	a = (1 - K2)/K1
	b = C0*C0
	c = -b*K2/K1
	call cubic_roots(a,b,c,r,n)
	if (n == 1) then
		Cin = r(1)
	else
		n = 0
		do i = 1,3
			if (r(i) > 0) then
				n = n+1
				Cin = r(i)
			endif
		enddo
		if (n > 1) then
			write(*,*) 'getCin_SS: two roots > 0: ',r
			stop
		endif
	endif
elseif (chemo(ichemo)%Hill_N == 1) then
	b = (K1*C0 - K2 + 1)/K1
	c = -C0*K2/K1
	D = sqrt(b*b - 4*c)
	Cin = (D - b)/2
endif
if (Cave < 0) then
	write(*,*) 'getCin_SS: Cave,Cin: ',Cave,Cin
	stop
endif
end function

!-------------------------------------------------------------------------------------------
! This steady-state solution is valid for OXYGEN only.
! No need to separate Kin and Kout
!-------------------------------------------------------------------------------------------
!function getCin_O2(Cextra) result(Cin)
!real(REAL_KIND) :: Cextra, Cin
!real(REAL_KIND) :: K1, K2, K2K1, C0, a, b, cc, D, r(3)
!integer :: ichemo, i, n
!
!ichemo = OXYGEN
!K1 = chemo(ichemo)%membrane_diff_in
!K2 = chemo(ichemo)%max_cell_rate
!K2K1 = K2/K1
!C0 = chemo(ichemo)%MM_C0
!if (chemo(ichemo)%Hill_N == 2) then
!	a = K2K1 - Cextra
!	b = C0*C0
!	cc = -b*Cextra
!	call cubic_roots(a,b,cc,r,n)
!	if (n == 1) then
!		Cin = r(1)
!	else
!		n = 0
!		do i = 1,3
!			if (r(i) > 0) then
!				n = n+1
!				Cin = r(i)
!			endif
!		enddo
!		if (n > 1) then
!			write(*,*) 'getCin_O2: two roots > 0: ',r
!			stop
!		endif
!	endif
!elseif (chemo(ichemo)%Hill_N == 1) then
!	b = K2K1 + C0 - Cextra
!	cc = -C0*Cextra
!	D = sqrt(b*b - 4*cc)
!	Cin = (D - b)/2
!endif
!if (Cextra < 0) then
!	write(*,*) 'getCin_O2: Cex,Cin: ',Cextra,Cin
!	stop
!endif
!end function
