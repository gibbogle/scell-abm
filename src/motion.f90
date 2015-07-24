! To simulate off-lattice cells as spheres, following Drasdo & Hohme 
!
! Biomechanics, includung cell motion, is treated by the methods described in:
!
! Metropolis N, Rosenbluth A W, Rosenbluth M N, Teller A H and Teller E 1953 Equation of state calculations by fast computing
! machines J. Chem. Phys. 21 1087–92
!
! Drasdo D, Kree R and McCaskill J S 1995 Monte-carlo approach to tissue-cell populations Phys. Rev. E 52 6635–57
!
! D. Drasdo and S. Hoehme, A single-cell based model to tumor growth in-vitro: Monolayers and spheroids. Phys. Biol. 2:133 (2005).
! 
!	The state S; of a cell i is described by its position, its age class, and its shape. The cell position is given
!	by the vector r, of its center of mass. For an I cell, the shape is characterized by its radius R, whereas for
!	an M cell it is the vector d(i), joining the centers of the corresponding dumbbell which fixes its shape. As the
!	radii R(1,i) = R(2,i) = R(i), of the spheres making up the dumbbell are equal and furthermore determined by the
!	length of d(i), we may also describe the shape by R(i) and a unit vector d(i), which fixes the orientation of an M cell in
!	space.
!	During a single updating step, only one randomly chosen cell is considered as active and tries to perform one
!	of the following actions: (i) action_I (migration, growth) for I cells, (ii) action_M (migration, rotation, deformation)
!	for M cells.
!	A single migration trial consists of a shift r(i) -> r(i) + dr(i), of cell i in a random direction with a step length
!	which is chosen at random out of an interval [O, dr_max].
!	A rotation trial is defined analogously. During a growth step, the radius of an I cell is increased by a random
!	amount dR in [O, dR_max] whereas during a deformation step the length of d(i) is increased and at the same time
!	the dumbbell radius is decreased such that the volume of the cell is kept constant. All these trials are described in
!	detail in the Appendix.
!	Whether cell i actually performs its chosen trial or not depends on its interaction with other cells. Following
!	the discussion of the preceding subsection, we model the interaction between two cells i and j by potentials
!	V(S(i), S(j)). Note that these interaction potentials may be different for different types of actions. We will discuss
!	the biological basis of this complication below.
!	The cell dynamics is based on the assumption that a trial is accepted according to the standard prescription
!	of the Metropolis algorithm [49], i.e., an action a is always accepted if the change of the state S(i) -> S^(i) of cell i
!	induced by the action, leads to a decrease of its interaction potential with other cells:
!       E(S(i)) = Sum {V(S(i),S(j))}, 
!   whereas it is only accepted with probability 
!       exp(—[E(S^(i))—E(S(i))]) if E(S(i)) is increased by the action. 

module motion

use global
implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine mover(ok)
logical :: ok
integer :: it, kcell, nits, k, i, axis, axismin
integer :: kpar=0
real(REAL_KIND) :: R1, c1(3), E0, E1, dE, v(3), dr(3), drot, p, Emin, drmin(3), drotmin, cc1(3,2)
logical :: Mphase
type(cell_type), pointer :: cp1

ok = .true.
if (ncells == 1) return
write(logmsg,*) 'mover: ',istep,ncells
call logger(logmsg)
do k = 1,ncells
	kcell = perm_index(k)
	cp1 => cell_list(kcell)
	Mphase = .not.cp1%Iphase
	write(logmsg,*) '--- cell: ',kcell,cp1%Iphase
	call logger(logmsg)
	E0 = Evalue(kcell,cp1%centre,.true.)
	if (E0 >= BIG/2) then
		write(logmsg,*) 'mover: bad E0: ',istep,kcell,E0
		call logger(logmsg)
		ok = .false.
!		ok = .true.	! try just going on, treat this cell as stuck - NOT A GOOD IDEA - cells pile up
		return
	endif
	Emin = 100*BIG
	do i = 1,ntries
		call get_random_vector3(v)
		dr = dr_max*v
		! For Mphase add random rotation here, e.g. about a randomly chosen axis (x,y,z)
		if (Mphase) then
			call get_random_drot(axis,drot)
		else
			axis = 0
			drot = 0
		endif
		write(logmsg,'(3f8.3,i6,f8.4)') dr,axis,drot
		call logger(logmsg)
		call transform_centres(cp1,dr,axis,drot,cc1)
		E1 = Evalue(kcell,cc1,.false.)
		if (E1 < Emin) then
			Emin = E1
			drmin = dr
			axismin = axis
			drotmin = drot
		endif
	enddo
	dE = Emin - E0
	write(logmsg,'(a,2f8.3)') 'E0, dE: ',E0,dE
	call logger(logmsg)
	if (dE >= 0) then
		p = exp(-dE)
		if (p > par_uni(kpar)) then
			call logger('no move')
			call logger('')
			cycle
		endif
	endif
	dr = drmin
	if (Mphase) then
		axis = axismin
		drot = drotmin
	else
		axis = 0
		drot = 0
	endif
	write(logmsg,'(a,3f6.2,i6,f8.4)') 'drmin,axismin,drotmin: ',dr,axis,drot
	call logger(logmsg)
	call transform_centres(cp1,dr,axis,drot,cp1%centre)
	call logger('')
	
	! Need to update contacts with neighbours here?
	call set_neighbour_contacts(cp1)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Need to sum over all interactions, max = 4
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function Evalue(kcell,cc1,current) result(E)
integer :: kcell
real(REAL_KIND) :: cc1(3,2)
logical :: current
logical :: contact(2,2)
type(cell_type), pointer :: cp1, cp2
real(REAL_KIND) :: R1, c1(3),R2, c2(3), v(3), d, dE
integer :: isphere1, nspheres1, isphere2, nspheres2
integer :: k, knbr
logical :: incontact, dbug

dbug = (istep > 10001750)
cp1 => cell_list(kcell)
if (cp1%Iphase) then
	nspheres1 = 1
else
	nspheres1 = 2
endif
if (dbug) write(*,*) 'kcell, nspheres1, nbrs: ',kcell,nspheres1,cp1%nbrs
E = 0
do k = 1,cp1%nbrs
	knbr = cp1%nbrlist(k)%indx
	cp2 => cell_list(knbr)
	if (cp2%Iphase) then
		nspheres2 = 1
	else
		nspheres2 = 2
	endif
	if (dbug) write(*,*) 'knbr,  nspheres2, nbrs: ',knbr,nspheres2,cp2%nbrs
	do isphere1 = 1,nspheres1
		R1 = cp1%radius(isphere1)
		c1 = cc1(:,isphere1)
		do isphere2 = 1,nspheres2
			R2 = cp2%radius(isphere2)
			c2 = cp2%centre(:,isphere2)
			incontact = cp1%nbrlist(k)%contact(isphere1,isphere2)
			v = c1 - c2
			d = sqrt(dot_product(v,v))
			if (dbug) write(*,'(2i4,8f6.1,L,f10.2)') isphere1,isphere2,R1,R2,c1,c2,incontact,d
            dE = Vfun(R1,c1,R2,c2,incontact)
            if (current .and. dE > BIG/2) then
				write(logmsg,'(a,L2,2x,2i3,2x,2i3,2f6.2,2x,3f6.2,2x,3f6.2,2x,3f6.2,2x,f6.2)')  &
					'Vfun: bad dE: ',current,kcell,knbr,isphere1,isphere2,R1,R2,c1,c2,v,d
				call logger(logmsg)
				stop
			endif
			E = E + dE
			if (dbug) write(*,'(a,e12.3)') 'dE: ',dE
!			contact(isphere1,isphere2) = incontact
		enddo
	enddo
enddo
end function

!-----------------------------------------------------------------------------------------
! dd is d~
! incontact conveys whether or not the neighbour was previously in contact.
! This affects the adhesion force, and interaction energy, when 1 < d/Rsum < k_detach
! Note that incontact may be changed here, both ways, implementing the hysteresis loop
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function Vfun(R1,c1,R2,c2,incontact)
real(REAL_KIND) :: R1, R2, c1(3), c2(3)
logical :: incontact
real(REAL_KIND) :: v(3), d, dd, Rsum, delta	!,Rc, d1, Ac

Vfun = 0
v = c2 - c1
d = sqrt(dot_product(v,v))
Rsum = R1+R2
dd = d - Rsum
delta = alpha_v*Rsum
if (d > k_detach*Rsum) then
	incontact = .false.
elseif (.not. incontact .and. dd < 0) then
	incontact = .true.
endif

if (dd < -delta) then
	Vfun = BIG
elseif (.not.incontact) then
	Vfun = 0
else
	Vfun = epsilon*((2*dd/delta + sqr_es_e)**2 - es_e) + shift
endif
!write(*,'(a,2f6.2)') 'd,dd: ',d,dd

end function

!-----------------------------------------------------------------------------------------
! Translate then rotate about the mid-point
!-----------------------------------------------------------------------------------------
subroutine transform_centres(cp1,dr,axis,drot,cc1)
type(cell_type), pointer :: cp1
real(REAL_KIND) :: dr(3), drot, cc1(3,2)
integer :: axis
real(REAL_KIND) :: cmid(3), c1(3), r2, r, theta, dx(3), alpha

cc1(:,1) = cp1%centre(:,1) + dr
if (.not.cp1%Iphase) then
	cc1(:,2) = cp1%centre(:,2) + dr
	if (drot == 0) return
	drot = 0
	cmid = (cc1(:,1) + cc1(:,2))/2
	c1 = cc1(:,1) - cmid			! centre location relative to mid-point
	! expand dumbbell
	r2 = dot_product(c1,c1)
	if (r2 == 0) return
	r = sqrt(r2)
	if (axis == 1) then		! rotate (y,z) about x axis
		theta = atan2(c1(3),c1(2))
		alpha = acos(c1(1)/r)
		r = r*sin(alpha)
!		write(logmsg,*) 'c1(2): ',c1(2),r*cos(theta)
		c1(2) = r*cos(theta + drot)
		c1(3) = r*sin(theta + drot)
	elseif (axis == 2) then		! rotate (x,z) about y axis
		theta = atan2(c1(1),c1(3))
		alpha = acos(c1(2)/r)
		r = r*sin(alpha)
!		write(logmsg,*) 'c1(3): ',c1(3),r*cos(theta)
		c1(3) = r*cos(theta + drot)
		c1(1) = r*sin(theta + drot)
	elseif (axis == 3) then		! rotate (x,y) about z axis
		theta = atan2(c1(2),c1(1))
		alpha = acos(c1(3)/r)
		r = r*sin(alpha)
!		write(logmsg,*) 'c1(1): ',c1(1),r*cos(theta)
		c1(1) = r*cos(theta + drot)
		c1(2) = r*sin(theta + drot)
!		dx = c1(1) - r*cos(theta)
	endif
!	call logger(logmsg)
	cc1(:,1) = cmid + c1
	cc1(:,2) = cmid - c1
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Whether or not a neighbour sphere is now in contact depends on whether or not it was
! previously in contact.
! Contact persists until d > k_detach*Rsum.
! New contact occurs when d < Rsum
!-----------------------------------------------------------------------------------------
subroutine set_neighbour_contacts(cp1)
type(cell_type), pointer :: cp1
type(cell_type), pointer :: cp2
integer :: isphere1, isphere2, nspheres1, nspheres2, k, knbr
real(REAL_KIND) :: R1, R2, c1(3), c2(3), v(3), Rsum, d
logical :: incontact

if (cp1%Iphase) then
	nspheres1 = 1
else
	nspheres1 = 2
endif
do k = 1,cp1%nbrs
	knbr = cp1%nbrlist(k)%indx
	cp2 => cell_list(knbr)
	if (cp2%Iphase) then
		nspheres2 = 1
	else
		nspheres2 = 2
	endif
	do isphere1 = 1,nspheres1
		R1 = cp1%radius(isphere1)
		c1 = cp1%centre(:,isphere1)
		do isphere2 = 1,nspheres2
			R2 = cp2%radius(isphere2)
			c2 = cp2%centre(:,isphere2)
			incontact = cp1%nbrlist(k)%contact(isphere1,isphere2)
			v = c2 - c1
			d = sqrt(dot_product(v,v))
			Rsum = R1+R2
			if (d > k_detach*Rsum) then
				incontact = .false.
			elseif (.not. incontact .and. d < Rsum) then
				incontact = .true.
			endif
			cp1%nbrlist(k)%contact(isphere1,isphere2) = incontact
		enddo
	enddo
enddo
end subroutine

end module