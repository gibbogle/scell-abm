module envelope

use global
implicit none

integer, parameter :: MAXN = 5000
real(REAL_KIND), parameter :: dz = 2*Raverage	! cm

integer, allocatable :: czlist(:,:)
integer, allocatable :: nzlist(:)
real(REAL_KIND), allocatable :: area(:)
real(REAL_KIND), allocatable :: r_inner(:)
real(REAL_KIND), allocatable :: r_outer(:)

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function getCentre() result (C)
integer :: kcell, k, n
real(REAL_KIND) :: csum(3), C(3)
type(cell_type), pointer :: cp

n = 0
csum = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
	cp => cell_list(kcell)
	do k = 1,cp%nspheres
		csum = csum + cp%centre(:,k)
	enddo
enddo
C = csum/n
end function

!-----------------------------------------------------------------------------------------
! Estimate blob radius from the volume.
!-----------------------------------------------------------------------------------------
function getRadius() result(R)
real(REAL_KIND) :: R
real(REAL_KIND) :: vol_cm3, area_cm2
!real(REAL_KIND) :: cntr(3), rng(3), diam

!call getBlobCentreRange(cntr,rng, R)
call getVolume(vol_cm3, area_cm2)
R = sqrt(area_cm2/PI)
end function

!---------------------------------------------------------------------------------------
! Returns total blob volume in units of cm3
! The region occupied by the blob is divided into layers dz thick, and the cells in each
! layer are identified.
!---------------------------------------------------------------------------------------
subroutine getVolume(volume,maxarea)
real(REAL_KIND) :: volume, maxarea
real(REAL_KIND) :: z, zmin, zmax, r, cntr(2)
integer :: kcell, iz, nzz, nbig, nzlmax
type(cell_type), pointer :: cp

zmin = 1.0e10
zmax = -1.0e10
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	z = cp%centre(3,1)
	zmin = min(zmin,z)
	zmax = max(zmax,z)
enddo
zmin = zmin - dz
zmax = zmax + dz
nzz = (zmax-zmin)/dz + 1
r = (3*Ncells/(4*PI))**(1./3.)
nbig = PI*r**2
nbig = 1.5*nbig
! We need an array to store the cells indicies
allocate(nzlist(nzz))
allocate(czlist(nzz,nbig))
allocate(area(nzz))

nzlist = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	z = cp%centre(3,1)
	iz = (z-zmin)/dz + 1
	nzlist(iz) = nzlist(iz) + 1
	if (nzlist(iz) > nbig) then
		write(*,*) 'Error: getVolume: nzlist(iz) too big: ',iz,nzlist(iz),nbig
		stop
	endif
	czlist(iz,nzlist(iz)) = kcell
enddo

nzlmax = 0
do iz = 1,nzz
!	write(*,*) iz,nzlist(iz)
	nzlmax = max(nzlmax,nzlist(iz))
	if (nzlist(iz) > 0) then
		call getArea(iz,area(iz),cntr)
	else
		area(iz) = 0
	endif
enddo
!write(*,*) 'nzz,nbig,nzlmax: ',nzz,nbig,nzlmax
volume = 0
maxarea = 0
do iz = 1,nzz
	volume = volume + area(iz)
	if (area(iz) > maxarea) maxarea = area(iz)
enddo
deallocate(nzlist)
deallocate(czlist)
deallocate(area)
volume = dz*volume		! cm3
r = sqrt(maxarea/PI)	! cm
!write(*,'(a,i6,3e10.3,f6.1)') 'Ncells,volume,area,diam (cm,um): ',Ncells,volume,maxarea,2*r,2*r*10000
end subroutine

!---------------------------------------------------------------------------------------
! nsmax = dimension of rad(:,2)
!---------------------------------------------------------------------------------------
subroutine getSlices(nslices,dzslice,nsmax,rad)
integer :: nslices, nsmax
real(REAL_KIND) :: dzslice, rad(:,:)
real(REAL_KIND) :: z, zmin, zmax, r, cntr(2)
integer :: kcell, iz, nzz, nbig, i, iz1, iz2
type(cell_type), pointer :: cp

zmin = 1.0e10
zmax = -1.0e10
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	z = cp%centre(3,1)
	zmin = min(zmin,z)
	zmax = max(zmax,z)
enddo
zmin = zmin - dz
zmax = zmax + dz
nzz = (zmax-zmin)/dz + 1
r = (3*Ncells/(4*PI))**(1./3.)
nbig = PI*r**2
nbig = 1.5*nbig
! We need an array to store the cells indicies
allocate(nzlist(nzz))
allocate(czlist(nzz,nbig))
allocate(area(nzz))
allocate(r_inner(nzz))
allocate(r_outer(nzz))

nzlist = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	z = cp%centre(3,1)
	iz = (z-zmin)/dz + 1
	nzlist(iz) = nzlist(iz) + 1
	if (nzlist(iz) > nbig) then
		write(*,*) 'Error: getSlices: nzlist(iz) too big: ',iz,nzlist(iz),nbig
		stop
	endif
	czlist(iz,nzlist(iz)) = kcell
enddo

do iz = 1,nzz
	if (nzlist(iz) > 0) then
		call getArea(iz,area(iz),cntr)
        r_outer(iz) = sqrt(area(iz)/PI)	! cm
        call getInnerRadius(iz,cntr,r_inner(iz))
	else
		area(iz) = 0
		r_outer(iz) = 0
		r_inner(iz) = 0
	endif
enddo
iz1 = 0
do iz = 1,nzz
    if (r_outer(iz) > 0) then
        if (iz1 == 0) iz1 = iz
        iz2 = iz
    endif
enddo
nslices = iz2 - iz1 + 1
nslices = min(nslices,nsmax)
do i = 1,nslices
    iz = iz1 + i - 1
    rad(i,1) = r_inner(iz)
    rad(i,2) = r_outer(iz)
enddo
dzslice = dz
deallocate(nzlist)
deallocate(czlist)
deallocate(area)
deallocate(r_inner)
deallocate(r_outer)

end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine getInnerRadius(iz,cntr,r_in)
integer :: iz
real(REAL_KIND) :: r_in, cntr(2)
integer :: i, kcell
real(REAL_KIND) :: r2, r2min, dx, dy
type(cell_type), pointer :: cp

r2min = 1.0e10
do i = 1,nzlist(iz)
	kcell = czlist(iz,i)
	cp => cell_list(kcell)
	dx = cp%centre(1,1) - cntr(1)
	dy = cp%centre(2,1) - cntr(2)
	r2 = dx**2 + dy**2
	r2min = min(r2,r2min)
enddo
r_in = sqrt(r2min)
end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine getArea(iz, area, cntr)
integer :: iz
real(REAL_KIND) :: area, cntr(2)
real(REAL_KIND), allocatable :: x(:), y(:), xv(:), yv(:)
real(REAL_KIND) :: ave(2), xc, yc, d2, d2max, dx(4), dy(4), dx_um
real(REAL_KIND) :: delta = 0.20*dz
integer, allocatable :: iwk(:), vertex(:)
integer :: kcell, n, nvert, i, i1, i2, k, kmax
type(cell_type), pointer :: cp

allocate(x(MAXN))
allocate(y(MAXN))
allocate(iwk(MAXN))
allocate(vertex(1000))
allocate(xv(1000))
allocate(yv(1000))

n = nzlist(iz)
do i = 1,n
	kcell = czlist(iz,i)
!	write(*,*) 'iz,i,kcell: ',iz,i,kcell
	cp => cell_list(kcell)
	x(i) = cp%centre(1,1)
	y(i) = cp%centre(2,1)
enddo
if (n <= 3) then
	area = n*dz*dz
	return
endif
call enveloper(x,y,n,vertex,nvert,iwk)
!write(*,*) 'get_diameter: nvert: ',nvert
ave = 0
do i = 1,nvert
	xv(i) = x(vertex(i))
	yv(i) = y(vertex(i))
!	write(*,'(i6,2f6.1)') i,xv(i),yv(i)
	ave(1) = ave(1) + xv(i)
	ave(2) = ave(2) + yv(i)
enddo
cntr = ave/nvert	! this is the approximate centre
! Now adjust (xv,yv) to the site corner most remote from the centre
dx = [-delta, delta, delta, -delta]
dy = [-delta, -delta, delta, delta]
do i = 1,nvert
	d2max = 0
	do k = 1,4
		xc = xv(i) + dx(k)
		yc = yv(i) + dy(k)
		d2 = (xc - cntr(1))**2 + (yc-cntr(2))**2
		if (d2 > d2max) then
			d2max = d2
			kmax = k
		endif
	enddo
	xv(i) = xv(i) + dx(kmax)
	yv(i) = yv(i) + dy(kmax)
!	write(*,'(i6,2f6.1)') i,xv(i),yv(i)
enddo
area = 0
do i1 = 1,nvert
	i2 = i1+1
	if (i2 > nvert) i2 = 1
	area = area + xv(i1)*yv(i2) - xv(i2)*yv(i1)
enddo
area = -area/2
!dx_um = DELTA_X*10000	! site size in um
!area = area*dx_um*dx_um	! area in um^2
!diam = 2*sqrt(area/PI)
!write(*,'(a,i4,2f8.1)') 'area: ',iz,area,diam

deallocate(xv)
deallocate(yv)
deallocate(x)
deallocate(y)
deallocate(iwk)
deallocate(vertex)

end subroutine

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
SUBROUTINE enveloper(x, y, n, vertex, nvert, iwk)

!  Find the vertices (in clockwise order) of a polygon enclosing
!  the points (x(i), y(i), i=1, ..., n.

!  On output, vertex(i), i=1, ..., nvert contains the numbers of the vertices.
!  iwk() is an integer work array which must have dimension at least n
!  in the calling program.

!  There is a limit of 1000 vertices imposed by the dimension of array next.

!  Programmer: Alan Miller
!  Latest revision - 12 September 1987
!  Fortran 90 version - 8 August 1996

INTEGER :: n, vertex(n), nvert, iwk(n)
real(REAL_KIND)    :: x(n), y(n)

!       Local variables

INTEGER :: next(1000), i, i1, i2, j, jp1, jp2, i2save, i3, i2next
real(REAL_KIND)    :: xmax, xmin, ymax, ymin, dist, dmax, dmin, x1, y1, dx, dy, x2, y2, &
           dx1, dx2, dmax1, dmax2, dy1, dy2, temp, zero = 0.0

IF (n < 2) RETURN

!  Choose the points with smallest & largest x- values as the
!  first two vertices of the polygon.

IF (x(1) > x(n)) THEN
  vertex(1) = n
  vertex(2) = 1
  xmin = x(n)
  xmax = x(1)
ELSE
  vertex(1) = 1
  vertex(2) = n
  xmin = x(1)
  xmax = x(n)
END IF

DO i = 2, n-1
  temp = x(i)
  IF (temp < xmin) THEN
    vertex(1) = i
    xmin = temp
  ELSE IF (temp > xmax) THEN
    vertex(2) = i
    xmax = temp
  END IF
END DO

!       Special case, xmax = xmin.

IF (xmax == xmin) THEN
  IF (y(1) > y(n)) THEN
    vertex(1) = n
    vertex(2) = 1
    ymin = y(n)
    ymax = y(1)
  ELSE
    vertex(1) = 1
    vertex(2) = n
    ymin = y(1)
    ymax = y(n)
  END IF

  DO i = 2, n-1
    temp = y(i)
    IF (temp < ymin) THEN
      vertex(1) = i
      ymin = temp
    ELSE IF (temp > ymax) THEN
      vertex(2) = i
      ymax = temp
    END IF
  END DO

  nvert = 2
  IF (ymax == ymin) nvert = 1
  RETURN
END IF

!  Set up two initial lists of points; those points above & those below the
!  line joining the first two vertices.    next(i) will hold the pointer to the
!  point furthest from the line joining vertex(i) to vertex(i+1) on the left
!  hand side.

i1 = vertex(1)
i2 = vertex(2)
iwk(i1) = -1
iwk(i2) = -1
dx = xmax - xmin
y1 = y(i1)
dy = y(i2) - y1
dmax = zero
dmin = zero
next(1) = -1
next(2) = -1

DO i = 1, n
  IF (i == vertex(1) .OR. i == vertex(2)) CYCLE
  dist = (y(i) - y1)*dx - (x(i) - xmin)*dy
  IF (dist > zero) THEN
    iwk(i1) = i
    i1 = i
    IF (dist > dmax) THEN
      next(1) = i
      dmax = dist
    END IF
  ELSE IF (dist < zero) THEN
    iwk(i2) = i
    i2 = i
    IF (dist < dmin) THEN
      next(2) = i
      dmin = dist
    END IF
  END IF
END DO

!  Ends of lists are indicated by pointers to -ve positions.

iwk(i1) = -1
iwk(i2) = -1
nvert = 2

j = 1

!  Start of main process.

!  Introduce new vertex between vertices j & j+1, if one has been found.
!  Otherwise increase j.   Exit if no more vertices.

40 IF (next(j) < 0) THEN
IF (j == nvert) RETURN
j = j + 1
GO TO 40
END IF

jp1 = j + 1
DO i = nvert, jp1, -1
  vertex(i+1) = vertex(i)
  next(i+1) = next(i)
END DO
jp2 = jp1 + 1
nvert = nvert + 1
IF (jp2 > nvert) jp2 = 1
i1 = vertex(j)
i2 = next(j)
i3 = vertex(jp2)
vertex(jp1) = i2

!  Process the list of points associated with vertex j.   New list at vertex j
!  consists of those points to the left of the line joining it to the new
!  vertex (j+1).   Similarly for the list at the new vertex.
!  Points on or to the right of these lines are dropped.

x1 = x(i1)
x2 = x(i2)
y1 = y(i1)
y2 = y(i2)
dx1 = x2 - x1
dx2 = x(i3) - x2
dy1 = y2 - y1
dy2 = y(i3) - y2
DMAX1 = zero
dmax2 = zero
next(j) = -1
next(jp1) = -1
i2save = i2
i2next = iwk(i2)
i = iwk(i1)
iwk(i1) = -1
iwk(i2) = -1

60 IF (i /= i2save) THEN
  dist = (y(i) - y1)*dx1 - (x(i) - x1)*dy1
  IF (dist > zero) THEN
    iwk(i1) = i
    i1 = i
    IF (dist > DMAX1) THEN
      next(j) = i
      DMAX1 = dist
    END IF
  ELSE
    dist = (y(i) - y2)*dx2 - (x(i) - x2)*dy2
    IF (dist > zero) THEN
      iwk(i2) = i
      i2 = i
      IF (dist > dmax2) THEN
        next(jp1) = i
        dmax2 = dist
      END IF
    END IF
  END IF
  i = iwk(i)
ELSE
  i = i2next
END IF

!  Get next point from old list at vertex j.

IF (i > 0) GO TO 60

!  End lists with -ve values.

iwk(i1) = -1
iwk(i2) = -1

GO TO 40
END SUBROUTINE enveloper

end module
