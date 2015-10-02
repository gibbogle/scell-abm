module sparse_map

use global
use real_kind_mod
use omp_lib

implicit none

integer, parameter :: nfmap = 20

!integer, parameter :: NX = 33 
!integer, parameter :: NY = NX
!integer, parameter :: NZ = NX
!integer, parameter :: NXB = 35
!integer, parameter :: NYB = NXB
!integer, parameter :: NZB = NXB

integer :: nrow, nnz, nrow_b, nnz_b
integer, allocatable :: ia(:), ia_b(:), ja(:), ja_b(:), amap(:,:), amap_b(:,:)
!real(REAL_KIND), allocatable :: a(:), a_b(:)

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine read_map_file(mapfile,is_fine,ok)
character*(*) :: mapfile
logical :: is_fine, ok
integer :: k, NX_f, NY_f, NZ_f, nrow_f

write(nflog,*) 'Read map file'
open(nfmap,file=mapfile,status='old')
read(nfmap,'(3i4)') NX_f, NY_f, NZ_f
read(nfmap,'(i8)') nrow_f
if (is_fine) then
	if (NX_f /= NX .or. NY_f /= NY .or. NZ_f /= NZ .or. nrow_f /= nrow) then
		write(logmsg,*) 'Error: read_map_file: ',mapfile
		call logger(logmsg)
		write(logmsg,*) 'inconsistent parameters: NZ, NY, NZ, nrow'
		call logger(logmsg)
		write(logmsg,*) 'this file: ',NZ_f, NY_f, NZ_f, nrow_f
		call logger(logmsg)
		write(logmsg,*) 'this run:  ',NZ, NY, NZ, nrow
		call logger(logmsg)
		write(logmsg,*) 'Delete ',mapfile, ' and run the program again'
		call logger(logmsg)
		ok = .false.
		return
	endif
	read(nfmap,'(i8)') nnz
	read(nfmap,*)
	read(nfmap,'(10i10)') ia(1:nrow+1)
	read(nfmap,*)
	read(nfmap,'(10i7)') ja(1:nnz)
	read(nfmap,*)
	do k = 1,nnz
		read(nfmap,'(4i6)') amap(k,:)
	enddo
else
	if (NX_f /= NXB .or. NY_f /= NYB .or. NZ_f /= NZB .or. nrow_f /= nrow_b) then
		write(logmsg,*) 'Error: read_map_file: ',mapfile
		call logger(logmsg)
		write(logmsg,*) 'inconsistent parameters: NZ, NY, NZ, nrow'
		call logger(logmsg)
		write(logmsg,*) 'this file: ',NZ_f, NY_f, NZ_f, nrow_f
		call logger(logmsg)
		write(logmsg,*) 'this run:  ',NZB, NYB, NZB, nrow_b
		call logger(logmsg)
		write(logmsg,*) 'Delete ',mapfile, ' and run the program again'
		call logger(logmsg)
		ok = .false.
		return
	endif
	read(nfmap,'(i8)') nnz_b
	read(nfmap,*)
	read(nfmap,'(10i10)') ia_b(1:nrow_b+1)
	read(nfmap,*)
	read(nfmap,'(10i7)') ja_b(1:nnz_b)
	read(nfmap,*)
	do k = 1,nnz_b
		read(nfmap,'(4i6)') amap_b(k,:)
	enddo
endif
close(nfmap)
write(nflog,*) 'nnz_b: ',nnz_b
ok = .true.
!write(*,'(10i6)') ja_b(1:100)
!write(*,'(10i6)') amap_b(1:100,0)
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine write_map_file(mapfile,is_fine)
character*(*) :: mapfile
logical :: is_fine
integer :: k

write(nflog,*) 'Write map file'
open(nfmap,file=mapfile,status='replace')
if (is_fine) then
	write(nfmap,'(3i4)') NX, NY, NZ
	write(nfmap,'(i8)') nrow
	write(nfmap,'(i8)') nnz
	write(nfmap,'(a)') 'ia'
	write(nfmap,'(10i10)') ia(1:nrow+1)
	write(nfmap,'(a)') 'ja'
	write(nfmap,'(10i7)') ja(1:nnz)
	write(nfmap,'(a)') 'amap'
	do k = 1,nnz
		write(nfmap,'(4i6)') amap(k,0:3)
	enddo
else
	write(nfmap,'(3i4)') NXB, NYB, NZB
	write(nfmap,'(i8)') nrow_b
	write(nfmap,'(i8)') nnz_b
	write(nfmap,'(a)') 'ia_b'
	write(nfmap,'(10i10)') ia_b(1:nrow_b+1)
	write(nfmap,'(a)') 'ja_b'
	write(nfmap,'(10i7)') ja_b(1:nnz_b)
	write(nfmap,'(a)') 'amap_b'
	do k = 1,nnz_b
		write(nfmap,'(4i6)') amap_b(k,0:3)
	enddo
endif
close(nfmap)
end subroutine

!-------------------------------------------------------------------------------------------
! This assumes impermeable (reflective) boundaries on all sides.
! Does not account for OXYGEN!!!!!
!-------------------------------------------------------------------------------------------
subroutine make_sparse_map(mapfile,is_fine,ok)
character*(*) :: mapfile
logical :: is_fine, ok
integer :: ix, iy, iz, k, i, nnz_t, nrow_t, nx_t, ny_t, nz_t, krow, kcol
integer :: idx, idy, idz, ixx, iyy, izz, chk(7), knt
real(REAL_KIND) :: Kd, Kr, aval(7)
integer, allocatable :: arow(:)
logical :: file_exists = .false., upper_bdry
integer, parameter :: m = 3

inquire(file=mapfile,exist=file_exists)
if (file_exists) then
	call read_map_file(mapfile,is_fine,ok)
	return
endif
if (is_fine) then
	nx_t = NX
	ny_t = NY
	nz_t = NZ
	nrow_t = nrow
else
	nx_t = NXB
	ny_t = NYB
	nz_t = NZB
	nrow_t = nrow_b
endif
allocate(arow(nrow_t))
nnz_t = 0
do ix = 1,nx_t
!	write(*,*) 'make_sparse_map: ix: ',ix
	do iy = 1,ny_t
		do iz = 1,nz_t
			upper_bdry = .false.
			arow = 0
			krow = (ix-1)*ny_t*nz_t + (iy-1)*nz_t + iz
			arow(krow) = 2*m
			
			chk = 0
			ixx = ix - 1
			if (ixx < 1) then
				ixx = ixx + 2
			endif
			kcol = krow + (ixx-ix)*ny_t*nz_t
			arow(kcol) = arow(kcol) - 1
			aval(1) = arow(kcol)
			chk(1) = kcol
			iyy = iy - 1
			if (iyy < 1) then
				iyy = iyy + 2
			endif
			kcol = krow + (iyy-iy)*nz_t
			arow(kcol) = arow(kcol) - 1
			aval(2) = arow(kcol)
			chk(2) = kcol
			izz = iz - 1
			if (izz < 1) then
				izz = izz + 2
			endif
			kcol = krow + (izz-iz)
			arow(kcol) = arow(kcol) - 1
			aval(3) = arow(kcol)
			chk(3) = kcol
			
			chk(4) = krow
			aval(4) = arow(krow)
			
			izz = iz + 1
			if (izz > nz_t) then
				izz = izz - 2
			endif
			kcol = krow + (izz-iz)
			if (arow(kcol) /= 0) chk(3) = 0
			arow(kcol) = arow(kcol) - 1
			aval(5) = arow(kcol)
			chk(5) = kcol
			
			iyy = iy + 1
			if (iyy > ny_t) then
				upper_bdry = .true.
				iyy = iyy - 2
			endif
			kcol = krow + (iyy-iy)*nz_t
			if (arow(kcol) /= 0) chk(2) = 0
			arow(kcol) = arow(kcol) - 1
			aval(6) = arow(kcol)
			chk(6) = kcol
			
			ixx = ix + 1
			if (ixx > nx_t) then
				ixx = ixx - 2
			endif
			kcol = krow + (ixx-ix)*ny_t*nz_t
			if (arow(kcol) /= 0) chk(1) = 0
			arow(kcol) = arow(kcol) - 1
			aval(7) = arow(kcol)
			chk(7) = kcol
			
!			write(*,'(7i8)') chk
			if (is_fine) then
				ia(krow) = nnz_t + 1
				do i = 1,7
					if (chk(i) /= 0) then
						nnz_t = nnz_t + 1
						ja(nnz_t) = chk(i)
						amap(nnz_t,0) = aval(i)
						amap(nnz_t,1:3) = [ix,iy,iz]
	!					if (upper_bdry) write(*,'(7i6,f8.3)') i,ix,iy,iz,nnz_t,krow,chk(i),aval(i)
					endif
				enddo
			else
				ia_b(krow) = nnz_t + 1
				do i = 1,7
					if (chk(i) /= 0) then
						nnz_t = nnz_t + 1
						ja_b(nnz_t) = chk(i)
						amap_b(nnz_t,0) = aval(i)
						amap_b(nnz_t,1:3) = [ix,iy,iz]
					endif
				enddo
			endif
		enddo
	enddo
enddo
if (is_fine) then
	ia(nrow+1) = nnz_t+1	
else
	ia_b(nrow_t+1) = nnz_t+1	
endif
if (is_fine) then
	nnz = nnz_t
else
	nnz_b = nnz_t
	write(nflog,*) 'nnz_b: ',nnz_b
!	write(*,'(10i5)') ja_b(1:100)
!	write(*,'(10i5)') amap_b(1:100,0)
endif
deallocate(arow)
call write_map_file(mapfile,is_fine)
ok = .true.
end subroutine


!-------------------------------------------------------------------------------------------
! This version is for the embedded fine grid.
!-------------------------------------------------------------------------------------------
subroutine read_emap_file(mapfile,is_fine,ok)
character*(*) :: mapfile
logical :: is_fine, ok
integer :: k, NX_f, NY_f, NZ_f, nrow_f

write(nflog,*) 'Read emap file'
if (.not.is_fine) then
	write(logmsg,*) 'emap is only for the embedded (fine) grid'
	call logger(logmsg)
	ok = .false.
	return
endif
open(nfmap,file=mapfile,status='old')
read(nfmap,'(3i4)') NX_f, NY_f, NZ_f
read(nfmap,'(i8)') nrow_f
if (is_fine) then
	if (NX_f+2 /= NX .or. NY_f+2 /= NY .or. NZ_f+1 /= NZ .or. nrow_f /= nrow) then
		write(logmsg,*) 'Error: read_map_file: ',mapfile
		call logger(logmsg)
		write(logmsg,*) 'inconsistent parameters: NZ, NY, NZ, nrow'
		call logger(logmsg)
		write(logmsg,*) 'this file: ',NZ_f+2, NY_f+2, NZ_f+1, nrow_f
		call logger(logmsg)
		write(logmsg,*) 'this run:  ',NZ, NY, NZ, nrow
		call logger(logmsg)
		write(logmsg,*) 'Delete ',mapfile, ' and run the program again'
		call logger(logmsg)
		ok = .false.
		return
	endif
	read(nfmap,'(i8)') nnz
	read(nfmap,*)
	read(nfmap,'(10i10)') ia(1:nrow+1)
	read(nfmap,*)
	read(nfmap,'(10i7)') ja(1:nnz)
	read(nfmap,*)
	do k = 1,nnz
		read(nfmap,'(4i6)') amap(k,:)
	enddo
else
	if (NX_f /= NXB .or. NY_f /= NYB .or. NZ_f /= NZB .or. nrow_f /= nrow_b) then
		write(logmsg,*) 'Error: read_map_file: ',mapfile
		call logger(logmsg)
		write(logmsg,*) 'inconsistent parameters: NZ, NY, NZ, nrow'
		call logger(logmsg)
		write(logmsg,*) 'this file: ',NZ_f, NY_f, NZ_f, nrow_f
		call logger(logmsg)
		write(logmsg,*) 'this run:  ',NZB, NYB, NZB, nrow_b
		call logger(logmsg)
		write(logmsg,*) 'Delete ',mapfile, ' and run the program again'
		call logger(logmsg)
		ok = .false.
		return
	endif
	read(nfmap,'(i8)') nnz_b
	read(nfmap,*)
	read(nfmap,'(10i10)') ia_b(1:nrow_b+1)
	read(nfmap,*)
	read(nfmap,'(10i7)') ja_b(1:nnz_b)
	read(nfmap,*)
	do k = 1,nnz_b
		read(nfmap,'(4i6)') amap_b(k,:)
	enddo
endif
close(nfmap)
ok = .true.
!write(*,*) 'nnz, nrow_f: ',nnz,nrow_f
!write(*,*) 'ia: '
!write(*,'(20i6)') ia
!write(*,*) 'ja: '
!write(*,'(20i6)') ja
!write(*,'(20i3)') amap(:,0)
end subroutine

!-------------------------------------------------------------------------------------------
! This version is for the embedded fine grid.
!-------------------------------------------------------------------------------------------
subroutine write_emap_file(mapfile,is_fine)
character*(*) :: mapfile
logical :: is_fine
integer :: k

write(nflog,*) 'Write emap file'
if (.not.is_fine) then
	write(logmsg,*) 'emap is only for the embedded (fine) grid'
	stop
endif
open(nfmap,file=mapfile,status='replace')
if (is_fine) then
	write(nfmap,'(3i4)') NX-2, NY-2, NZ-1
	write(nfmap,'(i8)') nrow
	write(nfmap,'(i8)') nnz
	write(nfmap,'(a)') 'ia'
	write(nfmap,'(10i10)') ia(1:nrow+1)
	write(nfmap,'(a)') 'ja'
	write(nfmap,'(10i7)') ja(1:nnz)
	write(nfmap,'(a)') 'amap'
	do k = 1,nnz
		write(nfmap,'(4i6)') amap(k,0:3)
	enddo
else
	write(nfmap,'(3i4)') NXB, NYB, NZB
	write(nfmap,'(i8)') nrow_b
	write(nfmap,'(i8)') nnz_b
	write(nfmap,'(a)') 'ia_b'
	write(nfmap,'(10i10)') ia_b(1:nrow_b+1)
	write(nfmap,'(a)') 'ja_b'
	write(nfmap,'(10i7)') ja_b(1:nnz_b)
	write(nfmap,'(a)') 'amap_b'
	do k = 1,nnz_b
		write(nfmap,'(4i6)') amap_b(k,0:3)
	enddo
endif
close(nfmap)
end subroutine

!-------------------------------------------------------------------------------------------
! This version is for the embedded fine grid.
! This assumes impermeable (reflective) boundaries on the iz=1 boundary.  All other boundaries
! will use prescribed concentrations, computed on the coarse grid.
!-------------------------------------------------------------------------------------------
subroutine make_sparse_emap(mapfile,is_fine,ok)
character*(*) :: mapfile
logical :: is_fine, ok
integer :: ix, iy, iz, k, i, nnz_t, nrow_t, nx_t, ny_t, nz_t, krow, kcol
integer :: idx, idy, idz, ixx, iyy, izz, chk(7), knt, asum
real(REAL_KIND) :: Kd, Kr, aval(7)
integer, allocatable :: arow(:)
logical :: file_exists = .false., upper_bdry
integer, parameter :: m = 3

inquire(file=mapfile,exist=file_exists)
if (file_exists) then
	call read_emap_file(mapfile,is_fine,ok)
	return
endif
if (is_fine) then
	nx_t = NX-2
	ny_t = NY-2
	nz_t = NZ-1
	nrow_t = nx_t*ny_t*nz_t	! now the number of rows is reduced!
else
	write(logmsg,*) 'emap is only for the embedded (fine) grid'
	call logger(logmsg)
	ok = .false.
	return
endif
allocate(arow(nrow_t))
nnz_t = 0
do ix = 2,NX-1
!	write(*,*) 'make_sparse_emap: ix: ',ix
	do iy = 2,NY-1
		do iz = 1,NZ-1
!			upper_bdry = .false.
			arow = 0
			krow = (ix-2)*ny_t*nz_t + (iy-2)*nz_t + iz
			arow(krow) = 2*m
			chk = 0

			ixx = ix - 1
!			if (ixx < 1) then
!				ixx = ixx + 2
!			endif
!			arow(kcol) = arow(kcol) - 1
			if (ixx /= 1) then
				kcol = krow + (ixx-ix)*ny_t*nz_t
				arow(kcol) = -1
				aval(1) = arow(kcol)
				chk(1) = kcol
			endif
			
			iyy = iy - 1
!			if (iyy < 1) then
!				iyy = iyy + 2
!			endif
!			arow(kcol) = arow(kcol) - 1
			if (iyy /= 1) then
				kcol = krow + (iyy-iy)*nz_t
				arow(kcol) = -1
				aval(2) = arow(kcol)
				chk(2) = kcol
			endif
			
			izz = iz - 1
			if (izz < 1) then
				izz = izz + 2
			endif
			kcol = krow + (izz-iz)
			arow(kcol) = arow(kcol) - 1
			aval(3) = arow(kcol)
			chk(3) = kcol
			
			chk(4) = krow
			aval(4) = arow(krow)
			
			izz = iz + 1
			kcol = krow + (izz-iz)
			if (izz /= NZ) then
				kcol = krow + (izz-iz)
				arow(kcol) = -1
				aval(5) = arow(kcol)
				chk(5) = kcol
			endif
			
			iyy = iy + 1
			if (iyy /= NY) then
				kcol = krow + (iyy-iy)*nz_t
				arow(kcol) = -1
				aval(6) = arow(kcol)
				chk(6) = kcol
			endif
			
			ixx = ix + 1
			if (ixx /= NX) then
				kcol = krow + (ixx-ix)*ny_t*nz_t
				arow(kcol) = -1
				aval(7) = arow(kcol)
				chk(7) = kcol
			endif
			
!			write(*,'(7i8)') chk
			ia(krow) = nnz_t + 1
			asum = 0
			do i = 1,7
				if (chk(i) /= 0) then
					nnz_t = nnz_t + 1
					ja(nnz_t) = chk(i)
					amap(nnz_t,0) = aval(i)
					amap(nnz_t,1:3) = [ix,iy,iz]
					asum = asum + aval(i)
				endif
			enddo
!			if (asum /= 0) then
!				write(*,'(7e11.3)') aval
!				write(*,*) 'bad asum: ',asum
!				stop
!			endif
		enddo
	enddo
enddo
ia(nrow+1) = nnz_t+1
nnz = nnz_t	
write(nflog,*) 'nnz: ',nnz
!write(*,'(10i5)') ja(1:nnz)
!write(*,'(10i5)') amap(1:nnz,0)
deallocate(arow)
call write_emap_file(mapfile,is_fine)
ok = .true.
end subroutine


end module