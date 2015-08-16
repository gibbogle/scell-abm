interface

	subroutine itsol_create_matrix(ic, nrow, nnz, a, ja, ia, ierr) bind(C,name='itsol_create_matrix')
	use :: iso_c_binding
	integer(c_int),value :: ic, nrow, nnz
	real(c_double) :: a(*)
	integer(c_int) :: ja(*), ia(*)
	integer(c_int) :: ierr
	end subroutine

end interface

interface

	subroutine itsol_create_precond_ILUK(ic, nfill_level, ierr) bind(C,name='itsol_create_precond_ILUK')
	use :: iso_c_binding
	integer(c_int),value :: ic, nfill_level
	integer(c_int) :: ierr
	end subroutine

end interface

interface

	subroutine itsol_create_precond_ILUT(ic, nfill_level, tol, ierr) bind(C,name='itsol_create_precond_ILUT')
	use :: iso_c_binding
	integer(c_int),value :: ic, nfill_level
	real(c_double),value :: tol
	integer(c_int) :: ierr
	end subroutine

end interface

interface

	subroutine itsol_create_precond_VBILUK(ic, nfill_level, ierr) bind(C,name='itsol_create_precond_VBILUK')
	use :: iso_c_binding
	integer(c_int),value :: ic, nfill_level
	integer(c_int) :: ierr
	end subroutine

end interface

interface

	subroutine itsol_create_precond_ARMS(ic, nfill_level, tol, ierr) bind(C,name='itsol_create_precond_ARMS')
	use :: iso_c_binding
	integer(c_int),value :: ic, nfill_level
	real(c_double),value :: tol
	integer(c_int) :: ierr
	end subroutine

end interface

interface

	subroutine itsol_free_precond_ILU(ic, ierr) bind(C,name='itsol_free_precond_ILU')
	use :: iso_c_binding
	integer(c_int),value :: ic
	integer(c_int) :: ierr
	end subroutine

end interface

interface

	subroutine itsol_solve_fgmr_ILU(ic, rhs, x, im_krylov, maxits, tol, iters, ierr) bind(C,name='itsol_solve_fgmr_ILU')
	use :: iso_c_binding
	real(c_double) :: rhs(*), x(*)
	integer(c_int),value :: ic, im_krylov, maxits
	real(c_double),value :: tol
	integer(c_int) :: iters, ierr
	end subroutine

end interface

interface

	subroutine itsol_free_matrix(ic, ierr) bind(C,name='itsol_free_matrix')
	use :: iso_c_binding
	integer(c_int),value :: ic
	integer(c_int) :: ierr
	end subroutine

end interface
