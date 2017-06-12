! Main parameter fitting program
! Calls scell.dll subroutines 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: This has changes in subroutine run() because scell has different arguments to execute()
!       vspheroid has the same changes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fitter_mod
use scell
use strings
implicit none

integer, parameter :: nfexpt = 21, nfparam = 22, nfitlog = 23
integer, parameter :: max_params = 10, max_vals = 20

type parameter_type
	character*(48) :: ID
	integer :: nvals
	real(REAL_KIND) :: minval, maxval
	real(REAL_KIND) :: val(0:max_vals)
	integer :: ival
end type	

type experiment_type
	integer :: nvars, ntimes
	character*(24), allocatable :: varID(:)		! this must agree with size in monolayer > transfer.f90 > subroutine get_values
	integer, allocatable :: varnum(:)
	real(REAL_KIND), allocatable :: t(:)
	real(REAL_KIND), allocatable :: y(:,:)
	real(REAL_KIND), allocatable :: ymax(:)
	logical, allocatable :: is_dose(:,:)
	integer, allocatable :: itsim(:)
	real(REAL_KIND), allocatable :: ysim(:,:)
	real(REAL_KIND), allocatable :: ysim_opt(:,:)
end type

integer :: ncpu, Niters, Nparams, Nsims, Nexpts
integer :: ivalmin(0:max_params),n(0:max_params)
real(REAL_KIND) :: val(0:max_params,0:max_vals)
character*(128) :: template_infile, infile, outfile, logfile, paramfile, exptfile

type(experiment_type), allocatable, target :: experiment(:)
type(parameter_type), target :: param(0:max_params)
real(REAL_KIND) :: weight(10)	! we assume nvars <= 10

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
function exp_var_dose(iexp,str) result(dose)
integer :: iexp
character*(48) :: str
real(REAL_KIND) :: dose
character*(48) :: expstr
integer :: k, kvar
logical :: is_drug

is_drug = .false.
expstr = str
if (str /= 'RADIATION') then
!	expstr = trim(str)//'_EC'
	expstr = 'DRUG_A_EC'
	is_drug = .true.
endif
kvar = 0
dose = 0
do k = 1,experiment(iexp)%nvars
	if (expstr == experiment(iexp)%varID(k)) then
		kvar = k
		exit
	endif
enddo
if (kvar == 0) return
do k = 1,experiment(iexp)%ntimes
	if (experiment(iexp)%y(k,kvar) > 0) then
		dose = experiment(iexp)%y(k,kvar)
		if (is_drug) then
			experiment(iexp)%is_dose(k,kvar) = .true.
		endif
		exit
	endif
enddo
end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
function get_param_num(str) result(kpar)
character*(48) :: str
integer :: kpar

do kpar = 0,Nparams-1
	if (str == param(kpar)%ID) then
		return
	endif
enddo
kpar = -1
end function
		
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine read_parameters
integer :: iparam

open(nfparam,file=paramfile,status='old')
read(nfparam,*) Nparams
do iparam = 0,Nparams-1
	read(nfparam,'(a)') param(iparam)%ID
	read(nfparam,*) param(iparam)%nvals
	read(nfparam,*) param(iparam)%minval
	read(nfparam,*) param(iparam)%maxval
enddo
close(nfparam)
end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine set_param_vals
integer :: iparam, k
real(REAL_KIND) :: dval

do iparam = 0,Nparams-1
	if (param(iparam)%nvals == 1) then
		param(iparam)%val(0) = param(iparam)%minval
	else
		dval = (param(iparam)%maxval - param(iparam)%minval)/(param(iparam)%nvals-1)
		do k = 0,param(iparam)%nvals-1
			param(iparam)%val(k) = param(iparam)%minval + k*dval
		enddo
	endif
enddo
end subroutine

!----------------------------------------------------------------------------
! NOT USED
!----------------------------------------------------------------------------
subroutine read_input(nfin,old_infile,new_infile)
integer :: nfin
character*(48) :: old_infile, new_infile
character*(128) :: line
character*(50) :: str50
character*(48) :: arg(16)
integer :: nargs, nfnew = 20

open(nfin,file=old_infile,status='old')
do
	read(nfin,'(a)',end=99) line
	str50 = line(1:50)
	call parse(str50,' ',arg,nargs)
	if (nargs == 2) then
		write(*,'(2a,4x,i3)') arg(1:nargs),len(arg(2))
	else
		write(*,'(a)') str50
	endif
enddo
99 continue
close(nfin)
end subroutine

!----------------------------------------------------------------------------
! Protocol:
!  DRUG dose comes at nprot=7, RADIATION at nprot=3
!----------------------------------------------------------------------------
subroutine create_input(iexp,nfin,old_infile,new_infile)
integer :: iexp, nfin
character*(48) :: old_infile, new_infile
character*(128) :: line
character*(50) :: str50
character*(12) :: str12
character*(48) :: arg(16)
character*(48) :: paramstr, event, drugstr, drugname
integer :: nargs, nprot, kpar, nfnew = 20
real(REAL_KIND) :: dose, pval
logical :: protocol, first_rad_dose, first_drug_dose, normal

open(nfin,file=old_infile,status='old')
open(nfnew,file=new_infile,status='replace')
protocol = .false.
first_rad_dose = .true.
first_drug_dose = .true.
nprot = 0
event = ' '
normal = .true.
do
	read(nfin,'(a)',end=99) line
	str50 = line(1:50)
	call parse(str50,' ',arg,nargs)
	if (nargs == 1) then
		protocol = .true.
		paramstr = arg(1)
	else
		paramstr = arg(2)
	endif
	! Now see if the paramstr is in the parameter list.  
	! In the case of protocol:
	! The parameter ID will be of the form "DRUG_A_EC", "DRUG_A_METAB1_EC" etc
	! In fact the protocol drug name will be something like "SN30000", "PR104A" etc
	! We need to be able to associate the drug name with the drug dose line in the expt data
	! If we assume that only a single drug will ever be used in these experiments, then
	! we just need to look for "DRUG_A_EC"
	if (protocol) then
		if (first_rad_dose .and. paramstr == 'RADIATION') then
			normal = .false.
			dose = exp_var_dose(iexp,paramstr)
			if (dose > 0) then
				nprot = 0
				event = 'RADIATION'
			else
				event = ' '
				normal = .true.
			endif
		endif
		if (first_drug_dose .and. paramstr == 'DRUG') then
			normal = .false.
			nprot = 0
			event = 'DRUG'
		endif
		if (event /= ' ') nprot = nprot + 1
		if (first_drug_dose .and. event == 'DRUG' .and. nprot == 2) then		! paramstr is the drug name
			dose = exp_var_dose(iexp,paramstr)
			if (dose <= 0) then
				event = ' '
				normal = .true.
			else
				drugstr = paramstr
			endif
		endif
		if (first_rad_dose .and. event == 'RADIATION' .and. nprot == 3) then
			!dose = radiation dose for experiment(iexp)
			write(paramstr,'(f9.3)') dose
			paramstr = adjustl(paramstr)
		endif
		if (first_drug_dose .and. event == 'DRUG' .and. nprot == 8) then
			!dose = drug dose for experiment(iexp)
			write(paramstr,'(f9.5)') dose
			paramstr = adjustl(paramstr)
		endif
		if (normal) then
			write(nfnew,'(a)') trim(str50)
		else
			if (first_rad_dose .and. event == 'RADIATION') then
				write(nfnew,'(a)') trim(adjustl(paramstr))
				if (nprot == 3) then
					first_rad_dose = .false.
					normal = .true.
					event = ' '
					write(*,*) 'did rad, normal'
				endif
			elseif (first_drug_dose .and. event == 'DRUG') then
				write(nfnew,'(a)') trim(adjustl(paramstr))
				if (nprot == 8) then
					first_drug_dose = .false.
					normal = .true.
					event = ' '
					write(*,*) 'did drug, normal'
				endif
			endif
		endif
	else
		kpar = get_param_num(paramstr)
		if (kpar >= 0) then
			pval = param(kpar)%val(param(kpar)%ival)
			write(str12,'(e12.3)') pval
			str50 = adjustl(str12) // trim(arg(2))
		endif
		write(nfnew,'(a)') adjustl(str50)
!		write(nfitlog,'(a,a)') 'parameter set: ',adjustl(str50)
	endif
enddo
99 continue
close(nfin)
close(nfnew)
!stop
end subroutine

!----------------------------------------------------------------------------
! Note that a variable could be in the treatment protocol, i.e. it could
! be a drug (parent) (or even radiation dose) in which case special treatment
! will be needed.
! Note: 
!	multiple spaces are treated as a single space
!	single space is treated as a comma ','
!	space + comma ' ,' is treated as a comma ','
!	fields are delimitted by ','
!	missing value ',,' is parsed to return the string ' ', which is given the floating point value -1
!----------------------------------------------------------------------------
subroutine read_expt(nf,filename)
integer :: nf
character*(48) :: filename
character*(128) :: line
character*(48) :: arg(16)
integer :: nargs, nvarsmax
integer :: iexp, it, i, k
real(REAL_KIND) :: val
type(experiment_type), pointer :: p

open(nf,file=filename,status='old')

read(nf,'(a)',end=99) line
call parse(line,', ',arg,nargs)
read(arg(1),'(i)') Nexpts
allocate(experiment(Nexpts))
nvarsmax = 0
do iexp = 1,Nexpts
	p => experiment(iexp)
	p%nvars = nargs - 2
	allocate(p%varID(p%nvars))
	p%ntimes = 0
	do i = 1,p%nvars
		p%varID(i) = arg(i+2)
	enddo
	nvarsmax = max(nvarsmax,p%nvars)
enddo
read(nf,*) weight(1:nvarsmax)
do
	read(nf,'(a)',end=99) line
	if (len(trim(line)) == 0) cycle
	call parse(line,', ',arg,nargs)
	read(arg(1),'(i)') iexp
	p => experiment(iexp)
	p%ntimes = p%ntimes + 1
enddo
99 continue
close(nf)

do iexp = 1,Nexpts
	p => experiment(iexp)
	allocate(p%t(p%ntimes))
	allocate(p%y(p%ntimes,p%nvars))
	allocate(p%ymax(p%nvars))
	allocate(p%is_dose(p%ntimes,p%nvars))
	allocate(p%itsim(p%ntimes))
	allocate(p%ysim(p%ntimes,p%nvars))
	allocate(p%ysim_opt(p%ntimes,p%nvars))
	p%is_dose = .false.
enddo

open(nf,file=filename,status='old')
read(nf,*)	! varIDs
read(nf,*)	! weights
it = 0
iexp = 0
do
	read(nf,'(a)',end=199) line
	if (len(trim(line)) == 0) cycle
!	write(*,*) trim(line)
	call parse(line,', ',arg,nargs)
!	write(*,*) 'nargs: ',nargs
	read(arg(1),'(i)') i
	if (i /= iexp) then
		it = 0
	endif
	iexp = i
	p => experiment(iexp)
	it = it+1
	read(arg(2),*) p%t(it)
!	write(*,'(8(a,1x))') (trim(arg(k)),k=1,nargs)
	do i = 1,p%nvars
		if (arg(i+2) == ' ') then
			val = -1
		else
			read(arg(i+2),*,end=199) val
		endif
		p%y(it,i) = val
	enddo
enddo
199 continue
close(nf)

do iexp = 1,Nexpts
	p => experiment(iexp)
	do i = 1,p%nvars
		p%ymax(i) = 0
		do it = 1,p%ntimes
			p%ymax(i) = max(p%ymax(i),p%y(it,i))
		enddo
		write(nfitlog,*) 'iexp,ivar,ymax: ',iexp,i,p%ymax(i)
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------
! The crucial step will be to evaluate the experimental data variables at
! the specified time points.
!----------------------------------------------------------------------------
subroutine run(iexp)
integer :: iexp
integer :: ncpu, inbuflen, outbuflen, res,i_hypoxia_cutoff, i_growth_cutoff, jstep, i, it, k, nsumm_interval
real(REAL_KIND) :: summarydata(100)
real(REAL_KIND) :: t1, t2
real(REAL_KIND) :: centre(3)	!!!!!!!!!!!!!!!!!
type(experiment_type), pointer :: pexp

pexp => experiment(iexp)
ncpu = 1
i_hypoxia_cutoff = 3
i_growth_cutoff = 1
inbuflen = len(infile)
outbuflen = len(outfile)
write(*,*) 'call execute'
write(*,*) 'infile: ',infile
write(*,*) 'outfile: ',outfile
res = 0	!!!!!!!!!!!!
call execute(ncpu,infile,inbuflen,outfile,outbuflen,centre)		! res)	!!!!!!!!!!!!!!
write(*,*) 'did execute: nsteps, DELTA_T: ',nsteps, DELTA_T
if (res /= 0) then
	write(*,*) 'Error: execute: ',res
	stop
endif
! Set up sampling time steps, nearest DELTA_T
do i = 1,pexp%ntimes
	it = (60*60*pexp%t(i))/DELTA_T + 0.5
	pexp%itsim(i) = it		!max(it,1)
enddo
! Run the simulation
t1 = wtime()
nsumm_interval = (60*60)/DELTA_T   ! number of time steps per hour
!write(*,*) 'nsumm_interval: ',nsumm_interval
do jstep = 0,Nsteps-1
	do i = 1,pexp%ntimes
		if (pexp%itsim(i) == jstep) then
			call get_values(pexp%nvars,pexp%varID(:),pexp%ysim(i,:))
!			write(nfitlog,'(a,2i4,f8.2)') 'time: ',i,jstep,pexp%itsim(i)*DELTA_T/(60*60)
!			do k = 1,pexp%nvars
!				write(nfitlog,'(i3,a32,e12.3)') k,pexp%varID(k),pexp%ysim(i,k)
!			enddo
		endif
	enddo
	call simulate_step(res)
	if (mod(jstep,nsumm_interval) == 0) then
		call get_summary(summarydata,i_hypoxia_cutoff,i_growth_cutoff)
	endif
	if (res /= 0) then
		write(*,*) 'Error exit: ',res
		stop
	endif
enddo

!if (simulate_colony) then
!	call make_colony_distribution(colony_days,dist,ddist,ndist)
!	do idist = 1,ndist
!		write(nfrun,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
!	enddo
!endif
call terminate_run(res)
t2 = wtime()
write(*,*) 'time: ',t2-t1
end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine quantify(iexp,dfobj)
integer :: iexp
real(REAL_KIND) :: dfobj
real(REAL_KIND) :: fscale, dv, dv2
integer :: it, ivar
type(experiment_type), pointer :: pexp
logical :: use_abs = .true.
logical :: use_ysim = .true.

pexp => experiment(iexp)

write(nfitlog,*) 'quantify: iexp: ',iexp
dfobj = 0
do it = 1,pexp%ntimes
	do ivar = 1,pexp%nvars
		if (weight(ivar) == 0) cycle
		if (pexp%y(it,ivar) == -1 .or. pexp%ysim(it,ivar) == -1) cycle	! missing value
		if (pexp%is_dose(it,ivar)) cycle
		fscale = pexp%ysim(it,ivar)/pexp%ymax(ivar)
		dv = pexp%y(it,ivar) - pexp%ysim(it,ivar)
		write(nfitlog,'(2i4,4e12.3)') it,ivar,pexp%y(it,ivar),pexp%ysim(it,ivar),dv,abs(dv/pexp%y(it,ivar))
		if (use_ysim) then		
			if (pexp%ysim(it,ivar) == 0) then
				if (pexp%y(it,ivar) /= 0) then
					if (use_abs) then
						dv2 = abs(dv/pexp%y(it,ivar))
					else
						dv2 = (dv/pexp%y(it,ivar))**2
					endif
				else
					dv2 = 0
				endif
			else
				if (use_abs) then
					dv2 = abs(dv/pexp%ysim(it,ivar))
				else
					dv2 = (dv/pexp%ysim(it,ivar))**2
				endif
			endif
		else
			if (pexp%y(it,ivar) == 0) then
				if (pexp%ysim(it,ivar) /= 0) then
					if (use_abs) then
						dv2 = abs(dv/pexp%ysim(it,ivar))
					else
!						dv2 = dv*dv/pexp%ysim(it,ivar)
						dv2 = (dv/pexp%ysim(it,ivar))**2
					endif
				else
					dv2 = 0
				endif
			else
				if (use_abs) then
					dv2 = abs(dv/pexp%y(it,ivar))
				else
!					dv2 = dv*dv/pexp%y(it,ivar)
					dv2 = (dv/pexp%y(it,ivar))**2
				endif
			endif
		endif
!		if (isnan(dv2)) then
!			write(*,'(2i4,3e12.3)') it,ivar,pexp%y(it,ivar),pexp%ysim(it,ivar),dv2
!			write(*,*) 'Bad dv2'
!			write(nfitlog,*) 'Bad dv2'
!			stop
!		endif
		dfobj = dfobj + fscale*weight(ivar)*dv2
	enddo
enddo
end subroutine

end module

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
program main
use fitter_mod
implicit none
integer :: iter, iparam, isim, iexp, istart, imin, i, k, res
real(REAL_KIND) :: dfobj, fobj, fmin
logical :: ok
type(experiment_type), pointer :: pexp
character*(24), allocatable :: head(:)

integer :: status, nlen, cnt
character*(128) :: b, c, progname

infile = 'simcase.inp'
outfile = 'fitter.out'
logfile = 'fitter.log'

template_infile = 'basecase.inp'
exptfile = 'expt.dat'
paramfile = 'parameter.dat'

call get_command (b, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
endif
!write (*,*) 'command line = ', b(1:len)
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c(1:len)
progname = c(1:nlen)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 4) then
    write(*,*) 'Use: ',trim(progname),' template_file expt_file param_file Niters'
    write(*,*) '     Niters = number of iterations'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
        template_infile = c(1:nlen)
        write(*,*) 'Input file: ',template_infile
    elseif (i == 2) then
        exptfile = c(1:nlen)
        write(*,*) 'Experiment file: ',exptfile
    elseif (i == 3) then
        paramfile = c(1:nlen)
        write(*,*) 'Parameter file: ',paramfile
    elseif (i == 4) then
		read(c(1:nlen),*) Niters
		write(*,*) 'Number of iterations: ',Niters
    endif
end do

open(nfitlog,file=logfile,status='replace')
write(nfitlog,'(a,a)') 'infile: ',trim(template_infile)
write(nfitlog,'(a,a)') 'outfile: ',trim(outfile)
ncpu = 1

call read_expt(nfexpt,exptfile)
write(*,*) 'did read_expt: Nexpts: ',Nexpts
!do iexp = 1,Nexpts
!	write(*,'(8(a,1x))') 'varID: ',(trim(experiment(iexp)%varID(k)),k=1,experiment(iexp)%nvars)
!	write(*,*) 'ntimes: ',experiment(iexp)%ntimes
!	do k = 1,experiment(iexp)%ntimes
!		write(*,'(i2,6f9.2)') k,experiment(iexp)%t(k),experiment(iexp)%y(k,1:experiment(iexp)%nvars)
!	enddo
!enddo

call read_parameters

!Niters = 2
n(Nparams-1) = 1
do i = Nparams-2,0,-1
	n(i) = param(i+1)%nvals*n(i+1)
enddo
Nsims = n(0)*param(0)%nvals
write(nfitlog,*) 'Nsims: ',Nsims
write(*,*) 'Nsims: ',Nsims

call disableTCP

do iter = 1,Niters
	! This generates all the run cases
	write(*,*) 'Iteration: ',iter
	call set_param_vals
	fmin = 1.0e10
	do isim = 0,Nsims-1
		write(*,*) '    isim: ',isim
		istart = 0
		do k = 0,Nparams-1
			param(k)%ival = (isim - istart)/n(k)
			istart = istart + param(k)%ival*n(k)
		enddo
		write(nfitlog,*) isim,param(0:Nparams-1)%ival
		fobj = 0
		do iexp = 1,Nexpts
			write(nfitlog,*) '        iexp: ',iexp
			! Need to set the parameter values appropriately in the infile
			! including protocol initial values if they vary with iexp
			call create_input(iexp,nfin,template_infile, infile)
			call run(iexp)
			! Evaluate the objective function dfobj
			call quantify(iexp,dfobj)
			if (isnan(dfobj)) then
				write(*,*) 'Bad dfobj'
				write(nfitlog,*) 'Bad dfobj'
				stop
			endif
			fobj = fobj + dfobj
			pexp => experiment(iexp)
			do k = 1,pexp%ntimes
				write(nfitlog,'(i2,f8.3,10e14.6)') iexp,pexp%t(k),(pexp%ysim(k,i),i=1,pexp%nvars)
			enddo
		enddo
		write(nfitlog,*) 'fobj: ',fobj
		if (isnan(fobj)) then
			write(*,*) 'Bad fobj'
			write(nfitlog,*) 'Bad fobj'
			stop
		endif
		if (fobj < fmin) then
			fmin = fobj
			imin = isim
			ivalmin(:) = param(:)%ival
			write(nfitlog,*) 'Better fmin: isim,ivalmin,fmin: ',isim,ivalmin,fmin
			do iexp = 1,Nexpts
				pexp => experiment(iexp)
				pexp%ysim_opt = pexp%ysim
			enddo
		endif
	enddo
	write(nfitlog,*)
	write(nfitlog,*) 'Iteration: ',iter
	write(nfitlog,*) 'imin, fmin: ',imin,fmin
	write(nfitlog,*) 'ivalmin: ',ivalmin(0:Nparams-1)
	write(nfitlog,*)
	if (iter < Niters) then
		write(nfitlog,*) 'Interim optimum parameter values:'
	else
		write(nfitlog,*) 'Final optimum parameter values:'
	endif
	do k = 0,Nparams-1
		write(nfitlog,'(i2,2x,a48,e12.3)') k,param(k)%ID,param(k)%val(ivalmin(k))
	enddo
	
	pexp => experiment(1)
	allocate(head(2*pexp%nvars + 3))
!	head(1) = 'experiment'
!	head(2) = 'point'
!	head(3) = 'hour'
	do i = 1,pexp%nvars
		head(2*i-1) = trim(pexp%varID(i))//'-exp'
		head(2*i) = trim(pexp%varID(i))//'-sim'
	enddo
	write(nfitlog,*)
	write(nfitlog,'(a,10(a,1x))') 'experiment point hour ',(trim(head(i)),i=1,2*pexp%nvars)
	do iexp = 1,Nexpts
		pexp => experiment(iexp)
		do k = 1,pexp%ntimes
			write(nfitlog,'(2i4,f6.2,10e12.3)') iexp,k,pexp%t(k),(pexp%y(k,i),pexp%ysim_opt(k,i),i=1,pexp%nvars)
		enddo
	enddo
	deallocate(head)
	
	do k = 0,Nparams-1
		param(k)%minval = param(k)%val(max(0,ivalmin(k)-1))
		param(k)%maxval = param(k)%val(min(param(k)%nvals-1,ivalmin(k)+1))
	enddo
enddo   
end program

