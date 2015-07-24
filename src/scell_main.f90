program main
use scell
use global
implicit none
integer :: ncpu, res, summarydata(100)
character*(128) :: infile, outfile, runfile
integer :: status, nlen, cnt, i, inbuflen, outbuflen
integer :: jstep, hour, irun, nsumm_interval
character*(128) :: b, c, progname
real :: vasc
logical :: update_vtk
real(REAL_KIND) :: centre(3)
real(8) :: t1, t2
integer count_0, count_1, count_rate, count_max


runfile = 'scell_main.out'
call disableTCP

outfile = 'scell_main.res'

call get_command (b, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
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
if (cnt < 2) then
    write(*,*) 'Use: ',trim(progname),' num_cpu input_file'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu 
        read(c(1:nlen),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:nlen)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        outfile = c(1:nlen)																! --> outfile
        write(*,*) 'Output file: ',outfile
    endif
end do

	inbuflen = len(infile)
	outbuflen = len(outfile)
	write(*,*) 'call execute'
	call execute(ncpu,infile,inbuflen,outfile,outbuflen,centre)
!	call cpu_time(t1)
!	t1 = wtime()
    call system_clock(count_0, count_rate, count_max)
    t1 = count_0*1.0/count_rate
	write(*,*) 'did execute: nsteps, DELTA_T: ',nsteps, DELTA_T
	do jstep = 1,Nsteps
		call simulate_step(res)
		if (res /= 0) stop
	enddo
	call terminate_run(res)
!	call cpu_time(t2)
!	t2 = wtime()
    call system_clock(count_1, count_rate, count_max)
    t2 = count_1*1.0/count_rate
	write(*,*) 'time: ',t2-t1
end 