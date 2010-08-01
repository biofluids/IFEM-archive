!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   main.fcm                                                             c
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   initialization for the fluid and general variables is done in main.f
!   initialization for the solid variables (and input file reading) in hypo.f

program main
use parseinput
use mpi_variables
  implicit none
  include 'mpif.h'
logical DebugWait

! initialize mpi
 	call mpi_init(ierror)
	call mpi_comm_size(mpi_comm_world,ncpus,ierror)
	call mpi_comm_rank(mpi_comm_world,myid,ierror)



DebugWait=.false.
!DebugWait=.true.
        do 10 while(DebugWait)
       write(*,*) 'I am waiting'
10 end do






 !...set standard values 
  call initialize  ! for fluids

 !...read configuration files
  call parseinput_fluid  !...reading input_fluid.in
  call parseinput_solid  !...reading input_solid.in
  call nondimension
 !...echos input data
  call echoinput
  

 !...switch to main routine  
  call hypo
call mpi_finalize(ierror)
end program main
