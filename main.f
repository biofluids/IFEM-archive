c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	main.fcm                                                             c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  program main
      include "global.h"
	  integer ierr
	  
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numproc, ierr)
	  call initialize
	  call parseinput
	  call nondimension
	  call echoinput
      call hypo
      call MPI_FINALIZE( ierr )
	  
	  end
	  
















