ccccccccccccccccccc
c  mainpost.f  c
c  G. Wagner    c
ccccccccccccccccccc

      program mainpost
      include "global.h"
      integer ierr

      call MPI_INIT(ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numproc, ierr)
      call initialize
      call parseinput
      call tecplotpost
c      call tecplotpost2
      call MPI_FINALIZE(ierr)
      end
      
