	subroutine getnorm(v1,v2,mn,norm)
	implicit none
      include "global.h"
	integer i,mn
	real*8 v1(mn),v2(mn),mynorm,norm,sdot
	integer ierr

      mynorm = 0.0
	do i=1,mn
	 mynorm =  mynorm + v1(i)*v2(i)
	enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_REDUCE (mynorm,norm,1,MPI_DOUBLE_PRECISION,
     &                   MPI_SUM,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(norm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	return
	end
