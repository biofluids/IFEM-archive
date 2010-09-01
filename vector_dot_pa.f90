subroutine vector_dot_pa(x,y,ndf,nn,nn_local,node_local,z)
! calculate vector dot z=x*y 
use mpi_variables
implicit none
include 'mpif.h'

real(8) x(ndf,nn)
real(8) y(ndf,nn)
integer ndf
integer nn
real(8) z
integer nn_local
integer node_local(nn_local)

integer i
real(8) tmp
integer j
integer node

z=0.0d0
tmp=0.0d0
do i=1,nn_local
	node=node_local(i)
	do j=1,ndf
		tmp=tmp+x(j,node)*y(j,node)
	end do
end do
	   z=tmp
!           call mpi_barrier(mpi_comm_world,ierror)
!           call mpi_reduce(tmp,z,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror) 
!           call mpi_bcast(z,1,mpi_double_precision,0,mpi_comm_world,ierror)
end subroutine vector_dot_pa
