subroutine getnorm_pa(x,ndf,nn,node_local,nn_local,norm)
! calculate the norm of x parrallely 
use mpi_variables
implicit none
include 'mpif.h'
integer ndf
integer nn
real(8) x(ndf*nn)
integer node_local(nn_local)
integer nn_local
real(8) norm

real(8) tmp
integer i
integer j
integer node

norm=0.0d0
tmp=0.0d0
do i=1,nn_local
	node=node_local(i)
	do j=1,ndf
		tmp=tmp+x(j+(node-1)*ndf)**2
	enddo
enddo

call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(tmp,norm,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

end subroutine getnorm_pa