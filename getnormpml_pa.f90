subroutine getnormpml_pa(x,ndf,nn,node_local,nn_local,norm)
! calculate the norm of x parrallely
use mpi_variables
use pml_variables
implicit none
include 'mpif.h'
integer ndf
integer nn
real(8) x(ndf*(nn+nn_PML))
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
	if (seqcPML(node) > 0) then
	    do j=1,ndf
	        tmp=tmp+x(j+(nn+seqcPML(node)-1)*ndf)**2
	    enddo
	endif
enddo

call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(tmp,norm,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

end subroutine getnormpml_pa