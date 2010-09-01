subroutine communicate_res(global_com,nn_global_com,local_com,nn_local_com &
			,res,ndf,nn)
! collect the Residual on the nodes which are shared by different Processors
! Sum up these parts of residual accordingly then distribute them back to their own Processors 
use mpi_variables
implicit none
include 'mpif.h'
integer nn_global_com
integer global_com(nn_global_com)  ! global node index for communication
integer nn_local_com
integer local_com(nn_local_com)  ! local index in the communication region on each processor
integer ndf
integer nn
real(8) res(ndf,nn) ! Residual vector needs to be communicated
real(8) bus(ndf,nn_global_com) ! communication bus
real(8) bus_rec(ndf,nn_global_com)
integer i
integer node


bus(:,:)=0.0d0
bus_rec(:,:)=0.0d0

do i=1,nn_global_com
	node=global_com(i)
	bus(1:ndf,i)=res(1:ndf,node)
end do 
! Put value needs to be communicated on the bus

           call mpi_barrier(mpi_comm_world,ierror)
           call mpi_allreduce(bus(1,1),bus_rec(1,1),ndf*nn_global_com,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
! Sum up all the residuals on bus to Proc 0 
!           call mpi_bcast(bus_rec(1,1),ndf*nn_global_com,mpi_double_precision,0,mpi_comm_world,ierror)
! Set bus_res back to all processors

do i=1,nn_local_com
	node=global_com(local_com(i))
	res(1:ndf,node)=bus_rec(1:ndf,local_com(i))
end do
! Take global residual from bus, and put them back onto global node position







end subroutine communicate_res
