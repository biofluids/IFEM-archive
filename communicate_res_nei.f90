subroutine communicate_res_nei(res,ndf,nn,nei_max,nn_local_com,global_nei,global_com,nn_global_com,local_com)
use mpi_variables
include 'mpif.h'
! inputs
integer ndf
real(8) res(ndf,nn)

integer nei_max
integer nn_local_com
integer global_nei(nn_global_com,nei_max)
integer global_com(nn_global_com)
integer nn_global_com
integer local_com(nn_local_com)
! inside variables
integer new_group
integer orig_group
integer new_comm

integer icount
integer jcount
integer kcount
integer node
integer flag
real(8) recv(ndf)
real(8) send(ndf)

integer new_rank



call MPI_COMM_GROUP(MPI_COMM_WORLD,orig_group,ierr)

do icount=1,nn_global_com
flag=-999
recv(:)=0.0d0

	do jcount=1,nei_max
		kcount=jcount
		if (global_nei(icount,jcount) .eq. myid) then
			flag=1
		end if 

		if (global_nei(icount,jcount) .eq. -999) then
			kcount=kcount-1
			goto 100
		end if
	end do
100 continue
	call MPI_GROUP_INCL(orig_group,kcount,global_nei(icount,1:kcount),new_group,ierror)
	call MPI_COMM_CREATE(MPI_COMM_WORLD,new_group,new_comm,ierror)
!	call MPI_GROUP_RANK(new_group,new_rank,ierror)
	if (flag .eq. 1) then
	send(1:ndf)=res(1:ndf,global_com(icount))
	call MPI_ALLREDUCE(send(1),recv(1),ndf,mpi_double_precision,mpi_sum,new_comm,ierror)
	res(1:ndf,global_com(icount))=recv(1:ndf)
!	write(*,*) 'myid',myid, 'new_rank',new_rank,'finished','kcount',kcount
	call MPI_GROUP_FREE(new_group,ierror)
	call MPI_COMM_FREE(new_comm,ierror)
	end if
end do

!do icount=1,nn_local_com
!	node=global_com(local_com(icount))
!	res(1:ndf,node)=recv(1:ndf,local_com(icount))
!end do

!if (myid == 0) then
!	write(*,*) 'in com_nei res'
!	do icount=1,nn
!		write(*,*) res(:,icount)
!	end do
!end if


end subroutine communicate_res_nei

