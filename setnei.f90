subroutine setnei(global_com,nn_global_com,local_com,nn_local_com, &
		  local_nei,nei_max,nei_global,ad_length)
! This subroutine is to get the index of processors who sharing the same node


use mpi_variables
implicit none
include 'mpif.h'
! inputs
integer nn_global_com
integer global_com(nn_global_com)
integer nn_local_com
integer local_com(nn_local_com)
integer nei_max

! output local neighbour matrix
integer local_nei(nn_local_com,nei_max+1)
integer ad_length

! inside subroutine variables
integer inn
integer icount
integer jcount
integer id
integer inl
integer flag
integer kcount
integer nei(nn_global_com,ncpus)
integer nei_rec(nn_global_com,ncpus)
integer nei_global(nn_global_com,nei_max)


flag=-999
nei(:,:)=0
nei_global(:,:)=-999


! Loop over nn_global_com, find nodes in this set contained by which processor
do inn=1,nn_global_com
	id=global_com(inn)
	icount=1
	flag=-999
	kcount=2

loop:	do while((icount .le. nn_local_com) .and. (flag .lt. 0))
		if (inn == local_com(icount)) then
			flag = myid
		end if
		icount=icount+1
		continue
	end do loop 

if (flag .ne. -999) then
	nei(inn,myid+1)=1 ! before communication just record the information in a huge matrix
end if

end do

           call mpi_barrier(mpi_comm_world,ierror)
           call mpi_reduce(nei(1,1),nei_rec(1,1),nn_global_com*ncpus,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
! Collect the information all onto master processor 
! Then reduce such huge matrix to a smaller table to save space on master processor
if (myid == 0) then
	do inn=1,nn_global_com
		jcount=1
		do icount=1,ncpus
		if (nei_rec(inn,icount) == 1) then
		nei_global(inn,jcount)=icount-1
		jcount=jcount+1
		end if
		end do
		id =global_com(inn)
!		write(*,*) 'node id', id, 'nei_global', nei_global(inn,:)
	end do
end if
! Pass the redused table to all processors and every processor take its own portion to set up the local table 
           call mpi_bcast(nei_global(1,1),nn_global_com*nei_max,mpi_integer,0,mpi_comm_world,ierror)

ad_length=0

do icount=1,nn_local_com
	local_nei(icount,1)= global_com(local_com(icount))
	local_nei(icount,2:nei_max+1)= nei_global(local_com(icount),1:nei_max)
	do jcount=2,nei_max+1
		if ((local_nei(icount,jcount) .ne. myid) .and. (local_nei(icount,jcount) .ne. -999 )) & 
		ad_length=ad_length+1
	end do
end do

!if (myid == 3) then
!do icount=1,nn_local_com
!write(*,*) 'local_nei', local_nei(icount,:)
!end do
!end if

!write(*,*) 'ad_length', ad_length, 'myid', myid
           call mpi_barrier(mpi_comm_world,ierror)





end subroutine setnei
