subroutine setnei_new(global_com,nn_global_com,local_com,nn_local_com, &
		  local_nei,nei_max,ad_length)
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
integer nei(ncpus)
integer nei_rec(ncpus)
integer nei_global(nn_global_com,nei_max)


flag=-999
nei(:)=0
nei_rec(:)=0

nei_global(:,:)=-999


! Loop over nn_global_com, find nodes in this set contained by which processor
do inn=1,nn_global_com
	id=global_com(inn)
	icount=1
	flag=-999
	kcount=2

    loop: do while((icount .le. nn_local_com) .and. (flag .lt. 0))
		if (inn == local_com(icount)) then
			flag = myid
		endif
		icount=icount+1
		continue
	enddo loop 

    if (flag .ne. -999) then
	    nei(myid+1)=1 ! before communication just record the information in a huge matrix
    endif

    call mpi_barrier(mpi_comm_world,ierror)
    call mpi_allreduce(nei(1),nei_rec(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)
! Collect the information all onto master processor
    jcount=1
    do icount=1,ncpus
        if (nei_rec(icount) == 1) then
            nei_global(inn,jcount)=icount-1
            jcount=jcount+1
        endif
    enddo
	nei(:)=0
	nei_rec(:)=0
! Reduce this matrix to nei_global on every processor
enddo

ad_length=0

do icount=1,nn_local_com
	local_nei(icount,1)= global_com(local_com(icount))
	local_nei(icount,2:nei_max+1)= nei_global(local_com(icount),1:nei_max)
	do jcount=2,nei_max+1
		if ((local_nei(icount,jcount) .ne. myid) .and. (local_nei(icount,jcount) .ne. -999 )) & 
		ad_length=ad_length+1
	enddo
enddo

call mpi_barrier(mpi_comm_world,ierror)

end subroutine setnei_new