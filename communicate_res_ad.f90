subroutine communicate_res_ad(res,ndf,nn,send_address,ad_length)
use mpi_variables
include 'mpif.h'
! inputs

integer ndf
integer nn
real(8) res(ndf,nn)
integer ad_length
integer send_address(ad_length,2)

! inside variables
real(8) recv_tmp(ndf+1,ad_length)
integer icount
integer tag
integer size
integer req(ad_length)
!integer req(ad_length*2)
integer status(mpi_status_size,ad_length)
real(8) sendbuf(ndf+1)
integer node

tag=1
      call mpi_barrier(mpi_comm_world,ierror)

do icount=1,ad_length
	call mpi_irecv(recv_tmp(1,icount),ndf+1,mpi_double_precision,&
		mpi_any_source,mpi_any_tag,mpi_comm_world,req(icount),ierror)
	
!	sendbuf(1:ndf)=res(1:ndf,send_address(icount,2))
!	sendbuf(ndf+1)=send_address(icount,2)
!       call mpi_isend(sendbuf,ndf+1,mpi_double_precision,&
!                 send_address(icount,1), tag, mpi_comm_world,req(icount*2),ierror)
end do

do icount=1,ad_length
!       write(*,*) 'send buf',res(1,send_address(icount,2))
        sendbuf(1:ndf)=res(1:ndf,send_address(icount,2))
        sendbuf(ndf+1)=send_address(icount,2)
       call mpi_rsend(sendbuf(1),ndf+1,mpi_double_precision,&
                 send_address(icount,1), tag, mpi_comm_world,ierror)
end do

call mpi_waitall(ad_length,req,status,ierror)



!write(*,*) 'ierror', ierror, 'myid', myid

do icount=1,ad_length
	node=nint(recv_tmp(ndf+1,icount))
	res(1:ndf,node)=res(1:ndf,node)+recv_tmp(1:ndf,icount)
end do

      call mpi_barrier(mpi_comm_world,ierror)

!if (myid ==0) then
!write(*,*) 'recv_tmp', res(1:ndf,:)
!end if

end subroutine communicate_res_ad 
