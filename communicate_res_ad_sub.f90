subroutine communicate_res_ad_sub(res,ndf,nn,send_address,ad_length)
use mpi_variables
implicit none
include 'mpif.h'
! inputs
integer ndf
integer nn
real(8) res(ndf,nn)
integer ad_length
integer send_address(ad_length,2)

! inside variables
real(8) recv_tmp(ndf*ad_length)
real(8) send_tmp(ndf*ad_length)
integer icount
integer tag
integer size
integer reqr(countrow)
integer reqs(countrow)
integer status(mpi_status_size,countrow)
integer stat(mpi_status_size)
real(8),allocatable,dimension(:,:) ::sendbuf
real(8),allocatable,dimension(:,:) ::recbuf
integer node

integer i,j,k1,k2
integer des

tag=1
      call mpi_barrier(mpi_comm_world,ierror)
!if (mod(myid,2) .eq. 0) then

do i=1,ad_length
do j=1,ndf
send_tmp(j+ (i-1)*ndf ) = res(j,send_address(i,2))
end do
end do
k1=0
k2=0




do icount=1,countrow
	i=sub_address(icount,2) ! get envelope size
	des=sub_address(icount,1) ! get envelope destination

!	allocate (sendbuf(ndf+1,i)) ! declare envelope
	k1=k2+1
	k2=k2+i


!	sendbuf(1:ndf,1:i)=res(1:ndf,send_address(k1:k2,2))
!        sendbuf(ndf+1,1:i)=send_address(k1:k2,2) ! make envelope

!       call mpi_isend(sendbuf(1,1),(ndf+1)*i,mpi_double_precision,&
!                 des, tag, mpi_comm_world,reqs(icount),ierror)  ! send envelope

!if (myid == 0) write(*,*) 'k1, k2', k1,k2

	call mpi_isend(send_tmp((k1-1)*ndf +1),ndf*i,mpi_double_precision,&
                        des, tag, mpi_comm_world,reqs(icount),ierror)  ! send envelope
! changed from "mpi_ibsend" to "mpi_isend" by Jubiao Yang on 03/19/2013


	call mpi_irecv(recv_tmp((k1-1)*ndf +1),ndf*i,mpi_double_precision,&
        		des,tag,mpi_comm_world,reqr(icount),ierror)  ! receive envelope

!	deallocate(sendbuf)

end do


!do icount=1,countrow
!        i=sub_address(icount,2) ! get envelope size
!	des=sub_address(icount,1) ! get envelope destination
!        allocate (recbuf(ndf+1,i)) ! declare envelope
!
!    call mpi_recv(recbuf(1,1),(ndf+1)*i,mpi_double_precision,&
!            des,tag,mpi_comm_world,stat,ierror)  ! receive envelope
!
!	do j=1,i
!	       node=nint(recbuf(ndf+1,j))
!	       res(1:ndf,node)=res(1:ndf,node)+recbuf(1:ndf,j)
!	end do
!
!	deallocate(recbuf)
!end do

call mpi_waitall(countrow,reqs,status,ierror)
call mpi_waitall(countrow,reqr,status,ierror)

do icount=1,ad_length
	node=send_address(icount,2)
	do j=1,ndf
	res(j,node)=res(j,node)+recv_tmp(j+(icount-1)*ndf)
	end do
end do


      call mpi_barrier(mpi_comm_world,ierror)


end subroutine communicate_res_ad_sub
