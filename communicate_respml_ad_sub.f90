subroutine communicate_respml_ad_sub(res,ndf,nn,send_address,ad_length)
use mpi_variables
use pml_variables
implicit none
include 'mpif.h'
! inputs
integer ndf
integer nn
real(8) res(ndf,nn+nn_PML)
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
integer node
integer i,j,k1,k2
integer des
!------------------ PML local variables ----------------
real(8) recv_pml(ndf*adpml_length), send_pml(ndf*adpml_length)
integer reqrpml(countrowpml), reqspml(countrowpml)
integer statuspml(mpi_status_size,countrowpml)
!-------------------------------------------------------

    tag=1
    call mpi_barrier(mpi_comm_world,ierror)

    !-------------- communicating velocity and pressure -----------
    !||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    do i=1,ad_length
        do j=1,ndf
            send_tmp(j+ (i-1)*ndf ) = res(j,send_address(i,2))
        enddo
    enddo
    k1=0
    k2=0

    do icount=1,countrow
        i=sub_address(icount,2) ! get envelope size
        des=sub_address(icount,1) ! get envelope destination
        k1=k2+1
        k2=k2+i
	    call mpi_isend(send_tmp((k1-1)*ndf +1),ndf*i,mpi_double_precision,&
                        des, tag, mpi_comm_world,reqs(icount),ierror)  ! send envelope
	    call mpi_irecv(recv_tmp((k1-1)*ndf +1),ndf*i,mpi_double_precision,&
        		        des,tag,mpi_comm_world,reqr(icount),ierror)  ! receive envelope
    enddo

    call mpi_waitall(countrow,reqs,status,ierror)
    call mpi_waitall(countrow,reqr,status,ierror)

    do icount=1,ad_length
        node=send_address(icount,2)
        do j=1,ndf
            res(j,node)=res(j,node)+recv_tmp(j+(icount-1)*ndf)
        enddo
    enddo
    !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    call mpi_barrier(mpi_comm_world,ierror)

    !------------ communicating PML auxiliary variables -----------
    !||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    if (adpml_length > 0) then
        do i=1,adpml_length
            do j=1,ndf
                send_pml(j+ (i-1)*ndf) = res(j,nn+seqcPML(sendpml_addrs(i,2)))
            enddo
        enddo
        k1=0
        k2=0

        do icount=1,countrowpml
            i=subpml_addrs(icount,2)         ! get envelope size
            des=subpml_addrs(icount,1)       ! get envelope destination
            k1=k2+1
            k2=k2+i
            call mpi_isend(send_pml((k1-1)*ndf +1),ndf*i,mpi_double_precision,&
                            des,tag,mpi_comm_world,reqspml(icount),ierror)  ! send envelope
            call mpi_irecv(recv_pml((k1-1)*ndf +1),ndf*i,mpi_double_precision,&
                            des,tag,mpi_comm_world,reqrpml(icount),ierror)  ! receive envelope
        enddo

        call mpi_waitall(countrowpml,reqspml,statuspml,ierror)
        call mpi_waitall(countrowpml,reqrpml,statuspml,ierror)

        do icount=1,adpml_length
            node=sendpml_addrs(icount,2)
            do j=1,ndf
                res(j,nn+seqcPML(node))=res(j,nn+seqcPML(node))+recv_pml(j+(icount-1)*ndf)
            enddo
        enddo
    endif
    !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    call mpi_barrier(mpi_comm_world,ierror)

end subroutine communicate_respml_ad_sub