!==================================================================
! J. Yang
! Rensselaer Polytechnic Institute
! This subroutine initializes the parameters for Perfectly Matched Layer
!------------------------------------------------------------------
! Latest update: Jack, 07/22/2014
!==================================================================
subroutine initializePMLparam(xloc,node_local,nn_local,send_adress,ad_length)
  use global_constants
  use fluid_variables
  use mpi_variables
  use pml_variables
  implicit none

  real(8) :: xloc(nsd,nn)
  real(8) :: DPML, distPML
  integer :: irng, inn, isd, iad
  integer :: nn_local, error_id, nodecount, node
  integer :: node_local(nn_local)
  integer ad_length
  integer send_adress(ad_length,2)
  integer b, j, count, tmp(1,2), kcount

    if (allocated(sigmaPML)) then
        deallocate(sigmaPML)
    endif
    if (allocated(qv)) then
        deallocate(qv)
    endif
    if (allocated(qvold)) then
        deallocate(qvold)
    endif
    if (allocated(seqcPML)) then
        deallocate(seqcPML)
    endif
    if (allocated(sendpml_addrs)) then
        deallocate(sendpml_addrs)
    endif
    if (allocated(subpml_addrs)) then
        deallocate(subpml_addrs)
    endif

    allocate(sigmaPML(nsd,nn),stat=error_id)
    allocate(qv(ndf,nn)      ,stat=error_id)
    allocate(qvold(ndf,nn)   ,stat=error_id)
    allocate(seqcPML(nn)     ,stat=error_id)

    sigmaPML(:,:)=0.0
    qv(:,:)=0.0
    qvold(:,:)=0.0
    seqcPML(:)=0
    DPML=nDPML*hmax
    nn_PML=0
    nn_PML_local=0

    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    if (sumNbcPML .ne. 0) then

        nodecount=0
        do irng=1,nrng
            isd=flagPML(irng)
            if (isd .ne. 0) then
                do inn=1,nn
                    distPML=abs( xloc(isd,inn)-xyzcPML(irng) )
                    if ( distPML < DPML ) then
                        sigmaPML(isd,inn) = sigmaMaxPML(isd) * (1.0 - distPML/DPML)**4.0
                        if (seqcPML(inn)==0) then
                            nodecount=nodecount+1
                            seqcPML(inn)=nodecount
                        endif
                    endif
                enddo
            endif
        enddo

        nn_PML=nodecount
        if (allocated(nodePML)) then
            deallocate(nodePML)
        endif
        allocate(nodePML(nn_PML),stat=error_id)
        nodePML(1:nn_PML)=0
        do inn=1,nn
            if (seqcPML(inn) .ne. 0) nodePML(seqcPML(inn))=inn
        enddo

        if (allocated(seqcPMLlocal)) then
            deallocate(seqcPMLlocal)
        endif
        allocate(seqcPMLlocal(nn_local),stat=error_id)
        seqcPMLlocal(1:nn_local)=0
        nodecount=0
        do inn=1,nn_local
            if (seqcPML(node_local(inn)) .ne. 0) then
                nodecount=nodecount+1
                seqcPMLlocal(inn)=nodecount
            endif
        enddo
        nn_PML_local=nodecount
        !if (myid==9) write(*,*) "myid=", myid, ", nn_PML_local=", nn_PML_local, "seqcPMLlocal=", seqcPMLlocal

        !-------------- communication info for pml --------------------
        !||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        adpml_length=0
        do iad=1,ad_length
            node=send_adress(iad,2)
            if (seqcPML(node) > 0) then
                adpml_length=adpml_length+1
            endif
        enddo

        if (adpml_length > 0) then

            allocate(sendpml_addrs(adpml_length,2))
            sendpml_addrs(:,:)=0
            kcount=0

            do iad=1,ad_length
                node=send_adress(iad,2)
                if (seqcPML(node) > 0) then
                    kcount=kcount+1
                    sendpml_addrs(kcount,1)=send_adress(iad,1)
                    sendpml_addrs(kcount,2)=send_adress(iad,2)
                endif
            enddo

            !--------------------- subpml_addrs -----------------------
            do iad=adpml_length,1,-1
                do j=1,iad-1
                    if (sendpml_addrs(j,1) .lt. sendpml_addrs(j+1,1)) then
                        tmp(1,:) = sendpml_addrs(j,:)
                        sendpml_addrs(j,:) = sendpml_addrs(j+1,:)
                        sendpml_addrs(j+1,:) = tmp(1,:)
                    endif
                enddo
            enddo

            b=sendpml_addrs(1,1)
            countrowpml=1
            do iad=1,adpml_length
                if (sendpml_addrs(iad,1) .ne. b) then
                    countrowpml=countrowpml+1
                    b=sendpml_addrs(iad,1)
                endif
            enddo
            
            allocate(subpml_addrs(countrowpml,2))

            b=sendpml_addrs(1,1)
            count=0
            countrowpml=1
            do iad=1,adpml_length
                if (sendpml_addrs(iad,1) .eq. b) then
                    count=count+1
                else
                    subpml_addrs(countrowpml,1)=b
                    subpml_addrs(countrowpml,2)=count
                    countrowpml=countrowpml+1
                    b=sendpml_addrs(iad,1)
                    count=1
                endif
            enddo
            subpml_addrs(countrowpml,1)=b
            subpml_addrs(countrowpml,2)=count
            !-------------------------------------------------------
            !if (myid==5) then
            !    write(*,*) "myid=", myid, ", sendpml_addrs=", sendpml_addrs
            !    write(*,*) "myid=", myid, ", subpml_addrs=", subpml_addrs
            !endif
        endif

    endif
    !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

return
end subroutine initializePMLparam