!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Jack Jubiao Yang, 07/10/2014
! Rensselaer Polytechnic Institute
! This module defines most variables used by PML algorithm
! This module doesn't affect the standard fluid features even though PML is not applied to any BC
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module pml_variables
    implicit none
    save

    ! the sizes of flagPML and xyzcPML was maxnsurf, here hardcoded as 21
    integer :: flagPML(21)      ! flag if PML is defined on this boundary
    integer :: sumNbcPML
    real(8) :: xyzcPML(21)      ! coordinate of the maximum coordinate

    integer :: nn_PML, nn_PML_local
    integer,allocatable :: nodePML(:)           ! nodePML(nn_PML)
    integer,allocatable :: seqcPMLlocal(:)      ! seqcPMLlocal(nn_local)
    integer,allocatable :: seqcPML(:)           ! seqcPML(nn)

    real(8),allocatable :: qv(:,:)              ! qv(ndf,nn)  ::::: qv(nsd-1,ndf,nn)
    real(8),allocatable :: qvold(:,:)              ! qv(ndf,nn)  ::::: qv(nsd-1,ndf,nn)
    real(8),allocatable :: sigmaPML(:,:)        ! sigmaPML(nsd,nn)

    integer,allocatable :: sendpml_addrs(:,:)   ! sendpml_addrs(adpml,2)
    integer,allocatable :: subpml_addrs(:,:)    ! subpml_addrs(countrowpml,2)
    integer adpml_length, countrowpml

end module pml_variables