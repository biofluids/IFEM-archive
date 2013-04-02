subroutine materialdepth2d()
! to calculate material depth in the very beginning of the simulation
! created by Jubiao Yang on Mar. 9, 2013
!!!!!!!!!!!!!!!!!!! NOT FINISHED, TO BE REPLACED BY evalmatdepth.f90!!!!!!!!!!!!!!!!!!!!!!!!!!

use solid_variables
use mpi_variables, only: myid
implicit none

integer inode,ibcnode           ! iteration variables

do inode=1,nn_solid
    do ibcel=1,ne_sbc
        ! iterate through all boundary segments
        bcelglb=ien_sbc(ibcel,1)
        do ielnode=1,nen_solid
            elnodeglb=ien_solid(ine,ielnode)
            x(1:nsd,ielnode)=x_solid(1:nsd,elnodeglb)
        enddo
        if (nsd==2) then
            if (nen_solid==3) then            ! triangular elements
                if (ien_sbc(ibcel,2)==0) then          ! nodes 2 & 3 are on the boundary
                    bcnode1=2
                    bcnode2=3
                elseif (ien_sbc(ibcel,3)==0) then      ! nodes 1 & 3 are on the boundary 
                    bcnode1=1
                    bcnode2=3
                elseif (ien_sbc(ibcel,3)==0) then      ! nodes 1 & 2 are on the boundary
                    bcnode1=1
                    bcnode2=2
                endif
            elseif (nen_solid==4) then        ! quadrilateral elements
                if (ien_sbc(ibcel,2)==1 .and. ien_sbc(ibcel,3)==1) then          ! nodes 1 & 2 are on the boundary
                    bcnode1=1
                    bcnode2=2
                elseif (ien_sbc(ibcel,3)==1 .and. ien_sbc(ibcel,4)==1) then      ! nodes 2 & 3 are on the boundary 
                    bcnode1=2
                    bcnode2=3
                elseif (ien_sbc(ibcel,4)==1 .and. ien_sbc(ibcel,5)==1) then      ! nodes 3 & 4 are on the boundary
                    bcnode1=3
                    bcnode2=4
                elseif (ien_sbc(ibcel,5)==1 .and. ien_sbc(ibcel,2)==1) then      ! nodes 4 & 1 are on the boundary
                    bcnode1=4
                    bcnode2=1
                endif
            endif
        endif
        veclineseg(1:nsd)=x(1:nsd,bcnode2)-x(1:nsd,bcnode1)
        vecp
    enddo
enddo




end subroutine

