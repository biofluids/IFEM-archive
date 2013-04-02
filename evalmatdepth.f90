subroutine evalmatdepth(matdepth, tfaction, basenode1, basenode2, solid_coor_init, ien_sbc, ien_solid, node_sbc)
! to calculate material depth in the very beginning of the simulation
! created by Jubiao Yang on 03/20/2013

use solid_variables
use fluid_variables, only: nsd

real(8) matdepth(nn_solid), tfaction(nn_solid)
integer basenode1(nn_solid), basenode2(nn_solid)
real(8) solid_coor_init(nsd_solid,nn_solid)

integer ien_sbc(ne_sbc, nen_solid+2)
integer ien_solid(ne_solid, nen_solid)
integer node_sbc(nn_sbc)

integer inode, ibs, idel
integer bcnode1(ne_sbc), bcnode2(ne_sbc)
real(8) tfaction_tmp(ne_sbc), dcandidate(ne_sbc)
real(8) x(nsd), xnode1(nsd), xnode2(nsd), vlineseg(nsd), vpointv1(nsd), vpointv2(nsd), vpointpj(nsd)

if (nsd==2) then                ! ONLY ABLE TO DEAL WITH 2D WITH THIS BLOCK
    matdepth(1:nn_solid)=0.0
    do inode=1,nn_solid
        do ibs=1,nn_sbc
            if (inode == node_sbc(ibs)) goto 888     ! no need to calculate for boundary nodes
        enddo

        x(1:nsd)=solid_coor_init(1:nsd,inode)
        do ibs=1,ne_sbc
            idel=ien_sbc(ibs,1)
            !----------------------------------------------------
            ! to find the boundary nodes in the element
            if (nen_solid==3) then            ! triangular elements
                if (ien_sbc(ibs,2)==0) then          ! nodes 2 & 3 are on the boundary
                    bcnode1(ibs)=2
                    bcnode2(ibs)=3
                elseif (ien_sbc(ibs,3)==0) then      ! nodes 1 & 3 are on the boundary
                    bcnode1(ibs)=1
                    bcnode2(ibs)=3
                elseif (ien_sbc(ibs,3)==0) then      ! nodes 1 & 2 are on the boundary
                    bcnode1(ibs)=1
                    bcnode2(ibs)=2
                endif
            elseif (nen_solid==4) then        ! quadrilateral elements
                if (ien_sbc(ibs,2)==1 .and. ien_sbc(ibs,3)==1) then          ! nodes 1 & 2 are on the boundary
                    bcnode1(ibs)=1
                    bcnode2(ibs)=2
                elseif (ien_sbc(ibs,3)==1 .and. ien_sbc(ibs,4)==1) then      ! nodes 2 & 3 are on the boundary 
                    bcnode1(ibs)=2
                    bcnode2(ibs)=3
                elseif (ien_sbc(ibs,4)==1 .and. ien_sbc(ibs,5)==1) then      ! nodes 3 & 4 are on the boundary
                    bcnode1(ibs)=3
                    bcnode2(ibs)=4
                elseif (ien_sbc(ibs,5)==1 .and. ien_sbc(ibs,2)==1) then      ! nodes 4 & 1 are on the boundary
                    bcnode1(ibs)=4
                    bcnode2(ibs)=1
                endif
            endif
            !-------------------------------------------------------
            bcnode1(ibs)=ien_solid(idel,bcnode1(ibs))
            bcnode2(ibs)=ien_solid(idel,bcnode2(ibs))
            xnode1(1:nsd)=solid_coor_init(1:nsd,bcnode1(ibs))
            xnode2(1:nsd)=solid_coor_init(1:nsd,bcnode2(ibs))
            vlineseg(1:nsd)=xnode2(1:nsd)-xnode1(1:nsd)
            vpointv1(1:nsd)=x(1:nsd)-xnode1(1:nsd)
            
            tfaction_tmp(ibs)=(vlineseg(1)*vpointv1(1)+vlineseg(2)*vpointv1(2))/(vlineseg(1)**2+vlineseg(2)**2)
            if (tfaction_tmp(ibs)<=0.0) then
                dcandidate(ibs)=sqrt(vpointv1(1)**2+vpointv1(2)**2)
                tfaction_tmp(ibs)=0.0
            elseif (tfaction_tmp(ibs)>=1.0) then
                vpointv2(1:nsd)=x(1:nsd)-xnode2(1:nsd)
                dcandidate(ibs)=sqrt(vpointv2(1)**2+vpointv2(2)**2)
                tfaction_tmp(ibs)=1.0
            elseif (tfaction_tmp(ibs)>0.0 .and. tfaction(ibs)<1.0) then
                vpointpj(1:nsd)=vpointv1(1:nsd)-tfaction_tmp(ibs)*vlineseg(1:nsd)
                dcandidate(ibs)=sqrt(vpointpj(1)**2+vpointpj(2)**2)
            endif
        enddo
        
        matdepth(inode)=minval(dcandidate(1:ne_sbc))
        do ibs=1,ne_sbc
            if (matdepth(inode)==dcandidate(ibs)) then
                basenode1(inode)=bcnode1(ibs)
                basenode2(inode)=bcnode2(ibs)
                tfaction(inode)=tfaction_tmp(ibs)
                goto 888
            endif
        enddo

        888 continue
    enddo
endif

return
end subroutine
