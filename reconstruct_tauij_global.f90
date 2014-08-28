!================================================

!================================================
hattauij(1:nsd,1:nsd,1:nn)=0.0
do ie=1,ne     ! loop over elements
    do inl=1,nen
        x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
        d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
    enddo

    do iq=1,nquad  ! loop over the quadrature points in each element 
        if (nsd==2) then
            if (nen.eq.3) then !calculate shape function at quad point
                include "sh2d3n.h"
            elseif (nen.eq.4) then
                include "sh2d4n.h"
            endif
        elseif (nsd==3) then
            if (nen.eq.4) then !calculate shape function at quad point
                include "sh3d4n.h"
            elseif (nen.eq.8) then
                include "sh3d8n.h"
            endif
        endif
        eft0 = abs(det) * wq(iq) ! calculate the weight at each quad pt

        dr(1:nsd,1:ndf)=0.0    ! dd/dxi
        mu=0.0

        do inl=1,nen
            tempc(1:nsd) = ama*d(1:nsd,inl)+oma*d_old(1:nsd,inl)
            mu = mu+sh(0,inl)*fvis(ien(inl,ie))
            do isd=1,nsd
                dr(isd,1:nsd) = dr(isd,1:nsd)+sh(isd,inl)*tempc(1:nsd)
            enddo
        enddo
        !mu=1.8e-4

        do isd = 1,nsd
            do jsd = 1,nsd
                tau(isd,jsd) = mu*(dr(isd,jsd) + dr(jsd,isd))
            enddo
            if (nsd==2) then
                tau(isd,isd)=tau(isd,isd)-2.0*mu*(dr(1,1)+dr(2,2))/3.0
            elseif (nsd==3) then
                tau(isd,isd)=tau(isd,isd)-2.0*mu*(dr(1,1)+dr(2,2)+dr(3,3))/3.0
            endif
        enddo

        do inl=1,nen ! loop over number of nodes in an element 
            node = ien(inl,ie)
            do isd=1,nsd
                do jsd=1,nsd
                    hattauij(isd,jsd,node)=hattauij(isd,jsd,node)+&
                                           eft0*sh(0,inl)*tau(isd,jsd)
                enddo
            enddo
            if (abs(xloc(ysd,node)-1.397) < 1e-5) then
                hattauij(1,2,node)=0.0
                hattauij(2,1,node)=0.0
            endif
        enddo
    enddo ! end of qudrature pts loop
enddo ! end of element loop
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! ||||||||||||||||||||||||||||||||||||||||||||||||
! evaluating \R{\tau_{ij}}

! evaluating \hat{\tau_{ij}} = \R{\tau_{ij}} / M_A
! ||||||||||||||||||||||||||||||||||||||||||||||||
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
flagdivision(1:nn)=0
do ie=1,ne
    do inl=1,nen
        node=ien(inl,ie)
        if (flagdivision(node)==0) then
            hattauij(1:nsd,1:nsd,node)=hattauij(1:nsd,1:nsd,node)/mpqlumpmass(node)
            flagdivision(node)=1
        endif
    enddo
enddo




