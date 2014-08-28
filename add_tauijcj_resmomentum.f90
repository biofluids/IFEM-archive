! tempc(1:nsd)=res_t(1:nsd)
do inl=1,nen
    node=ien(inl,ie)
    do isd=1,nsd
        do jsd=1,nsd
            res_t(isd) = res_t(isd) - sh(jsd,inl)*hattauij(isd,jsd,node)
        enddo
    enddo
enddo

!            if (node==5000) then
!                write (*,*) tempc(1:nsd)
!                write (*,*) "residual", res_t(1:nsd)
!            endif