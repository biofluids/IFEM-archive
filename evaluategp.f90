subroutine evaluategp(xcoortest,gpenetr,lambdacf,ien_sbc,ien_solid,xinit)

use solid_variables

real(8) xcoortest(nsd_solid,nn_solid)
real(8) xinit(nsd_solid,nn_solid)
real(8) gpenetr(nn_solid)
real(8) lambdacf(nn_solid)
integer ien_sbc(ne_sbc,nen_solid+2)
integer ien_solid(ne_solid,nen_solid)

real(8) yUpperLimit,diffUyLimit
integer iCase,ibs,idel,ibe,idnode,ibstg,ideltg,ibetg,idnodetg

real(8) vecmastg(nsd_solid)
real(8) outnorml(nsd_solid)
real(8) normvec,cosang,cosangmax
real(8) tmpx(nen_solid,nsd)
integer tmp_ien(nen_solid)

iCase=1
yUpperLimit=1.397

if (iCase==1) then
    !++++++++++++++++++++++++
    ! Half-space
    do ibs=1,ne_sbc
        idel=ien_sbc(ibs,1)
        do ibe=1,nen_solid
            idnode=ien_solid(idel,ibe)
            diffUyLimit=xcoortest(2,idnode)-yUpperLimit
!	    if (myid == 0) write(*,*) '(x,y)', xcoortest(1,idnode), xcoortest(2,idnode)
            if (diffUyLimit>0) then
                gpenetr(idnode)=diffUyLimit
            else
                gpenetr(idnode)=0.0
            endif
        enddo
    enddo
    !+++++++++++++++++++++++++
else
    !+++++++++++++++++++++++++
    ! Full-space, aka self-penetration
    ! !!!!!!!!!!!!!! Haven't identified the master and slave surfaces yet
    do ibs=1,ne_sbc
        idel=ien_sbc(ibs,1)
        !......................................
        tmp_ien(1:nen_solid)=ien_sbc(ibs,2:nen_solid+1)
        do isd=1,nsd
            do ibe=1,nen_solid
                idnode=ien_solid(idel,ibe)
                tmpx(ibe,isd)=xcoortest(isd,idnode)
            enddo
        enddo
        !......................................
        do ibe=1,nen_solid
            idnode=ien_solid(idel,ibe)
            !------------------ iterate around target bc nodes
            cosangmax=0.0
            do ibstg=1,ne_sbc
                ideltg=ien_sbc(ibstg,1)
                do ibetg=1,nen_solid
                    idnodetg=ien_solid(ideltg,ibetg)
                    vecmastg=xcoortest(:,idnode)-xcoortest(:,idnodetg)
                    call getnorm=(vecmastg,vecmastg,nsd_solid,normvec)
                    normvec=sqrt(normvec)
                    vecmastg=vecmastg/normvec
                    !...................
                    tmp_ien(1:nen_solid)=ien_sbc(ibs,2:nen_solid+1)
                    do isd=1,nsd
                        do jsd=1,nen_solid
                            tmpx(jsd,isd)=x(isd,jsd)
                        enddo
                    enddo
                    call outnormal_2d(tmpx,tmp_ien,out_norm,nsd,nen_solid,w)
                    !....................
                    call outnormal_2d( )
                    cosang=dot_product(vecmastg,outnorml)
                    if (cosang>cosangmax) then
                        cosangmax=cosang
                        gpenetr=normvec
                    endif
                enddo
            enddo
            !------------------
            if (cosangmax>0.0) then
                gpenetr(idnode)=gpenetr(idnode)/2.0
            else
                gpenetr(idnode)=0.0
            endif
        enddo
    enddo
    !+++++++++++++++++++++++++
endif


return
end

