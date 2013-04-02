subroutine evaluategp(xcoortest,gpenetr,lambdacf,ien_sbc,ien_solid)

use solid_variables

real(8) xcoortest(nsd_solid,nn_solid)
real(8) gpenetr(nn_solid)
real(8) lambdacf(nn_solid)
integer ien_sbc(ne_sbc,nen_solid+2)
integer ien_solid(ne_solid,nen_solid)

real(8) yUpperLimit,diffUyLimit
integer iCase,ibs,idel,ibe,idnode

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
endif



return
end

