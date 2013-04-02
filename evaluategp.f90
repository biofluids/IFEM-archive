subroutine evaluategp(xcoortest,gpenetr,dirpenetr,ien_sbc,ien_solid,matdepth,tfaction,basenode1,basenode2,node_sbc)

use solid_variables
use fluid_variables, only: maxconn, nsd

real(8) xcoortest(nsd_solid,nn_solid)
real(8) xsbc(nsd_solid)
real(8) gpenetr(nn_solid)
real(8) dirpenetr(nsd_solid,nn_solid)
integer ien_sbc(ne_sbc,nen_solid+2)
integer ien_solid(ne_solid,nen_solid)
integer node_sbc(nn_sbc)

real(8) Nshpfunc(nen_solid)
real(8) xn(nen_solid), yn(nen_solid), xx, yy
real(8) dirct(nsd_solid), norm_dirct

real(8) matdepth(nn_solid), tfaction(nn_solid)
integer basenode1(nn_solid), basenode2(nn_solid)

real(8) yUpperLimit,diffUyLimit
integer iCase,ibs,idel,ibe,idnode,isd, compnode, el_enclose
integer INFO
real(8) Amatrix(nen_solid,nen_solid), Bvector(nen_solid), IP(nen_solid)
real(8) sumshpfunc

iCase=2

iHalfDir=xdf
yUpperLimit=8.0      ! x-dir

if (iCase==1 .or. iCase==3) then
    !++++++++++++++++++++++++
    ! Half-space
    do ibs=1,ne_sbc
        idel=ien_sbc(ibs,1)
        do ibe=1,nen_solid
            idnode=ien_solid(idel,ibe)
            diffUyLimit=xcoortest(iHalfDir,idnode)-yUpperLimit
!	    if (myid == 0) write(*,*) '(x,y)', xcoortest(1,idnode), xcoortest(2,idnode)
            if (diffUyLimit>0.0) then
                gpenetr(idnode)=diffUyLimit
                dirpenetr(1:nsd_solid,idnode)=0.0
                dirpenetr(iHalfDir,idnode)=1.0
            else
                gpenetr(idnode)=0.0
            endif
        enddo
    enddo
    !+++++++++++++++++++++++++
endif
if (iCase==2 .or. iCase==3) then
    !+++++++++++++++++++++++++ added by Jubiao Yang on Mar. 13, 2013
    ! Self-contact
    do ibs=1,nn_sbc
        idnode=node_sbc(ibs)
        xsbc(1:nsd_solid)=xcoortest(1:nsd_solid,idnode)
        el_enclose=0
        call getinf_el_solid(el_enclose, xsbc, xcoortest, nn_solid, nsd_solid, ne_solid, nen_solid, ien_solid, maxconn)
        !-------------------------------------------------------
        if (el_enclose/=0) then
            ! calculate the depth based on information on el_enclose
            xx=xsbc(1)
            yy=xsbc(2)
            do ibe=1,nen_solid
                xn(ibe)=xcoortest(1,ien_solid(el_enclose,ibe))
                yn(ibe)=xcoortest(2,ien_solid(el_enclose,ibe))
            enddo
            ! to calculate shape functions
            if (nen_solid==3) then
                Nshpfunc(1)=(xx*yn(2) - xn(2)*yy - xx*yn(3) + xn(3)*yy + xn(2)*yn(3) - xn(3)*yn(2))/ &
                            (xn(1)*yn(2) - xn(2)*yn(1) - xn(1)*yn(3) + xn(3)*yn(1) + xn(2)*yn(3) - xn(3)*yn(2))
                Nshpfunc(2)=-(xx*yn(1) - xn(1)*yy - xx*yn(3) + xn(3)*yy + xn(1)*yn(3) - xn(3)*yn(1))/ &
                            (xn(1)*yn(2) - xn(2)*yn(1) - xn(1)*yn(3) + xn(3)*yn(1) + xn(2)*yn(3) - xn(3)*yn(2))
                Nshpfunc(3)=(xx*yn(1) - xn(1)*yy - xx*yn(2) + xn(2)*yy + xn(1)*yn(2) - xn(2)*yn(1))/ &
                            (xn(1)*yn(2) - xn(2)*yn(1) - xn(1)*yn(3) + xn(3)*yn(1) + xn(2)*yn(3) - xn(3)*yn(2))
            elseif (nen_solid==4) then
                ! this elseif block needs to be reviewed and corrected: bugs had been detected
                Bvector(1)=1.0
                Bvector(2)=xx
                Bvector(3)=yy
                Bvector(4)=xx*yy
                do ibe=1,nen_solid
                    Amatrix(1,ibe)=1.0
                    Amatrix(2,ibe)=xn(ibe)
                    Amatrix(3,ibe)=yn(ibe)
                    Amatrix(4,ibe)=xn(ibe)*yn(ibe)
                enddo
                call DGESV(4,4,Amatrix,4,IP,Bvector,4,INFO)
                Nshpfunc(1:nen_solid)=Bvector(1:nen_solid)
            endif
            gpenetr(idnode)=0.0
            dirpenetr(1:nsd_solid,idnode)=0.0
            do ibe=1,nen_solid
                compnode=ien_solid(el_enclose,ibe)
                gpenetr(idnode)=gpenetr(idnode)+Nshpfunc(ibe)*matdepth(compnode)
                if (matdepth(compnode)/=0.0 .and. basenode1(compnode)/=0 .and. basenode2(compnode)/=0) then
                    dirct(1:nsd_solid)=xcoortest(1:nsd_solid,compnode) &
                                 -(1-tfaction(compnode))*xcoortest(1:nsd_solid,basenode1(compnode)) &
                                 -tfaction(compnode)*xcoortest(1:nsd_solid,basenode2(compnode))
                endif
                norm_dirct=sqrt(dirct(1)**2+dirct(2)**2)
                dirct(1:nsd_solid)=dirct(1:nsd_solid)/norm_dirct
                dirpenetr(1:nsd_solid,idnode)=dirpenetr(1:nsd_solid,idnode)+Nshpfunc(ibe)*dirct(1:nsd_solid)
            enddo
            !sumshpfunc=sum(Nshpfunc)
            !if (sumshpfunc/=1.0 .and. sumshpfunc/=0.0) then
            !    write(*,*) 'SUM OF SHAPE FUNCTIONS', sumshpfunc
            !    write(*,*) 'ENCLOSING ELEMENT', el_enclose
            !endif
            gpenetr(idnode)=gpenetr(idnode)/2.0
        endif
        !-------------------------------------------------------
    enddo
    !+++++++++++++++++++++++++
endif

return
end

