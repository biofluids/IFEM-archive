!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  L. Zhang
!  Northwestern University
!  This subroutine solves for the residual for all degrees of freedom
!____________________________________________________________________
!  L. Zhang, 06/24/2004
!  Tulane University
!  Revised the subroutine to array
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine blockcvoutput(xloc, dloc, doloc, p, hk, ien, f_fluids,ne_local,ien_local,node_local,nn_local, &
                         fden,fvis,I_fluid,rngface)
  use global_constants
  use run_variables
  use fluid_variables
  use solid_variables, only: nn_solid
  use r_common, only: density_solid, vis_solid
  implicit none

  integer ien(nen,ne)
  real* 8 xloc(nsd,nn)                      ! geometry(nodes' locations)
  real* 8 dloc(ndf,nn),doloc(ndf,nn)        ! unknowns(current/previous)
  real* 8 p(ndf,nn),hk(ne)                  ! p: residual; hk: element characteristic length

  real* 8 x(nsd,nen)
  real* 8 d(ndf,nen),d_old(ndf,nen)

  real* 8 eft0,det,effd,effm,effc
  real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
  real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real* 8 drt(ndf),drs(ndf)
  real* 8 dr(nsd,ndf)
  real* 8 u,v,w,pp,ug
  real* 8 tau(nsd,nsd)
  real* 8 hg,taum,tauc,vel,ree, taul
  real* 8 res_c,res_a(nsd),res_t(nsd)
  real* 8 prs_c,prs_t(nsd),p_vec(3),prs_cc(nsd)
  real* 8 mu,nu,ro,g(nsd)
  real* 8 tempc(ndf),temp
  real* 8 dtinv,oma,ama
  integer inl, ie, isd, iq, node,jsd
  integer ieface,irng !,inface

  real* 8 f_fluids(nsd,nn)
  real* 8 fnode(nsd,nen),fq(nsd)
  integer nn_local
  integer node_local(nn_local)
!-------------------------------------------
  real(8) fden(nn)
  real(8) local_den(nen)
  real(8) fvis(nn)
  real(8) local_vis(nen)
!---------------------------------------------
  real(8) TC,ZC,RC,P0,dens0
  real(8) I_fluid(nn)
  real(8) kappa
  integer rngface(neface,ne)
  integer flagincv(nen),nodecv(2)
  real(8) vol,MomTerm(nsd,4),EneTerm(4)
  real(8) normvec(nsd),edgelen
!============================
! MPI varibalbes
  integer ne_local ! # of element on each processor
  integer ien_local(ne_local) ! subregion-->wholeregion element index
  integer ie_local ! loop parameter
  integer icount
!---------------------------------------------
  TC=273.15
  ZC=1.4
  RC=2.87058e6
  P0=1.01325e6
  dens0=1.2922e-3
  kappa=1.0e4
  dtinv = 1.0/dt
  if(steady) dtinv = 0.0
  oma   = 1.0 - alpha
  ama   = 1.0 - oma
  MomTerm(1:nsd,1:4)=0.0
  EneTerm(1:4)=0.0
!===================================================
do ie=1,ne     ! loop over elements
    flagincv(1:nen)=0
    nodecv(:)=0
    do inl=1,nen
        node=ien(inl,ie)
        x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
        if (x(1,inl)>=12.0 .and. x(1,inl)<=15.0 .and. I_fluid(node)==0.0) then
            flagincv(inl)=1
        endif
        fnode(:,inl)=0.0
        d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
        d_old(1:ndf,inl) = doloc(1:ndf,ien(inl,ie))
!-----------------------------------------------------
        local_den(inl)=(1.0+dloc(ndf,node)/(ZC*P0))*dens0*(1.0-I_fluid(node))+&
                                (density_solid+den_liq)*I_fluid(node)
        local_vis(inl)=fvis(ien(inl,ie))
    enddo

    hg = hk(ie)
    vol=0.0
  if (sum(flagincv(1:nen)) .ge. 1) then
    do iq=1,nquad  ! loop over the quadrature points in each element 
!...  calculate the shape function and the weight at quad point
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
!	write(*,*) 'shape functions=',sh(0,1),sh(0,2),sh(0,3),sh(0,4)
        eft0 = abs(det) * wq(iq) ! calculate the weight at each quad pt
        vol=vol+eft0
        ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0

!...  initialize d, dd/dx, dd/dy, dd/dz, dd/dt
        drs(:) = 0.0           ! d
        drt(:) = 0.0           ! dd/dt
        dr(1:nsd,1:ndf)=0.0    ! dd/dxi
        fq(:)=0.0              ! node force

!... calculate vi, dvi/dxj
        do inl=1,nen
            tempc(1:nsd) = ama*d(1:nsd,inl)+oma*d_old(1:nsd,inl)
            drs(1:nsd) = drs(1:nsd)+sh(0,inl)*tempc(1:nsd)
                do isd=1,nsd
                    dr(isd,1:nsd) = dr(isd,1:nsd)+sh(isd,inl)*tempc(1:nsd)
                enddo
            fq(:) = fq(:) + sh(0,inl)*fnode(:,inl)        
        enddo

        ro=0.0
        mu=0.0
!... calculate dvi/dt, p, dp/dxi
        do inl=1,nen
            drt(1:nsd)=drt(1:nsd)+sh(0,inl)*(d(1:nsd,inl)-d_old(1:nsd,inl))*dtinv
            drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)
            dr(1:nsd,pdf)=dr(1:nsd,pdf)+sh(1:nsd,inl)*d(pdf,inl)      
!----------------------------------------------------------------------------------------                   
            ro=ro+sh(0,inl)*local_den(inl) 
            mu=mu+sh(0,inl)*local_vis(inl)
        enddo

!... define u=v1, v=v2, w=v3, pp=p
        if (nsd==2) then
            u = drs(udf)
            v = drs(vdf)
            w = 0.0
        elseif (nsd==3) then
            u = drs(udf)
            v = drs(vdf)
            w = drs(wdf)
        endif
        pp = drs(pdf)

      if (sum(flagincv(1:nen))==nen) then       ! only calculate volume integral when in CV
        do isd = 1, nsd
            if (nsd==2) then
                MomTerm(isd,1)=MomTerm(isd,1)+ro*drt(isd)*eft0
                MomTerm(isd,2)=MomTerm(isd,2)+ro*(u*dr(1,isd)+v*dr(2,isd))*eft0
                MomTerm(isd,3)=MomTerm(isd,3)+dr(isd,pdf)*eft0
            elseif (nsd==3) then
                write(*,*) 'under work'
            endif
        enddo

        EneTerm(1)=EneTerm(1)+ro*(drt(xsd)*u+drt(ysd)*v)*eft0
        EneTerm(2)=EneTerm(2)+ro*((u*dr(1,1)+v*dr(2,1))*u+(u*dr(1,2)+v*dr(2,2))*v)*eft0
      endif                        ! endif only volume integral in CV

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    enddo ! end of qudrature pts loop

      if ( (sum(flagincv(1:nen))>1 .and. sum(flagincv(1:nen))<nen) .or. sum(rngface(:,ie))/=0) then
!      if ( sum(flagincv(1:nen))==nen ) then
        do isd = 1,nsd
            do jsd = 1,nsd
                tau(isd,jsd) = mu*(dr(isd,jsd) + dr(jsd,isd))
            enddo
            tau(isd,isd)=tau(isd,isd)-mu*dr(isd,isd)*2.0/3.0
        enddo

            do irng=1,nen
                    if (nen==4) then
                        if (irng==1) then
                            nodecv(1)=1
                            nodecv(2)=2
                        elseif (irng==2) then
                            nodecv(1)=2
                            nodecv(2)=3
                        elseif (irng==3) then
                            nodecv(1)=3
                            nodecv(2)=4
                        elseif (irng==4) then
                            nodecv(1)=4
                            nodecv(2)=1
                        endif
                    elseif (nen==3) then
                        if (irng==1) then
                            nodecv(1)=1
                            nodecv(2)=2
                        elseif (irng==2) then
                            nodecv(1)=2
                            nodecv(2)=3
                        elseif (irng==3) then
                            nodecv(1)=3
                            nodecv(2)=1
                        endif
                    endif

                if (rngface(irng,ie)/=0) then
!                if (.true.) then
                    call normal2d_gen(x,nodecv(1),nodecv(2),normvec,nsd,nen,edgelen)
                    do isd=1,nsd
                        MomTerm(isd,4)=MomTerm(isd,4)+&
                                   (tau(isd,1)*normvec(1)+tau(isd,2)*normvec(2))*edgelen
                    enddo
                    EneTerm(3)=EneTerm(3)+(ph(xsd,inl)*u+ph(ysd,inl)*v)*pp
                    EneTerm(4)=EneTerm(4)+ph(1,inl)*(tau(1,xsd)*u+tau(1,ysd)*v) + ph(2,inl)*(tau(2,xsd)*u+tau(2,ysd)*v) &
                                      - mu*(ph(xsd,inl)*u+ph(ysd,inl)*v)*(dr(1,1)+dr(2,2))*2.0/3.0
                elseif (sum(flagincv(1:nen))<nen .and. flagincv(nodecv(1))==1 .and. flagincv(nodecv(2))==1) then
                    call normal2d_gen(x,nodecv(1),nodecv(2),normvec,nsd,nen,edgelen)
                    normvec(1:nsd)=-normvec(1:nsd)
                    do isd=1,nsd
                        MomTerm(isd,4)=MomTerm(isd,4)+&
                                   (tau(isd,1)*normvec(1)+tau(isd,2)*normvec(2))*edgelen
                    enddo
                    EneTerm(3)=EneTerm(3)+(ph(xsd,inl)*u+ph(ysd,inl)*v)*pp
                    EneTerm(4)=EneTerm(4)+ph(1,inl)*(tau(1,xsd)*u+tau(1,ysd)*v) + ph(2,inl)*(tau(2,xsd)*u+tau(2,ysd)*v) &
                                      - mu*(ph(xsd,inl)*u+ph(ysd,inl)*v)*(dr(1,1)+dr(2,2))*2.0/3.0
                endif
            enddo

      endif  ! end of if the element is on CV boundary

  endif
enddo ! end of element loop

5007 format(6e14.5)

open(unit=28, file = 'poa.momentumX',status='unknown',position='append')
write(28,5007) tt, MomTerm(1,1), MomTerm(1,2), MomTerm(1,3), MomTerm(1,4), &
                   MomTerm(1,1)+ MomTerm(1,2)+ MomTerm(1,3)- MomTerm(1,4)
close(28)


open(unit=29, file = 'poa.energy',status='unknown',position='append')
write(29,5007) tt, EneTerm(1), EneTerm(2), EneTerm(3), EneTerm(4), &
                   EneTerm(1)+ EneTerm(2)+ EneTerm(3)+ EneTerm(4)
close(29)


return
end subroutine blockcvoutput

