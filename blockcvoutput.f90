! Jack Yang, Rensselaer Polytechnic Institute
! June 2014
! this subroutine evaluates the control volume analysis terms
subroutine blockcvoutput(xloc, dloc, doloc, p, hk, ien, f_fluids,ne_local,ien_local,node_local,nn_local, &
                         fden,fvis,I_fluid,I_fluid_old,rngface, &
                         node_sbc,x_solid)
    use global_constants
    use run_variables
    use fluid_variables
    use solid_variables, only: nn_solid,nn_sbc
    use r_common, only: density_solid, vis_solid
    use lumpedmass
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

    real* 8 drt(ndf),drs(ndf),drtavg(ndf),drsavg(ndf)
    real* 8 dr(nsd,ndf),dravg(nsd,ndf),drro(nsd)
    real* 8 u,v,w,pp,ug,uavg,vavg,ppavg,muavg
    real* 8 tau(nsd,nsd)
    real* 8 hg,taum,tauc,vel,ree, taul
    real* 8 res_c,res_a(nsd),res_t(nsd)
    real* 8 prs_c,prs_t(nsd),p_vec(3),prs_cc(nsd)
    real* 8 mu,nu,ro,g(nsd)
    real* 8 tempc(ndf),temp
    real* 8 dtinv,oma,ama
    integer inl, ie, isd, iq, node,jsd, inn
    integer ieface,irng !,inface
    integer node1,node2
    integer node_sbc(nn_sbc)
    real* 8 x_solid(nsd,nn_solid)

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
    !real(8) TC,ZC,RC,P0,dens0
    real(8) I_fluid(nn)
    real(8) I_fluid_old(nn)
    real(8) kappa
    integer rngface(neface,ne)
    integer flagincv(nen),nodecv(2),flagnearimg(nn),flagcnln(nen),flagclus(nen),flagclds(nen)
    real(8) vol,MomTerm(nsd,5),EneTerm(6),BerUps(5),BerDwn(5)
    real(8) XBmax,XBmin,XImax,YImax,Imax
    real(8) normvec(nsd),edgelen,absedgelen
    real(8) hattauijcj(nsd),hattauijquad(nsd,nsd)
!============================
! MPI varibalbes
    integer ne_local ! # of element on each processor
    integer ien_local(ne_local) ! subregion-->wholeregion element index
    integer ie_local ! loop parameter
    integer icount
!---------------------------------------------
    !TC=273.15
    !ZC=1.4
    !RC=2.87058e6
    !P0=1.01325e6
    !dens0=1.2922e-3
    kappa=1.0e4
    dtinv = 1.0/dt
    if(steady) dtinv = 0.0
    oma   = 1.0 - alpha
    ama   = 1.0 - oma
    MomTerm(1:nsd,1:5)=0.0
    EneTerm(1:6)=0.0
    BerUps(1:5)=0.0
    BerDwn(1:5)=0.0
    XBmin=99999.0
    XBmax=-99999.0
    XImax=0.0
    YImax=-99999.0
    Imax=0.0

    include "reconstruct_tauij_global.f90"

!===================================================
    do inn=1,nn          ! loop thru fluid nodes and look at indicator as well as coordinates
        if (I_fluid(inn)>0.0) then
            if ((xloc(2,inn)>YImax) .and. (abs(xloc(2,inn)-YImax)>=5e-6)) then
                YImax = xloc(2,inn)
                XImax = xloc(1,inn)
                Imax = I_fluid(inn)
            elseif (abs(xloc(2,inn)-YImax)<5e-6) then
                !if (I_fluid(inn)>Imax) then
                if (xloc(1,inn) < XImax) then
                    YImax = xloc(2,inn)
                    XImax = xloc(1,inn)
                    Imax = I_fluid(inn)
                endif
            endif
    
            if (abs(xloc(2,inn)-1.397)<5e-6) then
                if (xloc(1,inn)>XBmax) XBmax=xloc(1,inn)
                if (xloc(1,inn)<XBmin) XBmin=xloc(1,inn)
            endif
        endif
    enddo
    XImax = 0.0
    YImax = -99999.0
    write(*,*) "before"
    do inn=1,nn_sbc      ! loop thru solid boundary nodes and look at coordinates
        node = node_sbc(inn)
        if ( x_solid(2,node)>YImax ) then
            YImax = x_solid(2,node)
            XImax = x_solid(1,node)
        endif
    enddo
    write(*,*) "after"
    if ((XBmax==-99999.0) .and. (XBmin==99999.0)) then
        XBmin = XImax
        XBmax = XImax
    endif
!===================================================
    flagnearimg(1:nn)=1
    do ie=1,ne
        flagincv(1:nen)=0
        do inl=1,nen
            node=ien(inl,ie)
            if (I_fluid(node)==0.0 .and. I_fluid_old(node)==0.0) then
                flagincv(inl)=1
            endif
        enddo
        if (sum(flagincv(1:nen))<nen) then
            do inl=1,nen
                node=ien(inl,ie)
                flagnearimg(node)=flagnearimg(node)*0
            enddo
        endif
    enddo
!===================================================
    do ie=1,ne     ! loop over elements

        !flagincv(1:nen) = 0
        !flagcnln(1:nen) = 0
        !do inl=1,nen
        !    node = ien(inl,ie)
        !    if ((abs(xloc(2,node)-1.397)<5.0e-6) .and. (xloc(1,node)>=12.0) .and. (xloc(1,node)<=15.0)) then
        !        flagcnln(inl) = 1
        !    endif
        !    if (I_fluid(node)==0.0 .and. I_fluid_old(node)==0.0) flagincv = 1
        !    d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
        !    d_old(1:ndf,inl) = doloc(1:ndf,ien(inl,ie))
        !    local_den(inl)=(1.0+dloc(ndf,node)/(ZC*P0))*dens0*(1.0-I_fluid(node))+&
        !                            (density_solid+den_liq)*I_fluid(node)
        !    local_vis(inl)=fvis(ien(inl,ie))
        !enddo
!
        !if (sum(flagcnln(1:nen))==2 .and. sum(flagincv(1:nen))==nen) then
        !    do irng=1,nen
        !        if (nen==4) then
        !            if (irng==1) then
        !                nodecv(1)=1
        !                nodecv(2)=2
        !            elseif (irng==2) then
        !                nodecv(1)=2
        !                nodecv(2)=3
        !            elseif (irng==3) then
        !                nodecv(1)=3
        !                nodecv(2)=4
        !            elseif (irng==4) then
        !                nodecv(1)=4
        !                nodecv(2)=1
        !            endif
        !        elseif (nen==3) then
        !            if (irng==1) then
        !                nodecv(1)=1
        !                nodecv(2)=2
        !            elseif (irng==2) then
        !                nodecv(1)=2
        !                nodecv(2)=3
        !            elseif (irng==3) then
        !                nodecv(1)=3
        !                nodecv(2)=1
        !            endif
        !        endif
!
        !        node1 = ien(nodecv(1),ie)
        !        node2 = ien(nodecv(2),ie)
        !        edgelen = xloc(1,node1)-xloc(1,node2)
        !        absedgelen = abs(edgelen)
        !        if ( (flagcnln(nodecv(1))==1) .and. (flagcnln(nodecv(2))==1) .and. &
        !            (xloc(1,node1)>XBmax) .and. (xloc(1,node2)>XBmax) ) then
        !            if (xloc(1,node1)>xloc(1,node2)) then
        !                BerDwn(1) = BerDwn(1) + 0.5*d(udf,nodecv(1))**2 - 0.5*d(udf,nodecv(2))**2
        !                BerDwn(2) = BerDwn(2) + d(pdf,nodecv(1))/local_den(nodecv(1)) - &
        !                              d(pdf,nodecv(2))/local_den(nodecv(2))
        !            elseif (xloc(1,node1)<xloc(1,node2)) then
        !                BerDwn(1) = BerDwn(1) - 0.5*d(udf,nodecv(1))**2 + 0.5*d(udf,nodecv(2))**2
        !                BerDwn(2) = BerDwn(2) - d(pdf,nodecv(1))/local_den(nodecv(1)) + &
        !                              d(pdf,nodecv(2))/local_den(nodecv(2))
        !            endif
        !            BerDwn(3) = BerDwn(3) + (d(udf,nodecv(1))-d_old(udf,nodecv(1))+&
        !                                     d(udf,nodecv(2))-d_old(udf,nodecv(2)))*0.5*dtinv*absedgelen
        !            BerDwn(4) = BerDwn(4) + ((d(pdf,nodecv(1))+d(pdf,nodecv(2)))*0.5)/&
        !                                    (((local_den(nodecv(1))+local_den(nodecv(2)))*0.5)**2)*&
        !                                    (local_den(nodecv(1))-local_den(nodecv(2)))/edgelen*absedgelen
        !            BerDwn(5) = BerDwn(5) + 2.0/(local_den(nodecv(1))+local_den(nodecv(2)))*&
        !                                    (hattauij(1,1,node1)-hattauij(1,1,node2))/edgelen*absedgelen
!
        !        elseif ( (flagcnln(nodecv(1))==1) .and. (flagcnln(nodecv(2))==1) .and. &
        !                (xloc(1,node1)<XBmin) .and. (xloc(1,node2)<XBmin) ) then
        !            if (xloc(1,node1)>xloc(1,node2)) then
        !                BerUps(1) = BerUps(1) + 0.5*d(udf,nodecv(1))**2 -0.5*d(udf,nodecv(2))**2
        !                BerUps(2) = BerUps(2) + d(pdf,nodecv(1))/local_den(nodecv(1))- &
        !                              d(pdf,nodecv(2))/local_den(nodecv(2))
        !            elseif (xloc(1,node1)<xloc(1,node2)) then
        !                BerUps(1) = BerUps(1) - 0.5*d(udf,nodecv(1))**2 +0.5*d(udf,nodecv(2))**2
        !                BerUps(2) = BerUps(2) - d(pdf,nodecv(1))/local_den(nodecv(1))+ &
        !                              d(pdf,nodecv(2))/local_den(nodecv(2))
        !            endif
        !            BerUps(3) = BerUps(3) + (d(udf,nodecv(1))-d_old(udf,nodecv(1))+&
        !                                     d(udf,nodecv(2))-d_old(udf,nodecv(2)))*0.5*dtinv*absedgelen
        !            BerUps(4) = BerUps(4) + ((d(pdf,nodecv(1))+d(pdf,nodecv(2)))*0.5)/&
        !                                    (((local_den(nodecv(1))+local_den(nodecv(2)))*0.5)**2)*&
        !                                    (local_den(nodecv(1))-local_den(nodecv(2)))/edgelen*absedgelen
        !            BerUps(5) = BerUps(5) + 2.0/(local_den(nodecv(1))+local_den(nodecv(2)))*&
        !                                    (hattauij(1,1,node1)-hattauij(1,1,node2))/edgelen*absedgelen
        !        endif
        !    enddo
        !endif

        flagincv(1:nen)=0
        flagclus(1:nen)=0
        flagclds(1:nen)=0
        nodecv(:)=0
        do inl=1,nen
            node=ien(inl,ie)
            x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
            if (x(1,inl)>=12.0 .and. x(1,inl)<=15.0 .and. flagnearimg(node)==1) then
                      ! .and. I_fluid(node)==0.0 .and. I_fluid_old(node)==0.0) then
                flagincv(inl)=1
            endif
            if ((x(1,inl)>=12.7) .and. (x(1,inl)<=XBmin) .and. (x(2,inl)>=1.386) .and. (flagnearimg(node)==1)) then
                flagclus(inl) = 1
            endif
            if ((x(1,inl)<=15.0) .and. (x(1,inl)>=XBmax) .and. (x(2,inl)>=1.386) .and. (flagnearimg(node)==1)) then
                flagclds(inl) = 1
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
        drtavg(:)=0.0
        drsavg(:)=0.0
        dravg(:,:)=0.0
        uavg=0.0
        vavg=0.0
        ppavg=0.0
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
                drro(1:nsd)=0.0
!... calculate dvi/dt, p, dp/dxi
                do inl=1,nen
                    drt(1:ndf)=drt(1:ndf)+sh(0,inl)*(d(1:ndf,inl)-d_old(1:ndf,inl))*dtinv
                    drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)
                    dr(1:nsd,pdf)=dr(1:nsd,pdf)+sh(1:nsd,inl)*d(pdf,inl)      
                    ro=ro+sh(0,inl)*local_den(inl)
                    drro(1:nsd) = drro(1:nsd)+sh(1:nsd,inl)*local_den(inl)
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

!----------------------------------------------------------------------------------------                   
                drtavg(1:ndf)=drtavg(1:ndf)+drt(1:ndf)
                drsavg(1:ndf)=drsavg(1:ndf)+drs(1:ndf)
                dravg(1:nsd,1:ndf)=dravg(1:nsd,1:ndf)+dr(1:nsd,1:ndf)
                uavg=uavg+u
                vavg=vavg+v
                ppavg=ppavg+pp
                muavg=muavg+mu
!---------------------------------------------------------------------------------------- 
                hattauijcj(1:nsd)=0.0
                hattauijquad(1:nsd,1:nsd)=0.0

                do inl=1,nen
                    node=ien(inl,ie)
                    do isd=1,nsd
                        do jsd=1,nsd
                            hattauijcj(isd) = hattauijcj(isd)+&
                                            sh(jsd,inl)*hattauij(isd,jsd,node)
                            hattauijquad(isd,jsd) = hattauijquad(isd,jsd)+&
                                                  sh(0,inl)*hattauij(isd,jsd,node)
                        enddo
                    enddo
                enddo
!----------------------------------------------------------------------------------------

                if (sum(flagincv(1:nen))==nen) then       ! only calculate volume integral when in CV
                    do isd = 1, nsd
                        if (nsd==2) then
                            MomTerm(isd,1)=MomTerm(isd,1)+ro*drt(isd)*eft0
                            MomTerm(isd,2)=MomTerm(isd,2)+ro*(u*dr(1,isd)+v*dr(2,isd))*eft0
                            MomTerm(isd,3)=MomTerm(isd,3)+dr(isd,pdf)*eft0
                            MomTerm(isd,5)=MomTerm(isd,5)+hattauijcj(isd)*eft0
                        elseif (nsd==3) then
                            write(*,*) 'under work'
                        endif
                    enddo
  
                    if (nsd==2) then
                        EneTerm(1)=EneTerm(1)+ro*(drt(xsd)*u+drt(ysd)*v)*eft0
                        EneTerm(2)=EneTerm(2)+ro*((u*dr(1,1)+v*dr(2,1))*u+(u*dr(1,2)+v*dr(2,2))*v)*eft0
                        EneTerm(3)=EneTerm(3)+(dr(1,pdf)*u+dr(2,pdf)*v)*eft0
                        !EneTerm(5)=EneTerm(5)+mu*(0.5*((dr(1,1)+dr(1,1))**2+(dr(1,2)+dr(2,1))**2+(dr(2,1)+dr(1,2))**2+(dr(2,2)+dr(2,2))**2)- &
                        !                           (dr(1,1)+dr(2,2))**2*2.0/3.0)*eft0
                        EneTerm(5)=EneTerm(5)+(hattauijquad(1,1)*dr(1,1)+hattauijquad(1,2)*dr(1,2)+&
                                               hattauijquad(2,1)*dr(2,1)+hattauijquad(2,2)*dr(2,2))*eft0
                        EneTerm(6)=EneTerm(6)+(hattauijcj(1)*u+hattauijcj(2)*v)*eft0
                    elseif (nsd==3) then
                        write(*,*) 'under work'
                    endif
                endif                        ! endif only volume integral in CV

                if (sum(flagclus(1:nen))==nen) then
                    if (nsd==2) then
                        BerUps(1) = BerUps(1) + (u*dr(1,1)+v*dr(2,1))*eft0
                        !BerUps(2) = BerUps(2) + 1.09778e9*(1/(1418551.316+pp)*dr(1,pdf)-pp/(1418551.316+pp)**2*dr(1,pdf))*eft0
                        BerUps(2) = BerUps(2) + dr(1,pdf)*dens0*eft0/(ro**2)
                        BerUps(3) = BerUps(3) + drt(1)*eft0
                        BerUps(4) = BerUps(4) + pp*drro(1)*eft0/(ro**2)
                        BerUps(5) = BerUps(5) + hattauijcj(1)*eft0/ro
                    elseif (nsd==3) then
                        write(*,*) 'under work'
                    endif
                endif
                if (sum(flagclds(1:nen))==nen) then
                    if (nsd==2) then
                        BerDwn(1) = BerDwn(1) + (u*dr(1,1)+v*dr(2,1))*eft0
                        !BerDwn(2) = BerDwn(2) + 1.09778e9*(1/(1418551.316+pp)*dr(1,pdf)-pp/(1418551.316+pp)**2*dr(1,pdf))*eft0
                        BerDwn(2) = BerDwn(2) + dr(1,pdf)*dens0*eft0/(ro**2)
                        BerDwn(3) = BerDwn(3) + drt(1)*eft0
                        BerDwn(4) = BerDwn(4) + pp*drro(1)*eft0/(ro**2)
                        BerDwn(5) = BerDwn(5) + hattauijcj(1)*eft0/ro
                    elseif (nsd==3) then
                        write(*,*) 'under work'
                    endif
                endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            enddo ! end of qudrature pts loop

!----------------------------------------------------------------------------------------                   
            drtavg(1:ndf)=drtavg(1:ndf)/nquad
            drsavg(1:ndf)=drsavg(1:ndf)/nquad
            dravg(1:nsd,1:ndf)=dravg(1:nsd,1:ndf)/nquad
            uavg=uavg/nquad
            vavg=vavg/nquad
            ppavg=ppavg/nquad
            muavg=1.8e-4
!----------------------------------------------------------------------------------------                   
            if ((sum(flagincv(1:nen))>1 .and. sum(flagincv(1:nen))<nen) .or. sum(rngface(:,ie))/=0 ) then
!      if ( sum(flagincv(1:nen))==nen ) then
                do isd = 1,nsd
                    do jsd = 1,nsd
                        tau(isd,jsd) = muavg*(dravg(isd,jsd) + dravg(jsd,isd))
                    enddo
                    tau(isd,isd)=tau(isd,isd)-muavg*dravg(isd,isd)*2.0/3.0
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

                    if (sum(flagincv(1:nen))==nen .and. rngface(irng,ie)/=0) then
!                    if (.true.) then
                        call normal2d_gen(x,nodecv(1),nodecv(2),normvec,nsd,nen,edgelen)
                        do isd=1,nsd
                            MomTerm(isd,4)=MomTerm(isd,4)+&
                                       (tau(isd,1)*normvec(1)+tau(isd,2)*normvec(2))*edgelen
                        enddo
                        EneTerm(4)=EneTerm(4)+ (tau(1,1)*normvec(1)+tau(1,2)*normvec(2))*edgelen*uavg +&
                                               (tau(2,1)*normvec(1)+tau(2,2)*normvec(2))*edgelen*vavg
                    elseif (sum(flagincv(1:nen))<nen .and. flagincv(nodecv(1))==1 .and. flagincv(nodecv(2))==1) then
                        call normal2d_gen(x,nodecv(1),nodecv(2),normvec,nsd,nen,edgelen)
                        normvec(1:nsd)=-normvec(1:nsd)
                        do isd=1,nsd
                            MomTerm(isd,4)=MomTerm(isd,4)+&
                                       (tau(isd,1)*normvec(1)+tau(isd,2)*normvec(2))*edgelen
                        enddo
                        EneTerm(4)=EneTerm(4)+ (tau(1,1)*normvec(1)+tau(1,2)*normvec(2))*edgelen*uavg +&
                                               (tau(2,1)*normvec(1)+tau(2,2)*normvec(2))*edgelen*vavg
                    endif
                enddo

            endif  ! end of if the element is on CV boundary

        endif
    enddo ! end of element loop

    5007 format(7e14.5)
    5008 format(8e14.5)
    5009 format(9e14.5)

    open(unit=28, file = 'poa.momentumX',status='unknown',position='append')
    write(28,5007) tt, MomTerm(1,1), MomTerm(1,2), MomTerm(1,3), MomTerm(1,4), &
                       MomTerm(1,1)+ MomTerm(1,2)+ MomTerm(1,3)- MomTerm(1,4), &
                       MomTerm(1,5)
    close(28)

    open(unit=29, file = 'poa.energy',status='unknown',position='append')
    write(29,5008) tt, EneTerm(1), EneTerm(2), EneTerm(3), EneTerm(4), EneTerm(5), &
                       EneTerm(1)+ EneTerm(2)+ EneTerm(3)- EneTerm(4)+ EneTerm(5), &
                       EneTerm(6)
    close(29)

    open(unit=30, file = 'poa.bernoulliDwn', status='unknown', position='append')
    write(30,5009) tt, BerDwn(1), BerDwn(2), BerDwn(3), BerDwn(4), BerDwn(5), &
                       BerDwn(1)+BerDwn(2)+BerDwn(3)+BerDwn(4)-BerDwn(5), XBmin, XBmax
    close(30)

    open(unit=31, file = 'poa.bernoulliUps', status='unknown', position='append')
    write(31,5007) tt, BerUps(1), BerUps(2), BerUps(3), BerUps(4), BerUps(5), &
                       BerUps(1)+BerUps(2)+BerUps(3)+BerUps(4)-BerUps(5)
    close(31)

return
end subroutine blockcvoutput
