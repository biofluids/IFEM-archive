!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  L. Zhang
!  Northwestern University
!  This subroutine solves for the residual for all degrees of freedom
!____________________________________________________________________
!  L. Zhang, 06/24/2004
!  Tulane University
!  Revised the subroutine to array
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine block(xloc, dloc, doloc, p, q_p, hk, ien, f_fluids,rngface, f_stress,&
                 ne_local,ien_local,node_local,nn_local,fden,fvis,I_fluid)
    use global_constants
    use run_variables
    use fluid_variables
    use solid_variables, only: nn_solid
    use r_common, only: density_solid, vis_solid
    use mpi_variables
    use pml_variables
    use lumpedmass
    implicit none

    integer ien(nen,ne)
    real* 8 xloc(nsd,nn)
    real* 8 dloc(ndf,nn),doloc(ndf,nn)
    real* 8 p(ndf,nn+nn_PML),q(ndf,nn),hk(ne)
    integer flagpmlelm
    real(8) betaPML

    real* 8 x(nsd,nen)
    real* 8 d(ndf,nen),d_old(ndf,nen)
    real* 8 qvel(ndf,nen),qvoldel(ndf,nen)

    real* 8 eft0,det,effd,effm,effc
    real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
    real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
    real* 8 f_stress(nsd,nsd,nn)

    real* 8 qrs(ndf), qrt(ndf)  ! interpolation and derivatives of qv
    real* 8 drt(ndf),drs(ndf)
    real* 8 dr(nsd,ndf)
    real* 8 u,v,w,pp,ug
    real* 8 tau(nsd,nsd)
    real* 8 hg,taum,tauc,vel,ree, taul
    real* 8 res_c,res_a(nsd),res_t(nsd)
    real* 8 prs_c,prs_t(nsd),p_vec(3),prs_cc(nsd)
    real* 8 mu,nu,ro,g(nsd),sigmaPe(nsd)
    real* 8 tempc(ndf),temp
    real* 8 dtinv,oma,ama
    integer inl, ie, isd, iq, node,jsd
    integer ieface,irng, rngface(neface,ne)
    real* 8 res_pml_c, res_pml_a(nsd), res_pml_t(nsd)

    real* 8 f_fluids(nsd,nn)
    real* 8 I_fluid(nn)
    real* 8 fnode(nsd,nen),fq(nsd)
    integer nn_local
    integer node_local(nn_local)
    !======================================
    ! Defined by Chu
    ! updated by Jack, 05/29/2014, for PML
    real* 8 q_d(ndf,nen)
    real* 8 q_res_c(nen), q_pml_c(nen)
    real* 8 q_p(ndf,nn+nn_PML)
    real* 8 q_res_a(nsd,nen), q_pml_a(nsd,nen)
    real* 8 diag(12)
    !======================================
    ! Define by Xingshi
    ! MPI varialbes & implicit FSI force
    real(8) fden(nn)
    real(8) local_den(nen)
    real(8) fvis(nn)
    real(8) local_vis(nen)
    !======================================
    ! varibles for mpi implementation
    integer ne_local ! # of element on each processor
    integer ien_local(ne_local) ! subregion-->whole region element index
    integer ie_local ! loop parameter
    integer icount
    !--------------------------------------------------
    !  real*8 TC,ZC,RC,P0,dens0          ! defined in global_variables now
    real*8 kappa !compressibility for artificial fluids
    !==================================================
    q_res_a(1:nsd,1:nen) = 0
    q_p(1:ndf,1:nn+nn_PML) = 0
    q_d(1:ndf,1:nen) = 1 !set each ndf for each node as 1
    q_res_c(1:nen) = 0
    !---------------------------------------------------
    !  TC=273.15
    !  ZC=1.4
    !  RC=2.87058e6
    !  P0=1.01325e6
    !  dens0=1.2922e-3
    !................................
    ! Jack 05/20/2014
    ! Defined in global_variables now
    !---------------------------------------------------
    kappa=1.0e4
    !---------------------------------------------------
    dtinv = 1.0/dt
    if(steady) dtinv = 0.0
    oma   = 1.0 - alpha
    ama   = 1.0 - oma
    if (nn_PML_local > 0) betaPML=( uvwpAvgPML(udf)/sqrt(ZC*P0/dens0) )/( 1 - uvwpAvgPML(udf)**2/(ZC*P0/dens0) )
    !betaPML=0.0d0
    !=================================================
    do icount=1, nn_local
        node=node_local(icount)
        p(1:nsd,node)=p(1:nsd,node)+f_fluids(1:nsd,node)
    enddo

    include "reconstruct_tauij.f90"
    !===================================================
    ! f_fluids is actually the FSI force at fluid nodes,
    ! then we will just subscrib it from p(!:nsd) which is the 
    ! residuals for the momentum equations. 
    ! 1 let fnod=0 then fndoe and fq 's effects will be canceled out
    ! 2 do the subscribition after the elements loop
    ! Xingshi 09/15/2008
    !===================================================
    do ie_local=1,ne_local              ! loop over subregion elements
        ie=ien_local(ie_local)
        flagpmlelm=0
        do inl=1,nen
            x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
            fnode(:,inl) = 0.0
            d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
            d_old(1:ndf,inl) = doloc(1:ndf,ien(inl,ie))
            f_stress(1:nsd,1:nsd,ien(inl,ie)) = 0.0
            !----------------------------------------------
            node=ien(inl,ie)
            local_den(inl)=(dloc(ndf,node)*dens0/(ZC*P0)+dens0)*(1.0-I_fluid(node))+&
                           (density_solid+den_liq)*I_fluid(node)
            local_vis(inl)=fvis(ien(inl,ie))
            !------------------- check PML nodes ------------------------------
            if (seqcPML(node) .ne. 0) flagpmlelm=flagpmlelm+1
            !------------------------------------------------------------------
        enddo
        hg = hk(ie)
        !======================= LOOP THRU QUADRATURE POINTS ==========================
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
            eft0 = abs(det) * wq(iq) ! calculate the weight at each quad pt
            ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0      
            !...  initialize d, dd/dx, dd/dy, dd/dz, dd/dt
            drs(:) = 0.0
            drt(:) = 0.0
            dr(1:nsd,1:ndf)=0.0
            !-----PML-------------
            qrs(1:ndf) = 0.0
            qrt(1:ndf) = 0.0
            !---------------------
            fq(:)=0.0
            !... calculate vi, dvi/dxj
            do inl=1,nen
                node=ien(inl,ie)
                tempc(1:nsd) = ama*d(1:nsd,inl)+oma*d_old(1:nsd,inl)
                drs(1:nsd) = drs(1:nsd)+sh(0,inl)*tempc(1:nsd)
                do isd=1,nsd
                    dr(isd,1:nsd) = dr(isd,1:nsd)+sh(isd,inl)*tempc(1:nsd)
                enddo
                fq(:) = fq(:) + sh(0,inl)*fnode(:,inl)        
            enddo

            ro=0.0
            mu=0.0
            sigmaPe(1:nsd)=0.0
            !... calculate dvi/dt, p, dp/dxi
            do inl=1,nen
                node=ien(inl,ie)
                drt(1:ndf)=drt(1:ndf)+sh(0,inl)*(d(1:ndf,inl)-d_old(1:ndf,inl))*dtinv
                !------PML---------
                if (seqcPML(node) > 0) then
                    qrs(1:ndf) = qrs(1:ndf)+sh(0,inl)*qv(1:ndf,node)
                    qrt(1:ndf) = qrt(1:ndf)+sh(0,inl)*(qv(1:ndf,node)-qvold(1:ndf,node))*dtinv
                endif
                sigmaPe(1:nsd)=sigmaPe(1:nsd)+sh(0,inl)*sigmaPML(1:nsd,node)
                !------------------
                drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)
                dr(1:nsd,pdf)=dr(1:nsd,pdf)+sh(1:nsd,inl)*d(pdf,inl) 
                !----------------------------------------
                ro=ro+sh(0,inl)*local_den(inl)  ! local fluid density
                mu=mu+sh(0,inl)*local_vis(inl)  ! local fluid viscosity
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
            pp= drs(pdf)

            if(stokes) then ! if stokes flow
                u = 0.0
                v = 0.0
                w = 0.0
            endif
            !....  calculate liquid constant and gravity
            g(1:nsd)  = gravity(1:nsd)  ! gravatitional force

            ! believe nu is calculated only for turbulent model
            if (nsd==2) then
                nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*dr(1,1)**2+(dr(2,1)+dr(1,2))**2 &
                                                        +2*dr(2,2)**2)
            elseif (nsd==3) then
                nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*dr(1,1)**2+(dr(2,1)+dr(1,2))**2 &
                                                        +2*dr(2,2)**2+(dr(3,1)+dr(1,3))**2 &
                                                        +2*dr(3,3)**2+(dr(3,2)+dr(2,3))**2)
            endif
            mu = mu + nu*ro                   

            !....  calculate each term in the residual equation
            res_c = 0.0                ! Continuity residual
            res_pml_c = 0.0
            res_pml_a(1:nsd) = 0.0
            q_pml_c(1:nen) = 0.0
            q_pml_a(1:nsd,1:nen) = 0.0
            if (nsd==2) then
                do inl=1,nen
                    node=ien(inl,ie)
                    res_c = res_c + (ZC+pp/P0)*(sh(xsd,inl)*d(udf,inl)+ sh(ysd,inl)*d(vdf,inl))+ &
                            (u*sh(xsd,inl)*d(pdf,inl)+v*sh(ysd,inl)*d(pdf,inl))/P0*(1.0-I_fluid(node))
                    q_res_c(inl) = (sh(0,inl)*dtinv+(ZC+pp/P0)*(u*sh(xsd,inl)+v*sh(ysd,inl)))/P0*(1.0-I_fluid(node)) &
                                 +sh(0,inl)*dtinv/kappa*I_fluid(node)
                enddo
            elseif (nsd==3) then
                do inl=1,nen
                    node = ien(inl,ie)
                    res_c = res_c + ZC*(sh(xsd,inl)*d(udf,inl) &
                                  + sh(ysd,inl)*d(vdf,inl) &
                                  + sh(zsd,inl)*d(wdf,inl)) &
                                  + (u*sh(xsd,inl)*d(pdf,inl)+v*sh(ysd,inl)*d(pdf,inl)+w*sh(zsd,inl)*d(pdf,inl))&
                                  / P0*(1.0-I_fluid(node))
                    q_res_c(inl) = (sh(0,inl)*dtinv+ZC*(u*sh(xsd,inl)+v*sh(ysd,inl)+w*sh(zsd,inl)) )/P0*(1.0-I_fluid(node)) &
                                   + sh(0,inl)*dtinv/kappa*I_fluid(node)
                enddo
            endif

            do inl=1,nen
                node=ien(inl,ie)
                res_c=res_c+sh(0,inl)*( d(ndf,inl)-d_old(ndf,inl) )*dtinv* &
                      (1.0/P0*(1.0-I_fluid(node))+1.0/kappa*I_fluid(node))
            enddo  ! add dp/dt term

            do isd = 1, nsd
                if (nsd==2) then
                    res_a(isd)=ro*( drt(isd) + u*dr(1,isd) + v*dr(2,isd) &
                                  + drs(isd)*(dr(1,1) + dr(2,2)) - g(isd) ) &
                               + ( drs(isd)*(u*dr(1,pdf) + v*dr(2,pdf)) + drt(pdf)*drs(isd) )*dens0/(ZC*P0) &
                               - fq(isd)
                    do inl = 1, nen
                        q_res_a(isd,inl)=ro*(sh(0,inl)*q_d(isd,inl)*dtinv+u*sh(1,inl)*q_d(isd,inl)+v*sh(2,inl)*q_d(isd,inl) &
                                             + q_d(isd,inl)*(dr(1,1) + dr(2,2))) &
                                         + q_d(isd,inl)*(u*dr(1,pdf) + v*dr(2,pdf) + drt(pdf))*dens0/(ZC*P0)
                    enddo ! get res_a for u and v for momentum equation
                elseif (nsd==3) then
                    res_a(isd)=ro*(drt(isd)+u*dr(1,isd)+v*dr(2,isd)+w*dr(3,isd)-g(isd))-fq(isd)
                    do inl = 1, nen
                        q_res_a(isd,inl)=ro*(sh(0,inl)*dtinv+u*sh(1,inl)+v*sh(2,inl)+w*sh(3,inl))
                    enddo ! get res_a for u v and w for momentum equation
                endif
            enddo

            !c....   Calculate stress tau(ij) and pressure
            do isd = 1,nsd
                do jsd = 1,nsd
                    tau(isd,jsd) = mu*(dr(isd,jsd) + dr(jsd,isd))
                    do inl=1,nen
                        ! Compute fluid stress
                        f_stress(isd,jsd,ien(inl,ie))= f_stress(isd,jsd,ien(inl,ie)) + tau(isd,jsd)
                    enddo
                enddo
                if (nsd==2) then
                    tau(isd,isd)=tau(isd,isd)-2.0*mu*(dr(1,1)+dr(2,2))/3.0
                elseif (nsd==3) then
                    tau(isd,isd)=tau(isd,isd)-2.0*mu*(dr(1,1)+dr(2,2)+dr(3,3))/3.0
                endif
            enddo

!--------------- PML terms in Continuity/Momentum eqns --------------------------
            if (flagpmlelm > 0) then
                !------------ residual correction for NS eqns ----------------------
                if (nsd==2) then
                    do inl=1,nen
                        node=ien(inl,ie)
                        !q_pml_c(inl)=(sigmaPe(xsd)+sigmaPe(ysd))/P0
                        q_pml_a(xsd,inl)=(sigmaPe(xsd)+sigmaPe(ysd))*ro*sh(0,inl) + betaPML*sigmaPe(xsd)*ro*u*sh(0,inl)
                        q_pml_a(ysd,inl)=(sigmaPe(xsd)+sigmaPe(ysd))*ro*sh(0,inl) + betaPML*sigmaPe(xsd)*ro*u*sh(0,inl)
                    enddo
                    res_pml_c=( sigmaPe(xsd)*qrs(pdf) + sigmaPe(ysd)*(pp-qrs(pdf)) )/P0 + &
                                betaPML*sigmaPe(xsd)*((ZC+pp/P0)*u-ZC*uvwpAvgPML(udf))
                    res_pml_a(xsd)=ro*( sigmaPe(xsd)*qrs(udf)+sigmaPe(ysd)*(u-qrs(udf)) ) + betaPML*sigmaPe(xsd)*ro*u*u
                    res_pml_a(ysd)=ro*( sigmaPe(xsd)*qrs(vdf)+sigmaPe(ysd)*(v-qrs(vdf)) ) + betaPML*sigmaPe(xsd)*ro*u*v
                endif
                res_c=res_c+res_pml_c
                res_a(1:nsd)=res_a(1:nsd)+res_pml_a(1:nsd)
                !q_res_c(1:nen)=q_res_c(1:nen)+q_pml_c(1:nen)
                q_res_a(1:nsd,1:nen)=q_res_a(1:nsd,1:nen)+q_pml_a(1:nsd,1:nen)
                !------------ residual for auxiliary variables ----------------------
                ! res_pml_a & res_pml_c have been recycled to calculate dq/dt+\sigma*q+dF1/dx
                res_pml_c=0.0
                res_pml_a(:)=0.0
                q_pml_a(:,:)=0.0
                q_pml_c(:)=0.0
                if (nsd==2) then
                    ! PML Continuity
                    res_pml_c = ( qrt(pdf) + sigmaPe(xsd)*qrs(pdf) + (ZC*P0+pp)*dr(xsd,udf) + u*dr(xsd,pdf) )/P0 + &
                                betaPML*sigmaPe(xsd)*( (ZC+pp/P0)*u - ZC*uvwpAvgPML(udf) )
                    do inl=1,nen
                        node=ien(inl,ie)
                        q_pml_c(inl) = ( sh(0,inl)*dtinv + (ZC+pp/P0)*(u*sh(xsd,inl)) )/P0 + &
                                       betaPML*sigmaPe(xsd)*ZC*sh(0,inl)
                    enddo
                    ! PML Momentum
                    do isd = 1, nsd
                        res_pml_a(isd) = ro*( qrt(isd) + sigmaPe(xsd)*qrs(isd) + u*dr(1,isd) + drs(isd)*dr(1,1) ) + &
                                         ( drs(isd)*u*dr(1,pdf) + qrt(pdf)*drs(isd) )*dens0/(ZC*P0) + &
                                         betaPML*sigmaPe(xsd)*ro*u*drs(isd)
                        do inl = 1, nen
                            q_pml_a(isd,inl)=ro*( sh(0,inl)*q_d(isd,inl)*dtinv+u*sh(1,inl)*q_d(isd,inl) + q_d(isd,inl)*dr(1,1) ) &
                                             + q_d(isd,inl)*(u*dr(1,pdf) + qrt(pdf))*dens0/(ZC*P0) &
                                             + betaPML*sigmaPe(xsd)*ro*u*sh(0,inl)
                        enddo
                    enddo
                    res_pml_a(xsd) = res_pml_a(xsd) + betaPML*sigmaPe(xsd)*pp
                endif
            
                do inl=1,nen
                    node = ien(inl,ie)
                    if (seqcPML(node) > 0) then
                        p(pdf,nn+seqcPML(node)) = p(pdf,nn+seqcPML(node)) - ph(0,inl)*res_pml_c
                        q_p(pdf,nn+seqcPML(node)) = q_p(pdf,nn+seqcPML(node)) - ph(0,inl)*q_pml_c(inl)
                        p(1:nsd,nn+seqcPML(node)) = p(1:nsd,nn+seqcPML(node)) - ph(0,inl)*res_pml_a(1:nsd)
                        q_p(1:nsd,nn+seqcPML(node)) = q_p(1:nsd,nn+seqcPML(node)) - ph(0,inl)*q_pml_a(1:nsd,inl)
                        p(udf,nn+seqcPML(node)) = p(udf,nn+seqcPML(node)) + ph(xsd,inl)*(pp-tau(1,1))
                        p(vdf,nn+seqcPML(node)) = p(vdf,nn+seqcPML(node)) - ph(xsd,inl)*tau(2,1)
                        q_p(udf,nn+seqcPML(node))=q_p(1,nn+seqcPML(node)) + ph(xsd,inl)*mu*(sh(1,inl)*q_d(1,inl)*(2.0-2.0/3.0))
                        q_p(vdf,nn+seqcPML(node))=q_p(2,nn+seqcPML(node)) + ph(xsd,inl)*mu*(sh(1,inl)+sh(2,inl))*q_d(2,inl)
                    endif
                enddo
                !--------------------------------------------------------------------
            endif
!---------------------------------------------------------------------------------------
            res_t(1:nsd) = dr(1:nsd,pdf) + res_a(1:nsd)

            include "add_tauijcj_resmomentum.f90"

            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ! Mickael 02/01/2005
            ! TAUm, TAUc and TAUl (l for  lsic), See Tezduyar and Sathe, 2003
            vel  = sqrt(u*u+v*v+w*w)     ! magnitude of the velocity
            ree  = vel*hg*ro/(2.0*mu)    ! Reynolds number
            if(steady.or.(.not.taudt)) then !stablization, taum
                taum = 1.0/sqrt((2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
            else
                taum = 1.0/sqrt((2.0/dt)**2+(2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
            endif
            tauc = taum/ro

            if (ree.le.3.0) then
                taul = hg*vel*ree/6.0
            else
                taul = hg*vel/2.0
            endif

            prs_t(1:nsd) = res_t(1:nsd)*taum
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ! Mickael 02/01/2005
            ! TAUl (l for  lsic), See Tezduyar and Sathe, 2003
            prs_c = ro*res_c*taul
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            prs_cc(1:nsd) = res_t(1:nsd)*tauc 
            !.... calculate the residual at each degree of freedom
            do inl=1,nen ! loop over number of nodes in an element
                node = ien(inl,ie)
                if (nsd==2) then
                    temp = ro*(u*ph(xsd,inl)+v*ph(ysd,inl))
                elseif (nsd==3) then
                    temp = ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))
                endif
!=============================================================
                ! Continuty Equation
                p(pdf,node) = p(pdf,node)-ph(0,inl)*res_c
                q_p(pdf,node) = q_p(pdf,node)+ph(0,inl)*q_res_c(inl) ! get diag for continuity
                ! Momentum Equation (Euler Residual)
                ! why originally it is minus?
                p(1:nsd,node) = p(1:nsd,node)-ph(0,inl)*res_a(1:nsd)
                q_p(1:nsd,node) = q_p(1:nsd,node)+ph(0,inl)*q_res_a(1:nsd,inl)
!=============================================================
                ! Momentum Equation (C: Diffusion terms)
                if (nsd==2) then
                    do isd=1,nsd
                        p(isd,node)=p(isd,node) + ph(isd,inl)*pp -   &
                                                  ph(1,inl)*tau(1,isd) -  &
                                                  ph(2,inl)*tau(2,isd)
                        !p(isd,node)=p(isd,node) + mu*ph(isd,inl)*(dr(1,1)+dr(2,2))*2.0/3.0
                    enddo
                    q_p(1,node)=q_p(1,node)+ph(1,inl)*mu*(sh(1,inl)*q_d(1,inl)*(2.0-2.0/3.0)) + &
                                            ph(2,inl)*mu*sh(2,inl)*q_d(1,inl)
                    q_p(2,node)=q_p(2,node)+ph(1,inl)*mu*sh(1,inl)*q_d(2,inl) + &
                                            ph(2,inl)*(mu*sh(2,inl)*q_d(2,inl)*(2.0-2.0/3.0))
                elseif (nsd==3) then
                    do isd=1,nsd
                        p(isd,node)=p(isd,node) + ph(isd,inl)*pp -   &
                                                  ph(1,inl)*tau(1,isd) -  &
                                                  ph(2,inl)*tau(2,isd) -  &
                                                  ph(3,inl)*tau(3,isd)
                        !p(isd,node)=p(isd,node)+mu*ph(isd,inl)*(dr(1,1)+dr(2,2)+dr(3,3))*2.0/3.0
                    enddo
                    q_p(1,node)=q_p(1,node)+ph(1,inl)*mu*(sh(1,inl)*(2.0 - 2.0/3.0))+&
                                            ph(2,inl)*mu*sh(2,inl)+ph(3,inl)*mu*sh(3,inl)
                    q_p(2,node)=q_p(2,node)+ph(1,inl)*mu*sh(1,inl)+ph(2,inl)*mu*(sh(2,inl)*(2.0 -2.0/3.0))+&
                                            ph(3,inl)*mu*sh(3,inl)
                    q_p(3,node)=q_p(3,node)+ph(1,inl)*mu*sh(1,inl)+ph(2,inl)*mu*sh(2,inl)+&
                                            ph(3,inl)*mu*(sh(3,inl)*(2.0-2.0/3.0))
                endif
                ! Stablization with Tau_moment
                if (nsd==2) then
                    p(pdf,node) = p(pdf,node) - ph(xsd,inl)*prs_cc(udf)  &
                                              - ph(ysd,inl)*prs_cc(vdf)
                    q_p(pdf,node) = q_p(pdf,node)+tauc*(sh(1,inl)*sh(1,inl)+sh(2,inl)*sh(2,inl))*eft0
                elseif (nsd==3) then
                    p(pdf,node) = p(pdf,node) - ph(xsd,inl)*prs_cc(udf)  &
                                              - ph(ysd,inl)*prs_cc(vdf)  &
                                              - ph(zsd,inl)*prs_cc(wdf)
                    q_p(pdf,node) = q_p(pdf,node)+tauc*(sh(1,inl)**2+sh(2,inl)**2+sh(3,inl)**2)*eft0
                endif             ! Stablization with Tau_cont    
                p(1:nsd,node) = p(1:nsd,node) - prs_t(1:nsd)*temp - ph(1:nsd,inl)*prs_c
                q_p(1:nsd,node) = q_p(1:nsd,node)+taum*temp*q_res_a(1:nsd,inl)+ & 
                                         ph(1:nsd,inl)*taul*ro*sh(1:nsd,inl)
            enddo
        enddo ! end of qudrature pts loop
    enddo ! end of element loop
continue
!=====================================
! Apply boundary condition du/dx=0 on outedge
!call out2d4n(rngface,dloc,ien,xloc,p)
!======================================
end subroutine block