!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  L. Zhang
!  Northwestern University
!  This subroutine solves for the residual for all degrees of freedom
!____________________________________________________________________
!  L. Zhang, 06/24/2004
!  Tulane University
!  Revised the subroutine to array
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine blockgmresnew(xloc, dloc, doloc, p, hk, ien, f_fluids)
  use global_constants
  use run_variables
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real* 8 xloc(nsd,nn)
  real* 8 dloc(ndf,nn),doloc(ndf,nn)
  real* 8 p(ndf,nn),q(ndf,nn),hk(ne)

  real* 8 x(nsd,nen)
  real* 8 d(ndf,nen),d_old(ndf,nen)

  real* 8 eft0,det,effd,effm,effc
  real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
  real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real* 8 f_stress(nsd,nsd,nn)

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
  integer ieface,irng, rngface(neface,ne) !,inface

  real* 8 f_fluids(nsd,nn)
  real* 8 fnode(nsd,nen),fq(nsd)
  real* 8 com_kappa 

  com_kappa = 1.0e7

  dtinv = 1.0/dt
  if(steady) dtinv = 0.0
  oma   = 1.0 - alpha
  ama   = 1.0 - oma
 !=================================================
!f_fluids(:,:)=f_fluids(:,:)/(0.0625/6.0)
p(1:nsd,1:nn)=p(1:nsd,1:nn)+f_fluids(1:nsd,1:nn)
!=================================================
!===================================================
! f_fluids is actually the FSI force at fluid nodes,
! then we will just subscrib it from p(!:nsd) which is the 
! residuals for the momentum equations. 
! 1 let fnod=0 then fndoe and fq 's effects will be canceled out
! 2 do the subscribition after the elements loop
! Xingshi 09/15/2008
!===================================================
  do ie=1,ne		! loop over elements
     do inl=1,nen	
	     x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
!============================================================================
!		 fnode(1:nsd,inl) = f_fluids(1:nsd,ien(inl,ie))	
                 fnode(:,inl)=0.0
!============================================================================
		 d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
		 d_old(1:ndf,inl) = doloc(1:ndf,ien(inl,ie))
	 enddo

	 hg = hk(ie)


	 

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
!		stop
!	write(*,*) 'shape functions=',sh(0,1),sh(0,2),sh(0,3),sh(0,4)
	    eft0 = abs(det) * wq(iq) ! calculate the weight at each quad pt
!...  initialize d, dd/dx, dd/dy, dd/dz, dd/dt
		drs(:) = 0.0
		drt(:) = 0.0
		dr(1:nsd,1:ndf)=0.0
		fq(:)=0.0

!... calculate vi, dvi/dxj
        do inl=1,nen
		   tempc(1:nsd) = ama*d(1:nsd,inl)+oma*d_old(1:nsd,inl)
		   drs(1:nsd) = drs(1:nsd)+sh(0,inl)*tempc(1:nsd)
		   do isd=1,nsd
			 dr(isd,1:nsd) = dr(isd,1:nsd)+sh(isd,inl)*tempc(1:nsd)
		   enddo
		   fq(:) = fq(:) + sh(0,inl)*fnode(:,inl)        
	    enddo


!... calculate dvi/dt, p, dp/dxi
        do inl=1,nen
		   drt(1:ndf)=drt(1:ndf)+sh(0,inl)*(d(1:ndf,inl)-d_old(1:ndf,inl))*dtinv
		   drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)    		   
		   dr(1:nsd,pdf)=dr(1:nsd,pdf)+sh(1:nsd,inl)*d(pdf,inl)       
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
	    mu = vis_liq  ! liquid viscosity
	    ro = den_liq  ! liquid density
		g  = gravity  ! gravatitional force

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
	    res_c = 0.0
		if (nsd==2) then
		  do inl=1,nen
		   	 res_c = res_c+sh(xsd,inl)*d(udf,inl) &
	                    +sh(ysd,inl)*d(vdf,inl)
		  enddo
		elseif (nsd==3) then
		  do inl=1,nen
		     res_c = res_c+sh(xsd,inl)*d(udf,inl) &
	                    +sh(ysd,inl)*d(vdf,inl) &
	                    +sh(zsd,inl)*d(wdf,inl)
		  enddo
		endif

		res_c=res_c+drt(pdf)/com_kappa !add p_dot term for continuity equation

	    do isd = 1, nsd
			if (nsd==2) then
			   res_a(isd)=ro*(drt(isd)+u*dr(1,isd)+v*dr(2,isd)-g(isd))-fq(isd)
			  

			elseif (nsd==3) then
			   res_a(isd)=ro*(drt(isd)+u*dr(1,isd)+v*dr(2,isd)+w*dr(3,isd)-g(isd))-fq(isd)
			endif
	    enddo

		res_t(1:nsd) = dr(1:nsd,pdf) + res_a(1:nsd)

		!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! Mickael 02/01/2005
        ! TAUm, TAUc and TAUl (l for  lsic), See Tezduyar and Sathe, 2003

                                
        vel  = sqrt(u*u+v*v+w*w)  !magnitude of the velocity
                        
        ree  = vel*hg*ro/(2.0*mu)  !Reynolds number
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
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!.....   Density optimization

		ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0
	      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....   Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c....   Calculate stress tau(ij) and pressure
		do isd = 1,nsd
			do jsd = 1,nsd
				tau(isd,jsd) = mu*(dr(isd,jsd) + dr(jsd,isd))
				do inl=1,nen	
					! Compute fluid stress
					f_stress(isd,jsd,ien(inl,ie))= f_stress(isd,jsd,ien(inl,ie)) + tau(isd,jsd)
				enddo
			enddo
		enddo

		prs_t(1:nsd) = res_t(1:nsd)*taum


		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! Mickael 02/01/2005
        ! TAUl (l for  lsic), See Tezduyar and Sathe, 2003
  
        prs_c = ro*res_c*taul
        !
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
		! Continuty Equation
		   p(pdf,node) = p(pdf,node)-ph(0,inl)*res_c
		! Momentum Equation (Euler Residual)
!===========================================================
! why originally it is minus?
		   p(1:nsd,node) = p(1:nsd,node)-ph(0,inl)*res_a(1:nsd)
!=============================================================
		! Momentum Equation (C: Diffusion terms)
		    if (nsd==2) then
			  do isd=1,nsd
				p(isd,node)=p(isd,node) + ph(isd,inl)*pp -   &
										  ph(1,inl)*tau(1,isd) -  &
										  ph(2,inl)*tau(2,isd)
				p(isd,node)=p(isd,node)+mu*ph(isd,inl)*(dr(1,1)+dr(2,2))*2/3!add mu*2/3*du_k/dx_k
!				diag(inl)=p(1,node)
			  enddo

			elseif (nsd==3) then
			  do isd=1,nsd
				p(isd,node)=p(isd,node) + ph(isd,inl)*pp -   &
										  ph(1,inl)*tau(1,isd) -  &
										  ph(2,inl)*tau(2,isd) -  &
										  ph(3,inl)*tau(3,isd)
			  enddo
			endif

		! Stablization with Tau_moment
		   if (nsd==2) then
		   	 p(pdf,node) = p(pdf,node) - ph(xsd,inl)*prs_cc(udf)  &
	                                   - ph(ysd,inl)*prs_cc(vdf)
!                      diag(inl)=p(pdf,node)
	       elseif (nsd==3) then
		     p(pdf,node) = p(pdf,node) - ph(xsd,inl)*prs_cc(udf)  &
	                                   - ph(ysd,inl)*prs_cc(vdf)  &
	                                   - ph(zsd,inl)*prs_cc(wdf)
		   endif		! Stablization with Tau_cont    
		   p(1:nsd,node) = p(1:nsd,node) - prs_t(1:nsd)*temp - ph(1:nsd,inl)*prs_c
!                   diag(inl)=p(1,node)

	 enddo


	 enddo ! end of qudrature pts loop

!write(*,*)'diag=',diag(1:4)

  enddo ! end of element loop
!write(*,*)q_p(:,:)

 continue  
continue
!write(*,*)q_p(:,:)
!write(*,*)d(:,:)

end subroutine blockgmresnew

