!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  L. zhang
!  Northwestern Univeristy
!  blockgmres.f
!  called by gmres.f
!  this subroutine calculate the increment of the residual
!  L. Zhang, 06/25/2004
!  Tulane University
!  Revised the subroutine to array
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine blockgmres(xloc,dloc,doloc,qloc,p,hk,ien,fext)
  use global_constants
  use run_variables, only: dt
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real* 8 xloc(nsd,nn)
  real* 8 dloc(ndf,nn),doloc(ndf,nn)
  real* 8 p(ndf,nn),qloc(ndf,nn),hk(ne)

  real* 8 x(nsd,nen)
  real* 8 d(ndf,nen),d_old(ndf,nen),q(ndf,nen)

  real* 8 eft0,det
  real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
  real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real* 8 drs(ndf),qrt(ndf),qrs(ndf)
  real* 8 dr(nsd,ndf)
  real* 8 qr(nsd,ndf)
  real* 8 u,v,w,pp
  real* 8 tau(nsd,nsd)
  real* 8 hg,taum,tauc,vel,ree, taul
  real* 8 res_c,res_a(nsd),res_t(nsd)
  real* 8 prs_c,prs_t(nsd),prs_cc(nsd)
  real* 8 mu,nu,ro
  real* 8 tempc(ndf),temp
  real* 8 dtinv,oma,ama
  integer inl, ie, isd, iq, node, jsd
  real*8 fext(nsd,nn)
  real* 8 fnode(nsd,nen),fq(nsd)

!.....calculate 1/dt
  dtinv = 1.0/dt/alpha
  if(steady) dtinv = 0.0

!.....coefficient for the residuals
  oma   = 1.0 -alpha
  ama   = 1.0 - oma

  do ie=1,ne		! loop over the elements
!...    localize x and degrees of freedom in every node of the element
     do inl=1,nen
         x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
		 fnode(1:nsd,inl) = fext(1:nsd,ien(inl,ie))
		 q(1:ndf,inl) =  qloc(1:ndf,ien(inl,ie))
		 d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
		 d_old(1:ndf,inl) = doloc(1:ndf,ien(inl,ie))
	 enddo

	 hg = hk(ie)

	 do iq=1,nquad ! loop over the quad. pts in each element
!....    calculate shape function and weight at the quad pt
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

	      eft0 = abs(det) * wq(iq) * alpha
!....    initialize vi,dxi/dxj 

		drs(:) = 0.0
		dr(1:nsd,1:ndf)=0.0
		fq(:)=0.0

!....    calculate vi,dvi/dxj
        do inl=1,nen
		   tempc(1:nsd) = ama*d(1:nsd,inl)+oma*d_old(1:nsd,inl)
		   drs(1:nsd) = drs(1:nsd)+sh(0,inl)*tempc(1:nsd)
		   do isd=1,nsd
			 dr(isd,1:nsd) = dr(isd,1:nsd)+sh(isd,inl)*tempc(1:nsd)
		   enddo
		   fq(:) = fq(:) + sh(0,inl)*fnode(:,inl)        
	    enddo

!...     initialize delta_d, delta_p, d(delta_v)i/dxj, d(delta_v)i/dt

		qrs(:) = 0.0
		qrt(:) = 0.0
		qr(1:nsd,1:ndf)=0.0

!.....   calculate delta_v, d(delta_v)i/dxj
        do inl=1,nen
		   qrs(1:nsd) = qrs(1:nsd)+sh(0,inl)*q(1:nsd,inl)
		   do isd=1,nsd
			 qr(isd,1:nsd) = qr(isd,1:nsd)+sh(isd,inl)*q(1:nsd,inl)
		   enddo
	    enddo

!.....   calculate d(delta_v)/dt, delta_p, and d(delta_p)/dxj
        do inl=1,nen
		   qrt(1:nsd)=qrt(1:nsd)+sh(0,inl)*q(1:nsd,inl)*dtinv
		   qrs(pdf)=qrs(pdf)+sh(0,inl)*q(pdf,inl)    		   
		   qr(1:nsd,pdf)=qr(1:nsd,pdf)+sh(1:nsd,inl)*q(pdf,inl)       
	    enddo

!....    reset u=x_velocity, v=y_velocity, w=z_velocity, pp=pressure
		if (nsd==2) then
		    u = drs(udf)
			v = drs(vdf)
			w = 0.0
		elseif (nsd==3) then
			u = drs(udf)
		    v = drs(vdf)
		    w = drs(wdf)
		endif
	    pp= qrs(pdf)/alpha

	    if(stokes) then
		   u = 0.0
		   v = 0.0
		   w = 0.0
	    endif

!....    set liquid properties, density and viscosity
	    mu = vis_liq
	    ro = den_liq
              ! below is only used if turbulence is applied
		if (nsd==2) then
			nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*dr(1,1)**2+(dr(2,1)+dr(1,2))**2 &
	                                            +2*dr(2,2)**2)
		elseif (nsd==3) then
			nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*dr(1,1)**2+(dr(2,1)+dr(1,2))**2 &
	                                            +2*dr(2,2)**2+(dr(3,1)+dr(1,3))**2 &
	                                            +2*dr(3,3)**2+(dr(3,2)+dr(2,3))**2)
		endif
	    mu = mu + nu*ro                   

!.....   calculate the delta of residuals
	    res_c = 0.0
		if (nsd==2) then
		  do inl=1,nen
		   	 res_c = res_c+(sh(xsd,inl)*q(udf,inl) &
	                    +sh(ysd,inl)*q(vdf,inl))/alpha
		  enddo
		elseif (nsd==3) then
		  do inl=1,nen
		     res_c = res_c+(sh(xsd,inl)*q(udf,inl) &
	                    +sh(ysd,inl)*q(vdf,inl) &
	                    +sh(zsd,inl)*q(wdf,inl))/alpha
		  enddo
		endif

	    do isd = 1, nsd
			if (nsd==2) then
			   res_a(isd)=ro*(qrt(isd)+u*qr(1,isd)+v*qr(2,isd))
			elseif (nsd==3) then
			   res_a(isd)=ro*(qrt(isd)+u*qr(1,isd)+v*qr(2,isd)+w*qr(3,isd))
			endif
	    enddo

		res_t(1:nsd) = qr(1:nsd,pdf)/alpha + res_a(1:nsd)


		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0
	      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.....   Galerkin Terms (Look at notes)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.....   calculate stress and pressure terms
		do isd = 1,nsd
			do jsd = 1,nsd
				tau(isd,jsd) = mu*(qr(isd,jsd) + qr(jsd,isd))
			enddo
		enddo

		prs_t(1:nsd) = res_t(1:nsd)*taum

		!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		! Mickael 02/01/2005
		! TAUl (l for  lsic), See Tezduyar and Sathe, 2003

        prs_c = ro*res_c*taul
		!
		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	    prs_cc(1:nsd) = res_t(1:nsd)*tauc	      

!.... assemble the delta-residuals at nodes
	    do inl=1,nen
		   node = ien(inl,ie)
		   if (nsd==2) then
			   temp = ro*(u*ph(xsd,inl)+v*ph(ysd,inl))
		   elseif (nsd==3) then
			   temp = ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))
		   endif

!.....      Continuty Equation
		   p(pdf,node) = p(pdf,node)+ph(0,inl)*res_c

!.....      Momentum Equation (Euler Residual)
		   p(1:nsd,node) = p(1:nsd,node)+ph(0,inl)*res_a(1:nsd)

!.....      Momentum Equation (C: Diffusion terms)
		    if (nsd==2) then
			  do isd=1,nsd
				p(isd,node)=p(isd,node) - ph(isd,inl)*pp +   &
										  ph(1,inl)*tau(1,isd) +  &
										  ph(2,inl)*tau(2,isd)
			  enddo
			elseif (nsd==3) then
			  do isd=1,nsd
				p(isd,node)=p(isd,node) - ph(isd,inl)*pp +   &
										  ph(1,inl)*tau(1,isd) +  &
										  ph(2,inl)*tau(2,isd) +  &
										  ph(3,inl)*tau(3,isd)
			  enddo
			endif

!.....      Stablization with Tau_moment
		   if (nsd==2) then
		   	 p(pdf,node) = p(pdf,node) + ph(xsd,inl)*prs_cc(udf)  &
	                                   + ph(ysd,inl)*prs_cc(vdf)
	       elseif (nsd==3) then
		     p(pdf,node) = p(pdf,node) + ph(xsd,inl)*prs_cc(udf)  &
	                                   + ph(ysd,inl)*prs_cc(vdf)  &
	                                   + ph(zsd,inl)*prs_cc(wdf)
		   endif		! Stablization with Tau_cont 
		   p(1:nsd,node) = p(1:nsd,node) + prs_t(1:nsd)*temp + ph(1:nsd,inl)*prs_c

	    enddo


	    if(stokes) goto 500
!.....   Non-linear Terms (from momentum eq.) 

		   if (nsd==2) then
		       do isd = 1, nsd
		  	      res_a(isd)=ro*(qrs(1)*dr(1,isd)+qrs(2)*dr(2,isd))
			   enddo
		   elseif (nsd==3) then
		       do isd = 1, nsd
		  	      res_a(isd)=ro*(qrs(1)*dr(1,isd)+qrs(2)*dr(2,isd)+qrs(3)*dr(3,isd))
			   enddo
		   endif

	       do inl=1,nen
		      node = ien(inl,ie)
			  if (nsd==2) then
				temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl))
			  elseif (nsd==3) then
			    temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))
			  endif

!.......    Additional nonlinear term form Galerkin
			  p(1:nsd,node)=p(1:nsd,node)+ph(0,inl)*res_a(1:nsd)

!.......    Additional nonlinear term form SUPG
			  if (nsd==2) then
			      p(pdf,node) = p(pdf,node) + taum*(ph(xsd,inl)*res_a(udf)+ph(ysd,inl)*res_a(vdf))
			  elseif (nsd==3) then
			      p(pdf,node) = p(pdf,node) + taum*(ph(xsd,inl)*res_a(udf)+ph(ysd,inl)*res_a(vdf)+ph(zsd,inl)*res_a(wdf))
			  endif
		      p(1:nsd,node) = p(1:nsd,node) + res_a(1:nsd)*temp*taum
	       enddo

 500	continue 
     enddo  ! end of quad loop
  enddo  ! end of element loop

  return
end subroutine blockgmres
