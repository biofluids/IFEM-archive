!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  L. Zhang
!  Northwestern University
!  This subroutine solves for the residual for all degrees of freedom
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine block(xloc, dloc, doloc, p, q, hk, ien, f_fluids)
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

  real* 8 drt(ndf),drs(ndf)
!  real* 8 drx(ndf),dry(ndf),drz(ndf)
  real* 8 dr(nsd,ndf)
!  real* 8 u,v,w,pp,ug
  real* 8 u(1:nsd),pp,ug
!  real* 8 txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
  real* 8 hg,taum,tauc,vel,ree
  real* 8 res_c,res_a(nsd),res_t(nsd)
  real* 8 prs_c,prs_t(nsd)
  real* 8 mu,nu,ro,g(nsd)
!  real* 8 tempu,tempv,tempw,temp
  real* 8 dtinv,oma,ama
  integer inl, ie, isd, idf, iq, node,jsd
  real* 8 tempc(ndf),temp
  real* 8 t(nsd,nsd)


  real* 8 f_fluids(nsd,nn)
  real* 8 fnode(nsd,nen),fq(nsd)

  dtinv = 1.0/dt
  if(steady) dtinv = 0.0
  oma   = 1.0 - alpha
  ama   = 1.0 - oma

  do ie=1,ne		! loop over elements
     do inl=1,nen		
        do isd=1,nsd 
           x(isd,inl) = xloc(isd,ien(inl,ie))
		   fnode(isd,inl) = f_fluids(isd,ien(inl,ie))
	    enddo
	    do idf=1,ndf
		   d(idf,inl) =  dloc(idf,ien(inl,ie))
		   d_old(idf,inl) = doloc(idf,ien(inl,ie))
	    enddo
	 enddo

	 hg = hk(ie)

	 do iq=1,nquad  ! loop over the quadrature points in each element 
!...  calculate the shape function and the weight at quad point
	    if (nen.eq.4) then !calculate shape function at quad point
		   include "sh3d4n.h"
	    else if (nen.eq.8) then
		   include "sh3d8n.h"
	    endif

	    eft0 = abs(det) * wq(iq) ! calculate the weight at each quad pt
!...  initialize d, dd/dx, dd/dy, dd/dz, dd/dt
        do idf = 1,ndf
           drs(idf) = 0.0
!		   drx(idf) = 0.0
!		   dry(idf) = 0.0
!		   drz(idf) = 0.0
		   drt(idf) = 0.0
		   dr(:,:)=0.0
	    enddo
		fq(:)=0.0
!... calculate vi, dvi/dxj
        do inl=1,nen
		   tempc(1:nsd)=ama*d(1:nsd,inl)+oma*d_old(1:nsd,inl)
		   drs(isd)=drs(isd)+sh(0,inl)*tempc(isd)
		   do idf=1,ndf

!			 drx(isd)=drx(isd)+sh(1,inl)*tempc(isd)           
!		     dry(isd)=dry(isd)+sh(2,inl)*tempc(isd)           
!		     drz(isd)=drz(isd)+sh(3,inl)*tempc(isd)
			 dr(1:nsd,idf)=dr(1:nsd,idf)+sh(1:nsd,inl)*tempc(idf)
		   enddo
		   fq(:) = fq(:) + sh(0,inl)*fnode(:,inl)        
	    enddo

!... calculate dvi/dt, p, dp/dxi
        do inl=1,nen
		   drt(1:nsd)=drt(1:nsd)+sh(0,inl)*(d(1:nsd,inl)-d_old(1:nsd,inl))*dtinv
		   drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)
!		   drx(pdf)=drx(pdf)+sh(1,inl)*d(pdf,inl)      
!		   dry(pdf)=dry(pdf)+sh(2,inl)*d(pdf,inl)      
!		   drz(pdf)=drz(pdf)+sh(3,inl)*d(pdf,inl)      
		   dr(1:nsd,pdf)=dr(1:nsd,pdf)+sh(1:nsd,inl)*d(pdf,inl)
	    enddo

!... define u=v1, v=v2, w=v3, pp=p
!	    u = drs(udf)
!	    v = drs(vdf)
!	    w = drs(wdf)
		u=drs(1:nsd)
	    pp= drs(pdf)

	    if(stokes) then ! if stokes flow
!		   u = 0.0
!		   v = 0.0
!		   w = 0.0
		   u(1:nsd)=0.0
	    endif

!....  calculate liquid constant and gravity
	    mu = vis_liq  ! liquid viscosity
	    ro = den_liq  ! liquid density
	    g(1:nsd) = gravity(1:nsd) ! gravatitional force

	! believe nu is calculated only for turbulent model
!	    nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2 &
!	                                            +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2 &
!	                                            +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)
	    nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*dr(1,1)**2+(dr(2,1)+dr(1,2))**2	&
	                                            +2*dr(2,2)**2+(dr(3,1)+dr(1,3))**2	&
	                                            +2*dr(3,3)**2+(dr(3,2)+dr(2,3))**2)

!		nu=0
!		do isd=1,nsd
!			do jsd=1,nsd
!			   nu=nu+(dr(isd,jsd)+dr(jsd,isd))**2
!			enddo
!		enddo
!		nu=delta(4)*turb_kappa**2*hg**2*sqrt(nu)

	    mu = mu + nu*ro                   

!....  calculate each term in the residual equation
	    res_c = 0.0
	    do inl=1,nen
!		   res_c = res_c+sh(xsd,inl)*d(udf,inl) &
!	                    +sh(ysd,inl)*d(vdf,inl) &
!	                    +sh(zsd,inl)*d(wdf,inl)
			do isd=1,nsd
			  res_c=res_c+sh(isd,inl)*d(isd,inl)
			enddo
	    enddo

	    do isd = 1, nsd
!		   res_a(isd)=ro*(drt(isd)+u*drx(isd)+v*dry(isd)+w*drz(isd)-g(isd))-fq(isd)
		   res_a(isd)=ro*(drt(isd)+u(1)*dr(1,isd)+u(2)*dr(2,isd)+u(3)*dr(3,isd)-g(isd))-fq(isd)
!			temp=0
!			do jsd=1,nsd
!				temp = temp + u(jsd)*dr(jsd,isd)
!			enddo
!			res_a(isd)=ro*(drt(isd)+temp-g(isd))-fq(isd)
	    enddo

!	    res_t(xsd) = drx(pdf) + res_a(xsd) 
!	    res_t(ysd) = dry(pdf) + res_a(ysd) 
!	    res_t(zsd) = drz(pdf) + res_a(zsd)

		res_t(1:nsd)=dr(1:nsd,pdf)+res_a(1:nsd)

!.....  Stablization parameters, taum and tauc
!        vel  = sqrt(u*u+v*v+w*w)  !magnitude of the velocity
		vel=0
		do isd=1,nsd
			vel=vel+u(isd)*u(isd)
		enddo
		vel=sqrt(vel)

	    ree  = vel*hg/mu/12.0  !???
	    if(steady.or.(.not.taudt)) then !stablization, taum
		   taum = 1.0/sqrt((2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
	    else
		   taum = 1.0/sqrt((2.0/dt)**2+(2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
	    endif
	    taum = delta(1)* taum
	    tauc = delta(2)*hg*vel
	    if(ree.lt.1.0) tauc = tauc*ree

	    taum = taum/ro
	    tauc = tauc*ro 

!.....   Density optimization
        do inl=1,nen
		   ph(0,inl) = sh(0,inl)*eft0
!		   ph(1,inl) = sh(1,inl)*eft0
!		   ph(2,inl) = sh(2,inl)*eft0
!		   ph(3,inl) = sh(3,inl)*eft0
		   ph(1:nsd,inl)=sh(1:nsd,inl)*eft0
	    enddo
	      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....   Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c....   Calculate stress tau(ij) and pressure
!	    txx = mu*(drx(udf)+drx(udf))
!	    tyx = mu*(dry(udf)+drx(vdf))
!	    tzx = mu*(drz(udf)+drx(wdf))
!	    txy = mu*(drx(vdf)+dry(udf))
!	    tyy = mu*(dry(vdf)+dry(vdf))
!	    tzy = mu*(drz(vdf)+dry(wdf))
!	    txz = mu*(drx(wdf)+drz(udf))
!	    tyz = mu*(dry(wdf)+drz(vdf))
!	    tzz = mu*(drz(wdf)+drz(wdf))
!	    prs_t(udf) = res_t(udf) * taum
!	    prs_t(vdf) = res_t(vdf) * taum
!	    prs_t(wdf) = res_t(wdf) * taum
!	    prs_c      = res_c*tauc
		do isd=1,nsd
			do jsd=1,nsd
				t(isd,jsd)=mu*(dr(isd,jsd)+dr(jsd,isd))
			enddo
		enddo
		prs_t(1:nsd)=res_t(1:nsd)*taum
		prs_c		=res_c*tauc


	      
!.... calculate the residual at each degree of freedom
	    do inl=1,nen ! loop over number of nodes in an element
		 
		   node=ien(inl,ie)

!		   temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))
			do isd=1,nsd
				temp=temp+ro*(u(isd)*ph(isd,inl))
			enddo
		! Continuty Equation
		   p(pdf,node) = p(pdf,node)-ph(0,inl)*res_c

		! Momentum Equation (Euler Residual)
!		   p(udf,node) = p(udf,node)-ph(0,inl)*res_a(udf)
!		   p(vdf,node) = p(vdf,node)-ph(0,inl)*res_a(vdf)
!		   p(wdf,node) = p(wdf,node)-ph(0,inl)*res_a(wdf)
		   p(1:nsd,node)=p(1:nsd,node)-ph(0,inl)*res_a(1:nsd)

		! Momentum Equation (C: Diffusion terms)
!		   p(udf,node) = p(udf,node) + ph(xsd,inl)* pp -    &
!	                                   ph(xsd,inl)*txx -    &
!	                                   ph(ysd,inl)*tyx -    &
!	                                   ph(zsd,inl)*tzx 
!		   p(vdf,node) = p(vdf,node) + ph(ysd,inl)* pp -    &
!	                                   ph(xsd,inl)*txy -    &
!	                                   ph(ysd,inl)*tyy -    &
!	                                   ph(zsd,inl)*tzy 
!		   p(wdf,node) = p(wdf,node) + ph(zsd,inl)* pp -    &
!	                                   ph(xsd,inl)*txz -    &
!	                                   ph(ysd,inl)*tyz -    &
!	                                   ph(zsd,inl)*tzz 

		  do isd=1,nsd
		    do jsd=1,nsd
			    p(isd,node) = p(isd,node) + ph(isd,inl)*pp -    &
											ph(jsd,inl)*t(jsd,isd)
			enddo
		  enddo


		! Stablization with Tau_moment
!		   p(pdf,node) = p(pdf,node) - ph(xsd,inl)*prs_t(udf)  &
!	                                 - ph(ysd,inl)*prs_t(vdf)  &
!	                                 - ph(zsd,inl)*prs_t(wdf)
!		   p(udf,node) = p(udf,node) -prs_t(udf)*temp
!		   p(vdf,node) = p(vdf,node) -prs_t(vdf)*temp
!		   p(wdf,node) = p(wdf,node) -prs_t(wdf)*temp

		! Stablization with Tau_cont    
!		   p(udf,node) = p(udf,node) - ph(xsd,inl)*prs_c
!		   p(vdf,node) = p(vdf,node) - ph(ysd,inl)*prs_c
!		   p(wdf,node) = p(wdf,node) - ph(zsd,inl)*prs_c

		  do isd=1,nsd
		    p(pdf,node)=p(pdf,node)-ph(isd,inl)*prs_t(isd)
			p(isd,node)=p(isd,node)-prs_t(isd)*temp
			p(isd,node)=p(isd,node)-ph(isd,inl)*prs_c
		  enddo

	    enddo

	      ! Diagonal Preconditioner

	    effd =   mu*eft0*alpha
	    effm = taum*eft0
	    effc = tauc*eft0

	    do inl=1,nen

		   node=ien(inl,ie)
!		   ug = ro*(u*sh(1,inl)+v*sh(2,inl)+w*sh(3,inl))
			ug=0
		   do isd=1,nsd
		      ug=ug+ro*(u(isd)*sh(isd,inl))
		   enddo

		   temp = alpha*ug + sh(0,inl)*dtinv*ro
		 
!		   q(udf,node) = q(udf,node) + (alpha*ro*sh(0,inl)*drx(1)+temp)*(sh(0,inl)*eft0+ug*effm)
!		   q(vdf,node) = q(vdf,node) + (alpha*ro*sh(0,inl)*dry(2)+temp)*(sh(0,inl)*eft0+ug*effm)
!          q(wdf,node) = q(wdf,node) + (alpha*ro*sh(0,inl)*drz(3)+temp)*(sh(0,inl)*eft0+ug*effm)
		   do isd=1,nsd
		     q(isd,node)=q(isd,node)+(alpha*ro*sh(0,inl)*dr(isd,isd)+temp)*(sh(0,inl)*eft0+ug*effm)
		   enddo

!		   temp = sh(1,inl)**2+sh(2,inl)**2+sh(3,inl)**2

		   temp=0
		   do isd=1,nsd
		     temp=temp+sh(isd,inl)**2
		   enddo

!		   q(udf,node) = q(udf,node)+(sh(1,inl)**2+temp)*effd
!		   q(vdf,node) = q(vdf,node)+(sh(2,inl)**2+temp)*effd
!		   q(wdf,node) = q(wdf,node)+(sh(3,inl)**2+temp)*effd
!		   q(pdf,node) = q(pdf,node)+temp*effm
!		   q(udf,node) = q(udf,node)+sh(1,inl)**2*effc
!		   q(vdf,node) = q(vdf,node)+sh(2,inl)**2*effc
!		   q(wdf,node) = q(wdf,node)+sh(3,inl)**2*effc

		   q(1:nsd,node) = q(1:nsd,node)+(sh(1:nsd,inl)**2+temp)*effd   &
							+ sh(1:nsd,inl)**2 * effc
		   q(pdf,node) = q(pdf,node)+temp*effm

	    enddo

	 enddo ! end of qudrature pts loop
  enddo ! end of element loop

  return
end subroutine block

