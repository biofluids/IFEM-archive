!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  L. zhang
!  Northwestern Univeristy
!  blockgmres.f
!  called by gmres.f
!  this subroutine calculate the increment of the residual
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

  real* 8 eft0,det !,effd,effm,effc
  real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
  real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real* 8 drs(ndf),qrt(ndf),qrs(ndf)
  real* 8 drx(ndf),dry(ndf),drz(ndf)
  real* 8 qrx(ndf),qry(ndf),qrz(ndf)
  real* 8 u,v,w,pp !,ug
  real* 8 txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
  real* 8 hg,taum,tauc,vel,ree
  real* 8 res_c,res_q,res_a(nsd),res_t(nsd)
  real* 8 prs_c,prs_t(nsd)
  real* 8 mu,nu,ro
  real* 8 tempu,tempv,tempw,temp
  real* 8 dtinv,oma,ama
  integer inl, ie, isd, idf, iq, node
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
        do isd=1,nsd
           x(isd,inl) = xloc(isd,ien(inl,ie))
		   fnode(isd,inl) = fext(isd,ien(inl,ie))
	    enddo
	    do idf=1,ndf   
			 q(idf,inl) =  qloc(idf,ien(inl,ie))
			 d(idf,inl) =  dloc(idf,ien(inl,ie))
			 d_old(idf,inl) = doloc(idf,ien(inl,ie))
	    enddo
	 enddo

	 hg = hk(ie)

	 do iq=1,nquad ! loop over the quad. pts in each element
!....    calculate shape function and weight at the quad pt
	    if (nen.eq.4) then
	       include "sh3d4n.h"
	    else if (nen.eq.8) then
		   include "sh3d8n.h"
	    end if

	      eft0 = abs(det) * wq(iq) * alpha
!....    initialize vi,dxi/dxj 
	    do idf = 1,nsd
		   drs(idf) = 0.0
		   drx(idf) = 0.0
		   dry(idf) = 0.0
		   drz(idf) = 0.0
	    enddo
		fq(:)=0.0
!....    calculate vi,dvi/dxj
	    do inl=1,nen
			 tempu= ama*d(udf,inl) + oma*d_old(udf,inl)
			 tempv= ama*d(vdf,inl) + oma*d_old(vdf,inl)
			 tempw= ama*d(wdf,inl) + oma*d_old(wdf,inl)
			 drs(udf)=drs(udf)+sh(0,inl)*tempu            
			 drx(udf)=drx(udf)+sh(1,inl)*tempu           
			 dry(udf)=dry(udf)+sh(2,inl)*tempu          
			 drz(udf)=drz(udf)+sh(3,inl)*tempu         
			 drs(vdf)=drs(vdf)+sh(0,inl)*tempv            
			 drx(vdf)=drx(vdf)+sh(1,inl)*tempv           
			 dry(vdf)=dry(vdf)+sh(2,inl)*tempv          
			 drz(vdf)=drz(vdf)+sh(3,inl)*tempv         
			 drs(wdf)=drs(wdf)+sh(0,inl)*tempw            
			 drx(wdf)=drx(wdf)+sh(1,inl)*tempw           
			 dry(wdf)=dry(wdf)+sh(2,inl)*tempw          
			 drz(wdf)=drz(wdf)+sh(3,inl)*tempw
			 fq(:)=fq(:)+sh(0,inl)*fnode(:,inl)      
	    end do
!...     initialize delta_d, delta_p, d(delta_v)i/dxj, d(delta_v)i/dt
	    do idf=1,ndf
		   qrs(idf) = 0.0
		   qrx(idf) = 0.0
		   qry(idf) = 0.0
		   qrz(idf) = 0.0
		   qrt(idf) = 0.0
	    enddo
!.....   calculate delta_v, d(delta_v)i/dxj
	    do inl=1,nen
			 qrs(udf)=qrs(udf)+sh(0,inl)*q(udf,inl)      
			 qrx(udf)=qrx(udf)+sh(1,inl)*q(udf,inl)      
			 qry(udf)=qry(udf)+sh(2,inl)*q(udf,inl)      
			 qrz(udf)=qrz(udf)+sh(3,inl)*q(udf,inl)      
			 qrs(vdf)=qrs(vdf)+sh(0,inl)*q(vdf,inl)      
			 qrx(vdf)=qrx(vdf)+sh(1,inl)*q(vdf,inl)      
			 qry(vdf)=qry(vdf)+sh(2,inl)*q(vdf,inl)      
			 qrz(vdf)=qrz(vdf)+sh(3,inl)*q(vdf,inl)      
			 qrs(wdf)=qrs(wdf)+sh(0,inl)*q(wdf,inl)      
			 qrx(wdf)=qrx(wdf)+sh(1,inl)*q(wdf,inl)      
			 qry(wdf)=qry(wdf)+sh(2,inl)*q(wdf,inl)      
			 qrz(wdf)=qrz(wdf)+sh(3,inl)*q(wdf,inl)      
	    end do
!.....   calculate d(delta_v)/dt, delta_p, and d(delta_p)/dxj
	    do inl=1,nen
			 qrt(udf)=qrt(udf)+sh(0,inl)*q(udf,inl)*dtinv
			 qrt(vdf)=qrt(vdf)+sh(0,inl)*q(vdf,inl)*dtinv
			 qrt(wdf)=qrt(wdf)+sh(0,inl)*q(wdf,inl)*dtinv
			 qrs(pdf)=qrs(pdf)+sh(0,inl)*q(pdf,inl)      
			 qrx(pdf)=qrx(pdf)+sh(1,inl)*q(pdf,inl)      
			 qry(pdf)=qry(pdf)+sh(2,inl)*q(pdf,inl)      
			 qrz(pdf)=qrz(pdf)+sh(3,inl)*q(pdf,inl)      
	    enddo

!....    reset u=x_velocity, v=y_velocity, w=z_velocity, pp=pressure
	    u = drs(udf)
	    v = drs(vdf)
	    w = drs(wdf)
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
	    nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2  &
	                                            +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2  &
	                                            +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)
	    mu = mu + nu*ro                   

!.....   calculate the delta of residuals
	    res_c = 0.0
	    do inl=1,nen
		   res_c = res_c+(sh(xsd,inl)*q(udf,inl)  &
	                    + sh(ysd,inl)*q(vdf,inl)  &
	                    + sh(zsd,inl)*q(wdf,inl))/alpha
	    enddo

	    do isd = 1, nsd
		   res_a(isd)=ro*(qrt(isd)+u*qrx(isd)+v*qry(isd)+w*qrz(isd))
        enddo
        res_t(xsd) = qrx(pdf)/alpha + res_a(xsd) 
        res_t(ysd) = qry(pdf)/alpha + res_a(ysd) 
        res_t(zsd) = qrz(pdf)/alpha + res_a(zsd)

!....    calculate the stabilization parameters, taum and tauc
	    vel  = sqrt(u*u+v*v+w*w)
	    ree  = vel*hg/mu/12.0
	    if(steady.or.(.not.taudt)) then
		   taum = 1.0/sqrt((2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
	    else
		   taum = 1.0/sqrt((2.0/dt)**2+(2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
	    endif
	    taum = delta(1)* taum
	    tauc = delta(2)*hg*vel
	    if(ree.lt.1.0) tauc = tauc*ree

!.....   Density optimization
	    taum = taum/ro
	    tauc = tauc*ro 

	    do inl=1,nen
		   ph(0,inl) = sh(0,inl)*eft0
		   ph(1,inl) = sh(1,inl)*eft0
		   ph(2,inl) = sh(2,inl)*eft0
		   ph(3,inl) = sh(3,inl)*eft0
	    enddo
	      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.....   Galerkin Terms (Look at notes)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.....   calculate stress and pressure terms
	    txx = mu*(qrx(udf)+qrx(udf))
	    tyx = mu*(qry(udf)+qrx(vdf))
	    tzx = mu*(qrz(udf)+qrx(wdf))
	    txy = mu*(qrx(vdf)+qry(udf))
	    tyy = mu*(qry(vdf)+qry(vdf))
	    tzy = mu*(qrz(vdf)+qry(wdf))
	    txz = mu*(qrx(wdf)+qrz(udf))
	    tyz = mu*(qry(wdf)+qrz(vdf))
	    tzz = mu*(qrz(wdf)+qrz(wdf))
	    prs_t(udf) = res_t(udf) * taum
	    prs_t(vdf) = res_t(vdf) * taum
	    prs_t(wdf) = res_t(wdf) * taum
	    prs_c      = res_c*tauc
	      

!.... assemble the delta-residuals at nodes
	    do inl=1,nen
		   node = ien(inl,ie)
		   temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))
!.....      Continuty Equation
		   p(pdf,node) = p(pdf,node)+ph(0,inl)*res_c

!.....      Momentum Equation (Euler Residual)
		   p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
		   p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
		   p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)

!.....      Momentum Equation (C: Diffusion terms)
		   p(udf,node) = p(udf,node) - ph(xsd,inl)* pp + ph(xsd,inl)*txx + ph(ysd,inl)*tyx + ph(zsd,inl)*tzx 
		   p(vdf,node) = p(vdf,node) - ph(ysd,inl)* pp + ph(xsd,inl)*txy + ph(ysd,inl)*tyy + ph(zsd,inl)*tzy 
		   p(wdf,node) = p(wdf,node) - ph(zsd,inl)* pp + ph(xsd,inl)*txz + ph(ysd,inl)*tyz + ph(zsd,inl)*tzz 

!.....      Stablization with Tau_moment
		   p(pdf,node) = p(pdf,node) + ph(xsd,inl)*prs_t(udf) + ph(ysd,inl)*prs_t(vdf) + ph(zsd,inl)*prs_t(wdf)
		   p(udf,node) = p(udf,node) +prs_t(udf)*temp
		   p(vdf,node) = p(vdf,node) +prs_t(vdf)*temp
		   p(wdf,node) = p(wdf,node) +prs_t(wdf)*temp

!.....      Stablization with Tau_cont    
		   p(udf,node) = p(udf,node) + ph(xsd,inl)*prs_c
		   p(vdf,node) = p(vdf,node) + ph(ysd,inl)*prs_c
		   p(wdf,node) = p(wdf,node) + ph(zsd,inl)*prs_c

	    enddo


	    if(stokes) goto 500
!.....   Non-linear Terms (from momentum eq.) 

	       do isd = 1, nsd
	  	      res_a(isd)=ro*(qrs(1)*drx(isd)+qrs(2)*dry(isd)+qrs(3)*drz(isd))
	       enddo

	       do inl=1,nen
		      node = ien(inl,ie)
		      temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))

!.......    Additional nonlinear term form Galerkin
		      p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
		      p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
		      p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)

!.......    Additional nonlinear term form SUPG
		      p(pdf,node) = p(pdf,node) + taum*(ph(xsd,inl)*res_a(udf)+ph(ysd,inl)*res_a(vdf)+ph(zsd,inl)*res_a(wdf))

		      p(udf,node) = p(udf,node) + res_a(udf)*temp*taum
		      p(vdf,node) = p(vdf,node) + res_a(vdf)*temp*taum
		      p(wdf,node) = p(wdf,node) + res_a(wdf)*temp*taum
	       enddo

 500	continue 
     enddo  ! end of quad loop
  enddo  ! end of element loop

  return
end subroutine blockgmres
