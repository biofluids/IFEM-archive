c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  S. Aliabadi                                                          c
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine blockgmres(xloc,dloc,doloc,qloc,p,hk,ien)

	  implicit none
	  include "global.h"

      integer ien(nen,nec)
      real* 8 xloc(nsd,nn_loc)
      real* 8 dloc(ndf,nn_loc),doloc(ndf,nn_loc)
      real* 8 p(ndf,nn_loc),qloc(ndf,nn_loc),hk(nec)

      real* 8 x(nsdpad,nenpad)
      real* 8 d(ndfpad,nenpad),do(ndfpad,nenpad),q(ndfpad,nenpad)

      real* 8 eft0,det,effd,effm,effc
      real* 8 sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	  real* 8 drs(ndfpad),qrt(ndfpad),qrs(ndfpad)
	  real* 8 drx(ndfpad),dry(ndfpad),drz(ndfpad)
	  real* 8 qrx(ndfpad),qry(ndfpad),qrz(ndfpad)
	  real* 8 u,v,w,pp,ug
      real* 8 txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
	  real* 8 hg,taum,tauc,vel,ree
	  real* 8 res_c,res_q,res_a(nsdpad),res_t(nsdpad)
	  real* 8 prs_c,prs_t(nsdpad)
      real* 8 mu,nu,ro
	  real* 8 tempu,tempv,tempw,temp
	  real* 8 dtinv,oma,ama
	  integer inl, ie, isd, idf, iq, node

	  dtinv = 1.0/dt/alpha
	  if(steady) dtinv = 0.0
      oma   = 1.0 -alpha
      ama   = 1.0 - oma

	  do ie=1,nec 

        do inl=1,nen
		  do isd=1,nsd
			x(isd,inl) = xloc(isd,ien(inl,ie))
		  enddo
		  do idf=1,ndf
			q(idf,inl) =  qloc(idf,ien(inl,ie))
			d(idf,inl) =  dloc(idf,ien(inl,ie))
            do(idf,inl) = doloc(idf,ien(inl,ie))
		  enddo
		enddo

		hg = hk(ie)

		do iq=1,nquad

		  if (nen.eq.4) then
			include "sh3d4n.h"
		  else if (nen.eq.8) then
			include "sh3d8n.h"
		  end if

		  eft0 = abs(det) * wq(iq) * alpha

		  do idf = 1,nsd
			drs(idf) = 0.0
			drx(idf) = 0.0
			dry(idf) = 0.0
			drz(idf) = 0.0
		  enddo
		  do inl=1,nen
			tempu= ama*d(udf,inl) + oma*do(udf,inl)
			tempv= ama*d(vdf,inl) + oma*do(vdf,inl)
			tempw= ama*d(wdf,inl) + oma*do(wdf,inl)
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
		  end do

		  do idf=1,ndf
			qrs(idf) = 0.0
			qrx(idf) = 0.0
			qry(idf) = 0.0
			qrz(idf) = 0.0
			qrt(idf) = 0.0
		  enddo
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
		  do inl=1,nen
			qrt(udf)=qrt(udf)+sh(0,inl)*q(udf,inl)*dtinv
			qrt(vdf)=qrt(vdf)+sh(0,inl)*q(vdf,inl)*dtinv
			qrt(wdf)=qrt(wdf)+sh(0,inl)*q(wdf,inl)*dtinv
			qrs(pdf)=qrs(pdf)+sh(0,inl)*q(pdf,inl)      
			qrx(pdf)=qrx(pdf)+sh(1,inl)*q(pdf,inl)      
			qry(pdf)=qry(pdf)+sh(2,inl)*q(pdf,inl)      
			qrz(pdf)=qrz(pdf)+sh(3,inl)*q(pdf,inl)      
		  end do

		  u = drs(udf)
		  v = drs(vdf)
		  w = drs(wdf)
		  pp= qrs(pdf)/alpha

ccccccccccc
		  mu = vis_liq
		  ro = den_liq
ccccccccccc

		  nu = delta(4)*turb_kappa**2*hg**2
	1		   * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2
	2		   +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2
	3		   +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)
		  mu = mu + nu*ro                   

		  res_c = 0.0
		  do inl=1,nen
			res_c = res_c+(sh(xsd,inl)*q(udf,inl)
	1			 + sh(ysd,inl)*q(vdf,inl)
	2			 + sh(zsd,inl)*q(wdf,inl))/alpha
		  enddo

		  if(stokes) then
			u = 0.0
			v = 0.0
			w = 0.0
		  endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		  do isd = 1, nsd
			res_a(isd)=ro*(qrt(isd)+u*qrx(isd)+v*qry(isd)+w*qrz(isd))
		  enddo
		  res_t(xsd) = qrx(pdf)/alpha + res_a(xsd) 
		  res_t(ysd) = qry(pdf)/alpha + res_a(ysd) 
		  res_t(zsd) = qrz(pdf)/alpha + res_a(zsd)

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

c.....Density optimization
		  taum = taum/ro
		  tauc = tauc*ro 

		  do inl=1,nen
			ph(0,inl) = sh(0,inl)*eft0
			ph(1,inl) = sh(1,inl)*eft0
			ph(2,inl) = sh(2,inl)*eft0
			ph(3,inl) = sh(3,inl)*eft0
		  enddo
		  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Galerkin Terms (Look at notes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
		  
		  do inl=1,nen

			node = ien(inl,ie)

			temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))

c.....Continuty Equation
			p(pdf,node) = p(pdf,node)+ph(0,inl)*res_c

c.....Momentum Equation (Euler Residual)
			p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
			p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
			p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)

c.....Momentum Equation (C: Diffusion terms)
			p(udf,node) = p(udf,node) -
	1			 ph(xsd,inl)* pp +
	2			 ph(xsd,inl)*txx +
	3			 ph(ysd,inl)*tyx +
	4			 ph(zsd,inl)*tzx 
			p(vdf,node) = p(vdf,node) -
	1			 ph(ysd,inl)* pp +
	2			 ph(xsd,inl)*txy +
	3			 ph(ysd,inl)*tyy +
	4			 ph(zsd,inl)*tzy 
			p(wdf,node) = p(wdf,node) -
	1			 ph(zsd,inl)* pp +
	2			 ph(xsd,inl)*txz +
	3			 ph(ysd,inl)*tyz +
	4			 ph(zsd,inl)*tzz 

c.....Stablization with Tau_moment
			p(pdf,node) = p(pdf,node)
	1			 + ph(xsd,inl)*prs_t(udf)
	2			 + ph(ysd,inl)*prs_t(vdf)
	3			 + ph(zsd,inl)*prs_t(wdf)
			p(udf,node) = p(udf,node) +prs_t(udf)*temp
			p(vdf,node) = p(vdf,node) +prs_t(vdf)*temp
			p(wdf,node) = p(wdf,node) +prs_t(wdf)*temp

c.....Stablization with Tau_cont    
			p(udf,node) = p(udf,node) + ph(xsd,inl)*prs_c
			p(vdf,node) = p(vdf,node) + ph(ysd,inl)*prs_c
			p(wdf,node) = p(wdf,node) + ph(zsd,inl)*prs_c

		  enddo

		  if(stokes) goto 500
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Non-linear Terms (from momentum eq.) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		  do isd = 1, nsd
			res_a(isd)=ro*(qrs(1)*drx(isd)+qrs(2)*dry(isd)+qrs(3)*drz(isd))
		  enddo

		  do inl=1,nen

			node = ien(inl,ie)

			temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))

c.......Additional nonlinear term form Galerkin
			p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
			p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
			p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)

c.......Additional nonlinear term form SUPG
			p(pdf,node) = p(pdf,node) + taum*(ph(xsd,inl)*res_a(udf)+
	1			 ph(ysd,inl)*res_a(vdf)+ph(zsd,inl)*res_a(wdf))

			p(udf,node) = p(udf,node) + res_a(udf)*temp*taum
			p(vdf,node) = p(vdf,node) + res_a(vdf)*temp*taum
			p(wdf,node) = p(wdf,node) + res_a(wdf)*temp*taum
		  enddo

 500	  continue 
		enddo
	  enddo

	  return
	  end