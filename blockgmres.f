c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine blockgmres(xloc,dloc,doloc,floc,qloc,p,hk,ien,finv,jac,
	1    jaco,refvel, refvelo)

	implicit none
	include "global.h"
	
	integer ien(nen,nec)
	real* 8 xloc(nsd,nn_loc),floc(nn_loc)
	real* 8 dloc(ndf,nn_loc),doloc(ndf,nn_loc)
	real* 8 p(ndf,nn_loc),qloc(ndf,nn_loc),hk(nec)
	
	real* 8 x(nsdpad,nenpad),f(nenpad)
	real* 8 d(ndfpad,nenpad),do(ndfpad,nenpad),q(ndfpad,nenpad)
	
	real* 8 eft0,det,effd,effm,effc
	real* 8 sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
	real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	real* 8 drs(ndfpad),qrt(ndfpad),qrs(ndfpad),drs_ref(nsdpad)
	real* 8 drx(ndfpad),dry(ndfpad),drz(ndfpad)
	real* 8 qrx(ndfpad),qry(ndfpad),qrz(ndfpad)
	real* 8 qrx_hat,qry_hat,qrz_hat
	real* 8 u,v,w,pp,fi,ug,wx,wy,wz
	real* 8 txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
	real* 8 pxx,pxy,pxz,pyx,pyy,pyz,pzx,pzy,pzz
	real* 8 temprefu,temprefv,temprefw
	real* 8 finv(nsdpad,nsdpad,nquad,nec),jac(nquad,nec),jaco(nquad,nec)
	real* 8 t(nsdpad,nsdpad),dtemp(nsdpad,nsdpad)
	real* 8 hg,taum,tauc,vel,ree
	real* 8 res_c,res_q,res_a(nsdpad),res_t(nsdpad)
	real* 8 prs_c,prs_t(nsdpad)
	real* 8 mu,nu,ro,rodt
	real* 8 tempu,tempv,tempw,temp
	real* 8 dtinv,oma,ama
	integer inl, ie, isd, idf, iq, node,i,j,k
	real* 8 finv11,finv12,finv13,finv21,finv22,finv23,finv31
	real* 8 finv32,finv33
	real* 8 refvel(nsd,nquad,nec),refvelo(nsd,nquad,nec)
	real* 8 refv(nsd),refvo(nsd)
	
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
	      f(inl) = floc(ien(inl,ie))
	   enddo
	   hg = hk(ie)
	   do iq=1,nquad
	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if
	   
	      eft0 = abs(det) * wq(iq) * alpha

	      do isd=1,nsd
		 refv(isd)=refvel(isd,iq,ie)
		 refvo(isd)=refvelo(isd,iq,ie)
		 if(its.eq.1) refvo(isd)=refv(isd)
	      enddo

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

	      do isd = 1,nsd
		 drs_ref(isd) = 0.0
	      enddo

	      drs_ref(udf) = ama*refv(udf) + oma*refvo(udf)
	      drs_ref(vdf) = ama*refv(vdf) + oma*refvo(vdf)
	      drs_ref(wdf) = ama*refv(wdf) + oma*refvo(wdf)

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

	      wx = drs_ref(udf)
	      wy = drs_ref(vdf)
	      wz = drs_ref(wdf)

	      fi = 0.0
	      do inl=1,nen
		 fi = fi + sh(0,inl)*f(inl)
	      enddo
	      
ccccccccccc
c	      mu = (1.0-fi)*vis_gas+fi*vis_liq
c	      ro = (1.0-fi)*den_gas+fi*den_liq

	      mu = vis_liq
	      ro = den_liq
	      ro = ro*jac(iq,ie)
ccccccccccc
	     
	      nu = delta(4)*turb_kappa**2*hg**2
	1	   * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2
	2	   +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2
	3	   +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)
	      mu = mu + nu*ro                   

	      res_c=0
	      do inl=1,nen
c		 res_c = res_c+(sh(xsd,inl)*q(udf,inl)
c	1	      + sh(ysd,inl)*q(vdf,inl)
c	2	      + sh(zsd,inl)*q(wdf,inl))/alpha
		 do i=1,nsd
		    do j=1,nsd
		    res_c=res_c+sh(i,inl)*(ro*q(j,inl))*finv(i,j,iq,ie)/alpha
c		    res_c = res_c + sh(i,inl)*(q(j,inl))*finv(i,j,iq,ie)/alpha
		    res_c=res_c+2*sh(i,inl)*jac(iq,ie)*q(j,inl)*finv(i,j,iq,ie)/alpha
		    enddo
		 enddo
	      enddo
c		      if (ie.eq.10) then
c	      write(*,*) res_c
c	      stop
c	      endif      
	      if(stokes) then
		 u = 0.0
		 v = 0.0
		 w = 0.0
	      endif
	      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      do isd = 1, nsd
		 res_a(isd)=ro*(qrt(isd)+wx*qrx(isd)+wy*qry(isd)+wz*qrz(isd))
	      enddo
c	      if (ie.eq.10) then
c	      write(*,*) res_a(1),res_a(2),res_a(3)
c	      stop
c	      endif

	      finv11 = finv(1,1,iq,ie)
	      finv12 = finv(1,2,iq,ie)
	      finv13 = finv(1,3,iq,ie)
	      finv21 = finv(2,1,iq,ie)
	      finv22 = finv(2,2,iq,ie)
	      finv23 = finv(2,3,iq,ie)
	      finv31 = finv(3,1,iq,ie)
	      finv32 = finv(3,2,iq,ie)
	      finv33 = finv(3,3,iq,ie)
           qrx_hat=jac(iq,ie)*(finv11*qrx(pdf)+finv12*qry(pdf)+finv13*qrz(pdf))
           qry_hat=jac(iq,ie)*(finv21*qrx(pdf)+finv22*qry(pdf)+finv23*qrz(pdf))
           qrz_hat=jac(iq,ie)*(finv31*qrx(pdf)+finv32*qry(pdf)+finv33*qrz(pdf))

	      res_t(xsd) = qrx_hat/alpha + res_a(xsd) 
	      res_t(ysd) = qry_hat/alpha + res_a(ysd) 
	      res_t(zsd) = qrz_hat/alpha + res_a(zsd)
	      
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
	   
c.....   Density optimization
	      taum = taum/ro
	      tauc = tauc*ro 
	      do inl=1,nen
		 ph(0,inl) = sh(0,inl)*eft0
		 ph(1,inl) = sh(1,inl)*eft0
		 ph(2,inl) = sh(2,inl)*eft0
		 ph(3,inl) = sh(3,inl)*eft0
	      enddo
	      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....   Galerkin Terms (Look at notes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      dtemp(1,1) = qrx(udf)
	      dtemp(2,1) = qrx(vdf)
	      dtemp(3,1) = qrx(wdf)
	      dtemp(1,2) = qry(udf)
	      dtemp(2,2) = qry(vdf)
	      dtemp(3,2) = qry(wdf)
	      dtemp(1,3) = qrz(udf)
	      dtemp(2,3) = qrz(vdf)
	      dtemp(3,3) = qrz(wdf)
	      do i = 1,nsd
		 do j = 1,nsd
		    t(i,j) = 0
		    do k = 1,nsd
c		       t(i,j) = t(i,j) + finv(i,k,iq,ie)*dtemp(k,j)
c	1		               + finv(j,k,iq,ie)*dtemp(k,i)
		       t(i,j) = t(i,j) + dtemp(i,k)*finv(k,j,iq,ie)
	1		               + dtemp(j,k)*finv(k,i,iq,ie)
		    enddo
		 enddo
	      enddo

	      txx = mu * t(1,1)
	      tyx = mu * t(2,1)
	      tzx = mu * t(3,1)
	      txy = mu * t(1,2)
	      tyy = mu * t(2,2)
	      tzy = mu * t(3,2)
	      txz = mu * t(1,3)
	      tyz = mu * t(2,3)
	      tzz = mu * t(3,3)
	      
	      pxx = jac(iq,ie)*(finv11 *txx+finv12 *tyx+finv13 *tzx)
	      pyx = jac(iq,ie)*(finv21 *txx+finv22 *tyx+finv23 *tzx)
	      pzx = jac(iq,ie)*(finv31 *txx+finv32 *tyx+finv33 *tzx)
	      pxy = jac(iq,ie)*(finv11 *txy+finv12 *tyy+finv13 *tzy)
	      pyy = jac(iq,ie)*(finv21 *txy+finv22 *tyy+finv23 *tzy)
	      pzy = jac(iq,ie)*(finv31 *txy+finv32 *tyy+finv33 *tzy)
	      pxz = jac(iq,ie)*(finv11 *txz+finv12 *tyz+finv13 *tzz)
	      pyz = jac(iq,ie)*(finv21 *txz+finv22 *tyz+finv23 *tzz)
	      pzz = jac(iq,ie)*(finv31 *txz+finv32 *tyz+finv33 *tzz)

	      prs_t(udf) = res_t(udf) * taum
	      prs_t(vdf) = res_t(vdf) * taum
	      prs_t(wdf) = res_t(wdf) * taum
	      prs_c     = res_c*tauc

	      do inl=1,nen

		 node = ien(inl,ie)
		 temp=ro*(wx*ph(xsd,inl)+wy*ph(ysd,inl)+wz*ph(zsd,inl))
		 
c.....      Continuty Equation
		 p(pdf,node) = p(pdf,node)+ph(0,inl)*res_c
		 
c.....      Momentum Equation (Euler Residual)
		 p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
		 p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
		 p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)
c	   if (ie.eq.10) then
c	      write(*,*) p(1,node),p(2,node),p(3,node),p(4,node)
c	   stop
c	   endif
c.....      Momentum Equation (C: Diffusion terms)
c....       Piola Kirchoff Stress
		 p(udf,node) = p(udf,node) -
	1	      ph(xsd,inl)*jac(iq,ie)*finv11* pp +
	2	      ph(xsd,inl)*pxx +
	3	      ph(ysd,inl)*pyx +
	4	      ph(zsd,inl)*pzx 
		 p(vdf,node) = p(vdf,node) -
	1	      ph(ysd,inl)*jac(iq,ie)*finv22* pp +
	2	      ph(xsd,inl)*pxy +
	3	      ph(ysd,inl)*pyy +
	4	      ph(zsd,inl)*pzy 
		 p(wdf,node) = p(wdf,node) -
	1	      ph(zsd,inl)*jac(iq,ie)*finv33* pp +
	2	      ph(xsd,inl)*pxz +
	3	      ph(ysd,inl)*pyz +
	4	      ph(zsd,inl)*pzz 

c.....      Stablization with Tau_moment
		 p(pdf,node) = p(pdf,node)
	1	      + ph(xsd,inl)*prs_t(udf)
	2	      + ph(ysd,inl)*prs_t(vdf)
	3	      + ph(zsd,inl)*prs_t(wdf)
		 p(udf,node) = p(udf,node) +prs_t(udf)*temp
		 p(vdf,node) = p(vdf,node) +prs_t(vdf)*temp
		 p(wdf,node) = p(wdf,node) +prs_t(wdf)*temp
		 
c.....      Stablization with Tau_cont    
		 p(udf,node) = p(udf,node) + ph(xsd,inl)*prs_c
		 p(vdf,node) = p(vdf,node) + ph(ysd,inl)*prs_c
		 p(wdf,node) = p(wdf,node) + ph(zsd,inl)*prs_c
		 
	      enddo
c	   if (ie.eq.10) then
c	      write(*,*) p(1,node),p(2,node),p(3,node),p(4,node)
c	   stop
c	   endif
	      if(stokes) goto 500
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....   Non-linear Terms (from momentum eq.) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      do isd = 1, nsd
		res_a(isd)=ro*(qrs(1)*drx(isd)+qrs(2)*dry(isd)+qrs(3)*drz(isd))
	      enddo
	      
	      do inl=1,nen
		 
		 node = ien(inl,ie)		 
		 temp=ro*(wx*ph(xsd,inl)+wy*ph(ysd,inl)+wz*ph(zsd,inl))
		 
c.......    Additional nonlinear term form Galerkin
		 p(udf,node) = p(udf,node)+ph(0,inl)*res_a(udf)
		 p(vdf,node) = p(vdf,node)+ph(0,inl)*res_a(vdf)
		 p(wdf,node) = p(wdf,node)+ph(0,inl)*res_a(wdf)
		 
c.......    Additional nonlinear term form SUPG
		 p(pdf,node) = p(pdf,node) + taum*(ph(xsd,inl)*res_a(udf)+
	1	      ph(ysd,inl)*res_a(vdf)+ph(zsd,inl)*res_a(wdf)) 
		 p(udf,node) = p(udf,node) + res_a(udf)*temp*taum
		 p(vdf,node) = p(vdf,node) + res_a(vdf)*temp*taum
		 p(wdf,node) = p(wdf,node) + res_a(wdf)*temp*taum

	      enddo
	      
 500	      continue 
	   enddo
	enddo
	
	return
	end
