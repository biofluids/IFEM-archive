c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       S. K. ALIABADI
c       modified by L. Zhang for total ALE formulation, 7/21/99
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine block(xloc,dloc,doloc,floc,p,q,hk,ien,def,finv,jac,
	1    jaco, refvel, refvelo)
	implicit none
	include "global.h"

	integer ien(nen,ne)
	real* 8 xloc(nsd,nn),floc(nn)
	real* 8 dloc(ndf,nn),doloc(ndf,nn)
	real* 8 p(ndf,nn),q(ndf,nn),hk(ne)
	
	real* 8 x(nsdpad,nenpad),f(nenpad)
	real* 8 d(ndfpad,nenpad),do(ndfpad,nenpad)
	
	real* 8 eft0,det,effd,effm,effc
	real* 8 sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
	real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)
	
	real* 8 drt(ndfpad),drs(ndfpad)
	real* 8 drx(ndfpad),dry(ndfpad),drz(ndfpad)
	real* 8 drs_ref(nsdpad)
	real* 8 u,v,w,pp,fi,ug, wx,wy,wz

	real* 8 dtemp(nsd,nsd), t(nsd,nsd)
	real* 8 drx_hat,dry_hat,drz_hat
	real* 8 txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
	real* 8 pxx,pxy,pxz,pyx,pyy,pyz,pzx,pzy,pzz
	real* 8 hg,taum,tauc,vel,ree
	real* 8 res_c,res_a(nsdpad),res_t(nsdpad)
	real* 8 prs_c,prs_t(nsdpad)
	real* 8 mu,nu,ro,g(nsdpad)
	real* 8 tempu,tempv,tempw,temp
	real* 8 temprefu,temprefv,temprefw
	real* 8 dtinv,oma,ama
	integer inl, ie, isd, idf, iq, node
	real* 8 finv(nsd,nsd,nquad,ne),jac(nquad,ne),jaco(nquad,ne)
	real* 8 def(nsd,nsd,nquad,ne),tempf(nsd,nsd)
	real* 8 finv11,finv12,finv13,finv21,finv22,finv23,finv31
	real* 8 finv32,finv33, rodt
	real* 8 refvel(nsd,nquad,ne),refvelo(nsd,nquad,ne)
	real* 8 refv(nsd),refvo(nsd)
	real* 8 dp1(nsd),dp2(nsd,nsd,nsd),dp3(nsd),dfinv(nsd,nsd,nsd)
	real* 8 dp(nsd),dpp(nsd)
	integer i,j,k,l

	dtinv = 1.0/dt
	if(steady) dtinv = 0.0
	oma   = 1.0 -alpha
	ama   = 1.0 - oma

        do ie=1,ne 
	   do inl=1,nen
	      node=ien(inl,ie)
	      do isd=1,nsd
		 x(isd,inl) = xloc(isd,ien(inl,ie))
	      enddo
	      do idf=1,ndf
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
	      
	      eft0 = abs(det) * wq(iq)

	      do isd=1,nsd
		 refv(isd)=refvel(isd,iq,ie)
		 refvo(isd)=refvelo(isd,iq,ie)
		 if (its.eq.1) refvo(isd)=refv(isd)
	      enddo

	      do idf = 1,ndf
		 drs(idf) = 0.0
		 drx(idf) = 0.0
		 dry(idf) = 0.0
		 drz(idf) = 0.0
		 drt(idf) = 0.0
	      enddo
	      
	      do inl=1,nen
		 tempu = ama*d(udf,inl) + oma*do(udf,inl)
		 tempv = ama*d(vdf,inl) + oma*do(vdf,inl)
		 tempw = ama*d(wdf,inl) + oma*do(wdf,inl)
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
	      enddo
	      
	      do isd = 1,nsd
		 drs_ref(isd) = 0.0
	      enddo
	      
	      drs_ref(udf) = ama*refv(udf) + oma*refvo(udf)
	      drs_ref(vdf) = ama*refv(vdf) + oma*refvo(vdf)
	      drs_ref(wdf) = ama*refv(wdf) + oma*refvo(wdf)

	      do inl=1,nen
		 drt(udf)=drt(udf)+sh(0,inl)*(d(udf,inl)-do(udf,inl))*dtinv
		 drt(vdf)=drt(vdf)+sh(0,inl)*(d(vdf,inl)-do(vdf,inl))*dtinv
		 drt(wdf)=drt(wdf)+sh(0,inl)*(d(wdf,inl)-do(wdf,inl))*dtinv
		 drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)      
		 drx(pdf)=drx(pdf)+sh(1,inl)*d(pdf,inl)      
		 dry(pdf)=dry(pdf)+sh(2,inl)*d(pdf,inl)      
		 drz(pdf)=drz(pdf)+sh(3,inl)*d(pdf,inl)      
	      end do

	      u = drs(udf)
	      v = drs(vdf)
	      w = drs(wdf)
	      pp= drs(pdf)
	      
	      wx = drs_ref(udf)
	      wy = drs_ref(vdf)
	      wz = drs_ref(wdf)

	      fi = 0.0
	      do inl=1,nen
		 fi = fi + sh(0,inl)*f(inl)
	      enddo
cccccccccccccc
c	      mu = (1.0-fi)*vis_gas+fi*vis_liq
c	      ro = (1.0-fi)*den_gas+fi*den_liq
	      mu = vis_liq
	      ro = den_liq
	      ro = ro*jac(iq,ie) !ro(hat)
	      g(udf) = gravity(udf)
	      g(vdf) = gravity(vdf)
	      g(wdf) = gravity(wdf)
cccccccccccccc
	      
	      nu = delta(4)*turb_kappa**2*hg**2
	1	   * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2
	2	   +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2
	3	   +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)
	      mu = mu + nu*ro                   

	      res_c=0

	      do inl=1,nen
c		 res_c = res_c+sh(xsd,inl)*d(udf,inl)
c	1	      +sh(ysd,inl)*d(vdf,inl)
c	2	      +sh(zsd,inl)*d(wdf,inl)
		 do i=1,nsd
		    do j=1,nsd
		       res_c = res_c + sh(i,inl)*(ro*d(j,inl))*finv(i,j,iq,ie)
c		       res_c = res_c + sh(i,inl)*(d(j,inl))*finv(i,j,iq,ie)
		       res_c = res_c + 2*sh(i,inl)*jac(iq,ie)*d(j,inl)
	1		    *finv(i,j,iq,ie)
		    enddo
		 enddo
	      enddo
c		 if (ie.eq.10) then
c		    write(*,*) res_c
c		    stop
c		 endif
	      if(stokes) then
		 u = 0.0
		 v = 0.0
		 w = 0.0
	      endif
	      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      do isd = 1, nsd
		 res_a(isd) = ro*(drt(isd)+wx*drx(isd)+
	1	      wy*dry(isd)+wz*drz(isd)-g(isd))
	      enddo
c		 if (ie.eq.10) then
c		    write(*,*) wx,wy,wz
c		    write(*,*) res_a(1),res_a(2),res_a(3)
c		    stop
c		 endif

	      finv11 = finv(1,1,iq,ie)
	      finv12 = finv(1,2,iq,ie)
	      finv13 = finv(1,3,iq,ie)
	      finv21 = finv(2,1,iq,ie)
	      finv22 = finv(2,2,iq,ie)
	      finv23 = finv(2,3,iq,ie)
	      finv31 = finv(3,1,iq,ie)
	      finv32 = finv(3,2,iq,ie)
	      finv33 = finv(3,3,iq,ie)

	      dpp(1)=drx(pdf)
	      dpp(2)=dry(pdf)
	      dpp(3)=drz(pdf)

	      dtemp(1,1) = drx(udf)
	      dtemp(2,1) = drx(vdf)
	      dtemp(3,1) = drx(wdf)
	      dtemp(1,2) = dry(udf)
	      dtemp(2,2) = dry(vdf)
	      dtemp(3,2) = dry(wdf)
	      dtemp(1,3) = drz(udf)
	      dtemp(2,3) = drz(vdf)
	      dtemp(3,3) = drz(wdf)

	      do isd=1,nsd
		do i=1,nsd
		  do j=1,nsd
		    dfinv(i,j,isd)=0
		    do inl=1,nen
		      dfinv(i,j,isd)=dfinv(i,j,isd)+sh(isd,inl)*finv(i,j,iq,ie)
		    enddo
		   enddo
	         enddo
	      enddo

	      do i = 1,nsd
		 dp1(i) = 0
		 dp3(i) = 0
		 do j = 1,nsd
		    dp1(i) = dp1(i) + dpp(j)*jac(iq,ie)*finv(j,i,iq,ie)
c		    do k = 1,nsd
c		       dp2(k,i,j) = 0
c		       do l = 1,nsd		  
c			  dp2(k,i,j) = dp2(k,i,j) + dtemp(k,l)*dfinv(l,i,j)
c	1		        +dtemp(i,l)*dfinv(l,k,j)
c		       enddo
c		       dp3(i)=dp3(i)+mu*jac(iq,ie)*finv(j,k,iq,ie)*dp2(k,i,j)
c		    enddo
		 enddo
		 dp(i) = dp1(i) - dp3(i)
	      enddo

	      res_t(xsd) = dp(1) + res_a(xsd) 
	      res_t(ysd) = dp(2) + res_a(ysd) 
	      res_t(zsd) = dp(3) + res_a(zsd)

c	      res_t(xsd) = drx(pdf) + res_a(xsd) 
c	      res_t(ysd) = dry(pdf) + res_a(ysd) 
c	      res_t(zsd) = drz(pdf) + res_a(zsd)

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
	      prs_c      = res_c*tauc
	      do inl=1,nen
		 node=ien(inl,ie)
		 temp=ro*(wx*ph(xsd,inl)+wy*ph(ysd,inl)+wz*ph(zsd,inl))
c.....      Continuty Equation
		 p(pdf,node) = p(pdf,node)-ph(0,inl)*res_c

c.....      Momentum Equation (Euler Residual)
		 p(udf,node) = p(udf,node)-ph(0,inl)*res_a(udf)
		 p(vdf,node) = p(vdf,node)-ph(0,inl)*res_a(vdf)
		 p(wdf,node) = p(wdf,node)-ph(0,inl)*res_a(wdf)
c		 if (node.eq.10) then
c		    write(*,*) res_a(1),res_a(2),res_a(3)
c		    write(*,*) p(1,node),p(2,node),p(3,node)
c		    stop
c		 endif
c.....      Momentum Equation (C: Diffusion terms)
c.... Piola Kirchoff Stress
		 p(udf,node) = p(udf,node) +
	1	      ph(xsd,inl)* jac(iq,ie)*finv11 *pp - 
	2	      ph(xsd,inl)*pxx -
	3	      ph(ysd,inl)*pyx -
	4	      ph(zsd,inl)*pzx 
		 p(vdf,node) = p(vdf,node) +
	1	      ph(ysd,inl)* jac(iq,ie)*finv22 *pp -
	2	      ph(xsd,inl)*pxy -
	3	      ph(ysd,inl)*pyy -
	4	      ph(zsd,inl)*pzy 
		 p(wdf,node) = p(wdf,node) +
	1	      ph(zsd,inl)* jac(iq,ie)*finv33 *pp -
	2	      ph(xsd,inl)*pxz -
	3	      ph(ysd,inl)*pyz -
	4	      ph(zsd,inl)*pzz  

c.....      Stablization with Tau_moment
		 p(pdf,node) = p(pdf,node)
	1	      - ph(xsd,inl)*prs_t(udf)
	2	      - ph(ysd,inl)*prs_t(vdf)
	3	      - ph(zsd,inl)*prs_t(wdf)
		 p(udf,node) = p(udf,node) -prs_t(udf)*temp
		 p(vdf,node) = p(vdf,node) -prs_t(vdf)*temp
		 p(wdf,node) = p(wdf,node) -prs_t(wdf)*temp

c.....      Stablization with Tau_cont    
		 p(udf,node) = p(udf,node) - ph(xsd,inl)*prs_c
		 p(vdf,node) = p(vdf,node) - ph(ysd,inl)*prs_c
		 p(wdf,node) = p(wdf,node) - ph(zsd,inl)*prs_c
	      enddo
c	if (node.eq.10) write(*,*) p(1,node),p(2,node),p(3,node),p(4,node)  

c.....   Diagonal Preconditioner
	      
	      effd = mu*eft0*alpha
	      effm = taum*eft0
	      effc = tauc*eft0
	      
	      do inl=1,nen
		 
		 node=ien(inl,ie)
		 ug = ro*(u*sh(1,inl)+v*sh(2,inl)+w*sh(3,inl))
		 temp = alpha*ug + sh(0,inl)*dtinv*ro
		 
		 q(udf,node) = q(udf,node)+
	1	     (alpha*ro*sh(0,inl)*drx(1)+temp)*(sh(0,inl)*eft0+ug*effm)
		 q(vdf,node) = q(vdf,node)+
	1	     (alpha*ro*sh(0,inl)*dry(2)+temp)*(sh(0,inl)*eft0+ug*effm) 
		 q(wdf,node) = q(wdf,node)+
	1	     (alpha*ro*sh(0,inl)*drz(3)+temp)*(sh(0,inl)*eft0+ug*effm)
		 temp = sh(1,inl)**2+sh(2,inl)**2+sh(3,inl)**2
		 
		 q(udf,node) = q(udf,node)+(sh(1,inl)**2+temp)*effd
		 q(vdf,node) = q(vdf,node)+(sh(2,inl)**2+temp)*effd
		 q(wdf,node) = q(wdf,node)+(sh(3,inl)**2+temp)*effd
		 q(pdf,node) = q(pdf,node)+temp*effm
		 q(udf,node) = q(udf,node)+sh(1,inl)**2*effc
		 q(vdf,node) = q(vdf,node)+sh(2,inl)**2*effc
		 q(wdf,node) = q(wdf,node)+sh(3,inl)**2*effc
	      enddo
	   enddo
	enddo
	return
	end
	
