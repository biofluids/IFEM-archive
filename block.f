c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  S. K. ALIABADI
c  - modified for RKPM by G. Wagner
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine block(xloc,shrk,don,doon,p,q,hk,ien,rng,
	1	   cnn,ncnn)

	  implicit none
	  include "global.h"

      integer ien(nen,nec),rng(neface,nec)
      integer cnn(maxconn,nqdc),ncnn(nqdc)
      real* 8 xloc(nsd,nn_loc)
	  real* 8 shrk(0:nsd,maxconn,nquad*nec)
      real* 8 don(ndf,nn_on),doon(ndf,nn_on)
      real* 8 p(ndf,nn_on),q(ndf,nn_on),hk(nec)

      real* 8 x(nsdpad,nenpad)
      real* 8 d(ndfpad,maxconn),do(ndfpad,maxconn)

      real* 8 eft0,det,effd,effm,effc
      real* 8 sh(0:nsdpad,nenpad)
      real* 8 ph(0:nsdpad,maxconn)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	  real* 8 drt(ndfpad),drs(ndfpad)
	  real* 8 drx(ndfpad),dry(ndfpad),drz(ndfpad)
	  real* 8 u,v,w,pp,fi,ug
      real* 8 txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
	  real* 8 hg,taum,tauc,vel,ree
	  real* 8 res_c,res_a(nsdpad),res_t(nsdpad)
	  real* 8 a1,a2,a3,b1,b2,b3,c1,c2,c3
	  real* 8 prs_c,prs_t(nsdpad)
      real* 8 mu,nu,ro,g(nsdpad)
	  real* 8 tempu,tempv,tempw,temp
	  real* 8 dtinv,oma,ama
	  real* 8 res_g,area,gg,coef
	  integer inl,ie,ieface,irng,isd,idf,inface,iq,node
	  integer qp,qb,node1,node2,node3,node4
	  integer ierr

c	  if (myid.eq.6) write (*,*) myid,"flag0!"

	  dtinv = 1.0/dt
	  if(steady) dtinv = 0.0
	  oma   = 1.0 - alpha
	  ama   = 1.0 - oma
c  coef = vis_liq/sqrt(hmin*hmin + hmax*hmax)
	  coef = 1e3

	  qp = 0
	  qb = 0

c	  if (myid.eq.6) write (*,*) myid,"flag1!"

	  do ie=1,nec 

		do inl=1,nen
		  do isd=1,nsd
			x(isd,inl) = xloc(isd,ien(inl,ie))
		  enddo
		enddo

		hg = hk(ie)

		do iq=1,nquad
		  qp = qp + 1

		  do inl=1,ncnn(qp)
			do idf = 1,ndf
			  d(idf,inl)  = don(idf,cnn(inl,qp))
			  do(idf,inl) = doon(idf,cnn(inl,qp))
			enddo
		  enddo

		  if (nen.eq.4) then
			include "sh3d4n.h"
		  else if (nen.eq.8) then
			include "sh3d8n.h"
		  end if

		  eft0 = abs(det) * wq(iq)

c  Now use RKPM shape functions

		  do idf = 1,ndf
			drs(idf) = 0.0
			drx(idf) = 0.0
			dry(idf) = 0.0
			drz(idf) = 0.0
			drt(idf) = 0.0
		  enddo

		  do inl=1,ncnn(qp)
c			stop
			tempu = ama*d(udf,inl) + oma*do(udf,inl)
			tempv = ama*d(vdf,inl) + oma*do(vdf,inl)
			tempw = ama*d(wdf,inl) + oma*do(wdf,inl)
			drs(udf)=drs(udf)+shrk(0,inl,qp)*tempu           
			drx(udf)=drx(udf)+shrk(1,inl,qp)*tempu           
			dry(udf)=dry(udf)+shrk(2,inl,qp)*tempu           
			drz(udf)=drz(udf)+shrk(3,inl,qp)*tempu           
			drs(vdf)=drs(vdf)+shrk(0,inl,qp)*tempv           
			drx(vdf)=drx(vdf)+shrk(1,inl,qp)*tempv           
			dry(vdf)=dry(vdf)+shrk(2,inl,qp)*tempv           
			drz(vdf)=drz(vdf)+shrk(3,inl,qp)*tempv           
			drs(wdf)=drs(wdf)+shrk(0,inl,qp)*tempw           
			drx(wdf)=drx(wdf)+shrk(1,inl,qp)*tempw           
			dry(wdf)=dry(wdf)+shrk(2,inl,qp)*tempw           
			drz(wdf)=drz(wdf)+shrk(3,inl,qp)*tempw           
			drt(udf)=drt(udf)+shrk(0,inl,qp)*(d(udf,inl)-do(udf,inl))*dtinv
			drt(vdf)=drt(vdf)+shrk(0,inl,qp)*(d(vdf,inl)-do(vdf,inl))*dtinv
			drt(wdf)=drt(wdf)+shrk(0,inl,qp)*(d(wdf,inl)-do(wdf,inl))*dtinv
			drs(pdf)=drs(pdf)+shrk(0,inl,qp)*d(pdf,inl)      
			drx(pdf)=drx(pdf)+shrk(1,inl,qp)*d(pdf,inl)      
			dry(pdf)=dry(pdf)+shrk(2,inl,qp)*d(pdf,inl)      
			drz(pdf)=drz(pdf)+shrk(3,inl,qp)*d(pdf,inl)      
		  end do

c		  if (myid.eq.6) write (*,*) myid,"flag2!"

		  u = drs(udf)
		  v = drs(vdf)
		  w = drs(wdf)
		  pp= drs(pdf)

		  mu = vis_liq
		  ro = den_liq
		  g(udf) = gravity(udf)
		  g(vdf) = gravity(vdf)
		  g(wdf) = gravity(wdf)

		  nu = delta(4)*turb_kappa**2*hg**2
	1		   * sqrt(2*drx(udf)**2+(dry(udf)+drx(vdf))**2
	2		   +2*dry(vdf)**2+(drz(udf)+drx(wdf))**2
	3		   +2*drz(wdf)**2+(drz(vdf)+dry(wdf))**2)
		  mu = mu + nu*ro                   


		  res_c = 0.0
		  do inl=1,ncnn(qp)
			res_c = res_c+ shrk(xsd,inl,qp)*d(udf,inl)
	1			         + shrk(ysd,inl,qp)*d(vdf,inl)
	2			         + shrk(zsd,inl,qp)*d(wdf,inl)
		  enddo

		  if(stokes) then
			u = 0.0
			v = 0.0
			w = 0.0
		  endif
	
c		  stop

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		  do isd = 1, nsd
			res_a(isd)=ro*(drt(isd)+u*drx(isd)+v*dry(isd)+w*drz(isd)-g(isd))
		  enddo

		  res_t(xsd) = drx(pdf) + res_a(xsd) 
		  res_t(ysd) = dry(pdf) + res_a(ysd) 
		  res_t(zsd) = drz(pdf) + res_a(zsd)

c		  if (myid.eq.6) write (*,*) myid,"flag3!"

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

		  do inl=1,ncnn(qp)
			ph(0,inl) = shrk(0,inl,qp)*eft0
			ph(1,inl) = shrk(1,inl,qp)*eft0
			ph(2,inl) = shrk(2,inl,qp)*eft0
			ph(3,inl) = shrk(3,inl,qp)*eft0
		  enddo

			
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Galerkin Terms (Look at notes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

		  txx = mu*(drx(udf)+drx(udf))
		  tyx = mu*(dry(udf)+drx(vdf))
		  tzx = mu*(drz(udf)+drx(wdf))
		  txy = mu*(drx(vdf)+dry(udf))
		  tyy = mu*(dry(vdf)+dry(vdf))
		  tzy = mu*(drz(vdf)+dry(wdf))
		  txz = mu*(drx(wdf)+drz(udf))
		  tyz = mu*(dry(wdf)+drz(vdf))
		  tzz = mu*(drz(wdf)+drz(wdf))
		  prs_t(udf) = res_t(udf) * taum
		  prs_t(vdf) = res_t(vdf) * taum
		  prs_t(wdf) = res_t(wdf) * taum
		  prs_c      = res_c*tauc

c		  if (myid.eq.6) write (*,*) myid,"flag4!"

		  do inl=1,ncnn(qp)
			
			node=cnn(inl,qp)

			temp=ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))

c			stop

c.....Continuty Equation
			p(pdf,node) = p(pdf,node)-ph(0,inl)*res_c

c.....Momentum Equation (Euler Residual)
			p(udf,node) = p(udf,node)-ph(0,inl)*res_a(udf)
			p(vdf,node) = p(vdf,node)-ph(0,inl)*res_a(vdf)
			p(wdf,node) = p(wdf,node)-ph(0,inl)*res_a(wdf)

c.....Momentum Equation (C: Diffusion terms)
			p(udf,node) = p(udf,node) +
	1			 ph(xsd,inl)* pp -
	2			 ph(xsd,inl)*txx -
	3			 ph(ysd,inl)*tyx -
	4			 ph(zsd,inl)*tzx 
			p(vdf,node) = p(vdf,node) +
	1			 ph(ysd,inl)* pp -
	2			 ph(xsd,inl)*txy -
	3			 ph(ysd,inl)*tyy -
	4			 ph(zsd,inl)*tzy 
			p(wdf,node) = p(wdf,node) +
	1			 ph(zsd,inl)* pp -
	2			 ph(xsd,inl)*txz -
	3			 ph(ysd,inl)*tyz -
	4			 ph(zsd,inl)*tzz 

c.....Stablization with Tau_moment
			p(pdf,node) = p(pdf,node)
	1			 - ph(xsd,inl)*prs_t(udf)
	2			 - ph(ysd,inl)*prs_t(vdf)
	3			 - ph(zsd,inl)*prs_t(wdf)
			p(udf,node) = p(udf,node) -prs_t(udf)*temp
			p(vdf,node) = p(vdf,node) -prs_t(vdf)*temp
			p(wdf,node) = p(wdf,node) -prs_t(wdf)*temp

c.....Stablization with Tau_cont    
			p(udf,node) = p(udf,node) - ph(xsd,inl)*prs_c
			p(vdf,node) = p(vdf,node) - ph(ysd,inl)*prs_c
			p(wdf,node) = p(wdf,node) - ph(zsd,inl)*prs_c
		  enddo

c		  stop
c
c.....Diagonal Preconditioner

c		  if (myid.eq.6) write (*,*) myid,"flag5!"

		  effd =   mu*eft0*alpha
		  effm = taum*eft0
		  effc = tauc*eft0

		  do inl=1,ncnn(qp)

			node=cnn(inl,qp)
			ug = ro*(u*shrk(1,inl,qp)+v*shrk(2,inl,qp)+w*shrk(3,inl,qp))
			temp = alpha*ug + shrk(0,inl,qp)*dtinv*ro

			q(udf,node) = q(udf,node)+
	1			 (alpha*ro*shrk(0,inl,qp)*drx(1)+temp)*(shrk(0,inl,qp)*eft0+ug*effm)

			q(vdf,node) = q(vdf,node)+
	1			 (alpha*ro*shrk(0,inl,qp)*dry(2)+temp)*(shrk(0,inl,qp)*eft0+ug*effm)

			q(wdf,node) = q(wdf,node)+
	1			 (alpha*ro*shrk(0,inl,qp)*drz(3)+temp)*(shrk(0,inl,qp)*eft0+ug*effm)
			
c			stop

			temp = shrk(1,inl,qp)**2+shrk(2,inl,qp)**2+shrk(3,inl,qp)**2

			q(udf,node) = q(udf,node)+(shrk(1,inl,qp)**2+temp)*effd
			q(vdf,node) = q(vdf,node)+(shrk(2,inl,qp)**2+temp)*effd
			q(wdf,node) = q(wdf,node)+(shrk(3,inl,qp)**2+temp)*effd
			q(pdf,node) = q(pdf,node)+temp*effm
			q(udf,node) = q(udf,node)+shrk(1,inl,qp)**2*effc
			q(vdf,node) = q(vdf,node)+shrk(2,inl,qp)**2*effc
			q(wdf,node) = q(wdf,node)+shrk(3,inl,qp)**2*effc
		  enddo

		enddo

	  enddo
	  
	  return
	  end




