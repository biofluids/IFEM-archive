c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  S. K. ALIABADI
c  - modified for RKPM by G. Wagner
c  - modified for elasticity equation by L. Zhang
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine blockm(xloc,shrk,don,p,w,ien,rng,cnn,ncnn)

	implicit none
	include "global.h"

	integer ien(nen,nec),rng(neface,nec)
	integer cnn(maxconn,nqdc),ncnn(nqdc)
	real* 8 xloc(nsd,nn_loc)
	real* 8 shrk(0:nsd,maxconn,nquad*nec)
	real* 8 don(nsd,nn_on)
	real* 8 p(nsd,nn_on),w(nsd,nn_on)
	
	real* 8 x(nsdpad,nenpad)
	real* 8 d(nsdpad,maxconn)

	real* 8 eft0,det
	real* 8 sh(0:nsdpad,nenpad)
	real* 8 ph(0:nsdpad,maxconn)
	real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	real* 8 drx(nsdpad),dry(nsdpad),drz(nsdpad)
	real* 8 txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz,ttt
	real* 8 mu,la, landa_over_mu
	integer inl,ie,ieface,irng,isd,inface,iq,node
	integer ierr,qp

	qp=0
	mu = 1.0
	la = mu * landa_over_mu
	do ie = 1,nec 
	   do inl=1,nen
	      do isd=1,nsd
		 x(isd,inl) = xloc(isd,ien(inl,ie))
	      enddo
	   enddo

	   do iq = 1,nquad
	      qp = qp + 1

	      do inl = 1,ncnn(qp)
		 do isd = 1,nsd
		    d(isd,inl) = don(isd,cnn(inl,qp))
		 enddo
	      enddo

	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if

	      eft0 = abs(det) * wq(iq)

c......Calculate local Stiffness Matrix k (kd=f)
	      do inl = 1,ncnn(qp)
		 txx = shrk(1,inl,qp)**2
		 tyy = shrk(2,inl,qp)**2
		 tzz = shrk(3,inl,qp)**2
		 ttt = txx+tyy+tzz
		 node = cnn(inl,qp)
		 w(xsd,node) = w(xsd,node)+mu*(ttt+txx)*eft0
		 w(ysd,node) = w(ysd,node)+mu*(ttt+tyy)*eft0
		 w(zsd,node) = w(zsd,node)+mu*(ttt+tzz)*eft0
		 
		 w(xsd,node) = w(xsd,node)+la*txx*eft0
		 w(ysd,node) = w(ysd,node)+la*tyy*eft0
		 w(zsd,node) = w(zsd,node)+la*tzz*eft0
	      enddo

c  Now use RKPM shape functions

	      do isd = 1,nsd
		 drx(isd) = 0.0
		 dry(isd) = 0.0
		 drz(isd) = 0.0
	      enddo

	      do inl = 1,ncnn(qp)
		 do isd = 1,nsd
		    drx(isd) = drx(isd)+shrk(1,inl,qp)*d(isd,inl)           
		    dry(isd) = dry(isd)+shrk(2,inl,qp)*d(isd,inl)           
		    drz(isd) = drz(isd)+shrk(3,inl,qp)*d(isd,inl)           
		 enddo
	      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      do inl = 1,ncnn(qp)
		 ph(0,inl) = shrk(0,inl,qp)*eft0
		 ph(1,inl) = shrk(1,inl,qp)*eft0
		 ph(2,inl) = shrk(2,inl,qp)*eft0
		 ph(3,inl) = shrk(3,inl,qp)*eft0
	      enddo
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....Galerkin Terms (Look at notes)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      ttt = la*(drx(xsd)+dry(ysd)+drz(zsd))
	      txx = mu*(drx(xsd)+drx(xsd))
	      tyx = mu*(dry(xsd)+drx(ysd))
	      tzx = mu*(drz(xsd)+drx(zsd))
	      txy = mu*(drx(ysd)+dry(xsd))
	      tyy = mu*(dry(ysd)+dry(ysd))
	      tzy = mu*(drz(ysd)+dry(zsd))
	      txz = mu*(drx(zsd)+drz(xsd))
	      tyz = mu*(dry(zsd)+drz(ysd))
	      tzz = mu*(drz(zsd)+drz(zsd))

	      do inl = 1,ncnn(qp)	
		 node = cnn(inl,qp)

c.....Elasticity Equation (residual term:  r=kd=p)
		 p(xsd,node) = p(xsd,node) -
	1	      ph(xsd,inl)*ttt -
	2	      ph(xsd,inl)*txx -
	3	      ph(ysd,inl)*tyx -
	4	      ph(zsd,inl)*tzx 
		 p(ysd,node) = p(ysd,node) -
	1	      ph(ysd,inl)*ttt -
	2	      ph(xsd,inl)*txy -
	3	      ph(ysd,inl)*tyy -
	4	      ph(zsd,inl)*tzy 
		 p(zsd,node) = p(zsd,node) -
	1	      ph(zsd,inl)*ttt -
	2	      ph(xsd,inl)*txz -
	3	      ph(ysd,inl)*tyz -
	4	      ph(zsd,inl)*tzz 
		 
	      enddo
	   enddo
	enddo
	  
	return
	end




