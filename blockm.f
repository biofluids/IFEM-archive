c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       S. K. ALIABADI
c       modified for linear elastic equation, Lucy Zhang 4/22/99
c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine blockm(xloc, eloc, dloc, w, p, ien)

	implicit none
	include "global.h"

	integer ien(nen,nec)
	real* 8 xloc(nsd,nn_loc), dloc(nsd,nn_loc), eloc(nsd,nn_loc)
	real* 8 x(nsdpad,nenpad), d(nsdpad,nenpad), e(nsdpad,nenpad)
	real* 8 p(nsd,nn_loc)

	real* 8 eft0,det,eft1
	real* 8 sh(0:nsdpad,nenpad),ph(0:nsdpad,nenpad)
	real* 8 xr(nsdpad,nsdpad), cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	real* 8 drx(nsdpad),dry(nsdpad),drz(nsdpad)
	real* 8 ttt,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz

	real* 8 mu,la
	real* 8 erx,ery,erz,ers

c....   stiffness matrix
	real* 8 w(nsd,nn_loc)  
	integer inl, ie, isd, iq, node


	mu = 1.0
	la = mu * landa_over_mu
        do ie=1,nec 
	   do inl=1,nen
	      do isd=1,nsd
c		 e(isd,inl) = eloc(isd,ien(inl,ie))
		 x(isd,inl) = xloc(isd,ien(inl,ie))
		 d(isd,inl) = dloc(isd,ien(inl,ie))
	      enddo
	   enddo

	   do iq=1,nquad
	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if

	      eft0 = abs(det) * wq(iq)

c.... no jacobian
c	      eft0 = wq(iq)   
c	      write(*,*) 'eft0,eft1=',eft0,eft1
c............ Calculate local Stiffness Matrix k (kd=f)
	      do inl=1,nen
		 txx = sh(1,inl)**2
		 tyy = sh(2,inl)**2
		 tzz = sh(3,inl)**2
		 ttt = txx + tyy + tzz 
		 node = ien(inl,ie)
		 w(xsd,node)=w(xsd,node)+mu*(ttt+txx)*eft0
		 w(ysd,node)=w(ysd,node)+mu*(ttt+tyy)*eft0
		 w(zsd,node)=w(zsd,node)+mu*(ttt+tzz)*eft0  
		 
		 w(xsd,node)=w(xsd,node)+la*txx*eft0
		 w(ysd,node)=w(ysd,node)+la*tyy*eft0
		 w(zsd,node)=w(zsd,node)+la*tzz*eft0                  
	      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        CALCULATE kd-f
c                  in this problem f=0
c        ONLY CALCULATE  kd in local coordinates
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.........  initialize variables
	      erx = 0.0
	      ery = 0.0
	      erz = 0.0

	      do isd = 1,nsd
		 drx(isd) = 0.0
		 dry(isd) = 0.0
		 drz(isd) = 0.0
	      enddo

	      do inl=1,nen
c		 erx = erx + sh(0,inl)*e(1,inl)
c		 ery = ery + sh(0,inl)*e(2,inl)
c		 erz = erz + sh(0,inl)*e(3,inl)

c............... calculate the first derivative
		 do isd=1,nsd
		    drx(isd)=drx(isd)+sh(1,inl)*d(isd,inl)      
		    dry(isd)=dry(isd)+sh(2,inl)*d(isd,inl)      
		    drz(isd)=drz(isd)+sh(3,inl)*d(isd,inl)
		 enddo
	      end do
	    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      do inl=1,nen
		 ph(1,inl) = sh(1,inl)*eft0
		 ph(2,inl) = sh(2,inl)*eft0
		 ph(3,inl) = sh(3,inl)*eft0
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
       
	      do inl=1,nen
		 node=ien(inl,ie)
c		 write(*,*) node
c.....Elastic Equation (calculate residual: r=kd=p)
		 p(xsd,node) = p(xsd,node) -
	1	      ph(xsd,inl) * ttt -
	2	      ph(xsd,inl) * txx -
	3	      ph(ysd,inl) * tyx -
	4	      ph(zsd,inl) * tzx  
c	5	     +sh(xsd,inl) * erx
		 p(ysd,node) = p(ysd,node) -
	1	      ph(ysd,inl) * ttt -
	2	      ph(xsd,inl) * txy -
	3	      ph(ysd,inl) * tyy -
	4	      ph(zsd,inl) * tzy 
c	5	     +sh(ysd,inl) * ery
		 p(zsd,node) = p(zsd,node) -
	1	      ph(zsd,inl) * ttt -
	2	      ph(xsd,inl) * txz -
	3	      ph(ysd,inl) * tyz -
	4	      ph(zsd,inl) * tzz 
c	5	     +sh(zsd,inl) * erz


	      enddo
	   enddo
	enddo
      return
      end





