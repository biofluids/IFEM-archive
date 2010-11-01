!c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       S. K. ALIABADI
!c       modified for linear elastic equation, Lucy Zhang 4/22/99
!c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine blockm(xloc, eloc, dloc, w, p, ien,jac)
	use fluid_variables, only: nsd,nen,ne,nn,nquad,wq,sq
	use global_constants
	use ale_variables
	implicit none

	integer ien(nen,ne)
	real* 8 xloc(nsd,nn), dloc(nsd,nn), eloc(nsd,nn)
	real* 8 x(nsd,nen), d(nsd,nen), e(nsd,nen)
	real* 8 p(nsd,nn)
	real*8 jac(nquad,ne)

	real* 8 eft0,det,eft1
	real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
	real* 8 xr(nsd,nsd), cf(nsd,nsd),sx(nsd,nsd)

	real* 8 drx(nsd),dry(nsd),drz(nsd)
	real* 8 ttt,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz

	real* 8 mu,la
	real* 8 erx,ery,erz,ers

!c....   stiffness matrix
	real* 8 w(nsd,nn)  
	integer inl, ie, isd, iq, node


	mu = 1.0
	la = mu * landa_over_mu
        do ie=1,ne 
	   do inl=1,nen
	      do isd=1,nsd
!c		 e(isd,inl) = eloc(isd,ien(inl,ie))
		 x(isd,inl) = xloc(isd,ien(inl,ie))
		 d(isd,inl) = dloc(isd,ien(inl,ie))
	      enddo
	   enddo


	   do iq=1,nquad
	if (nsd == 3) then ! 3-D case
	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if
	end if
	if (nsd == 2) then ! 2-D case
	      if (nen.eq.3) then !calculate shape function at quad point
                 include "sh2d3n.h"
              elseif (nen.eq.4) then
                 include "sh2d4n.h"
              endif
	end if
	      eft0 = abs(det) * wq(iq)

!=============================================
! Find out when the mesh is distored 
	if (jac(iq,ie) .lt. 0.0) then
	write(*,*) 'Moving mesh is distorted'
	stop
	end if
! D=D/J stiffer for larger deformation elements

	mu =1.0
	mu = 1.0/jac(iq,ie)
	la = mu*landa_over_mu

!c.... no jacobian
!c	      eft0 = wq(iq)   
!c............ Calculate local Stiffness Matrix k (kd=f)
! Use the diagonal to be the pre-conditioner, see nodes
	      do inl=1,nen
	if (nsd == 3) then ! 3-D case
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
	end if

	if (nsd == 2) then ! 2-D case
		node = ien(inl,ie)
		w(xsd,node)=w(xsd,node)+(sh(1,inl)**2)*(la+2*mu)+(sh(2,inl)**2)*mu
		w(ysd,node)=w(ysd,node)+(sh(1,inl)**2)*mu+(sh(2,inl)**2)*(la+2*mu) 
	end if                
	      enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        CALCULATE kd-f
!c                  in this problem f=0
!c        ONLY CALCULATE  kd in local coordinates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.........  initialize variables
	      erx = 0.0
	      ery = 0.0
	      erz = 0.0

	      do isd = 1,nsd
		 drx(isd) = 0.0
		 dry(isd) = 0.0
	if (nsd == 3) then ! To be justified for both  2, 3 - D cases
		 drz(isd) = 0.0
	end if 
	      enddo

	      do inl=1,nen
!c		 erx = erx + sh(0,inl)*e(1,inl)
!c		 ery = ery + sh(0,inl)*e(2,inl)
!c		 erz = erz + sh(0,inl)*e(3,inl)

!c............... calculate the first derivative
		 do isd=1,nsd
		    drx(isd)=drx(isd)+sh(1,inl)*d(isd,inl)      
		    dry(isd)=dry(isd)+sh(2,inl)*d(isd,inl)     
        if (nsd == 3) then ! To be justified for both  2, 3 - D cases 
		    drz(isd)=drz(isd)+sh(3,inl)*d(isd,inl)
	end if
		 enddo
	      end do
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      do inl=1,nen
		 ph(:,inl) = sh(:,inl)*eft0
!		 ph(1,inl) = sh(1,inl)*eft0
!		 ph(2,inl) = sh(2,inl)*eft0
!		 ph(3,inl) = sh(3,inl)*eft0
	      enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	if (nsd == 3) then ! 3-D case
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
!c.....Elastic Equation (calculate residual: r=kd=p)
		 p(xsd,node) = p(xsd,node) - &
		      ph(xsd,inl) * ttt - &
		      ph(xsd,inl) * txx - &
		      ph(ysd,inl) * tyx - &
		      ph(zsd,inl) * tzx  
!c		     +sh(xsd,inl) * erx
		 p(ysd,node) = p(ysd,node) - &
		      ph(ysd,inl) * ttt - &
		      ph(xsd,inl) * txy - &
		      ph(ysd,inl) * tyy - &
		      ph(zsd,inl) * tzy 
!c		     +sh(ysd,inl) * ery
		 p(zsd,node) = p(zsd,node) - &
		      ph(zsd,inl) * ttt - &
		      ph(xsd,inl) * txz - &
		      ph(ysd,inl) * tyz - &
		      ph(zsd,inl) * tzz   
!c		     +sh(zsd,inl) * erz

	      enddo
	end if 

	if (nsd == 2) then ! 2-D case, plain strain model
	   do inl=1,nen
		node=ien(inl,ie)
		p(xsd,node)=p(xsd,node) - &
			ph(xsd,inl)*(la+2*mu)*drx(xsd) - &
			ph(ysd,inl)*mu*dry(xsd) - &
			ph(xsd,inl)*la*dry(ysd) - &
			ph(ysd,inl)*mu*drx(ysd)

		p(ysd,node)=p(ysd,node) - &
			ph(ysd,inl)*la*drx(xsd) - &
			ph(xsd,inl)*mu*dry(xsd) - &
			ph(ysd,inl)*(la+2*mu)*dry(ysd) - &
			ph(xsd,inl)*mu*drx(ysd)
	  end do
	end if

	   enddo ! Gausian point loop
	enddo ! element loop
      return
      end





