!c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c	S. Aliabadi                                                          c
!c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine blockgmresm(xloc,dloc,p,ien,jac, &
			ne_local,ien_local)
	use fluid_variables, only: nsd,nn,nen,ne,nquad,wq,sq
	use ale_variables
	use global_constants
	implicit none

	integer ien(nen,nen)
	real* 8 xloc(nsd,nn)
	real* 8 dloc(nsd,nn)
	real* 8 p(nsd,nn)
	real* 8 jac(nquad,ne)

	real* 8 x(nsd,nen)
	real* 8 d(nsd,nen)

	real* 8 eft0,det,eft1
	real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
	real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

	real* 8 drx(nsd),dry(nsd),drz(nsd)
	real* 8 ttt,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz
        real* 8 mu,la
!c	real* 8 mu,la,landa_over_mu
	integer inl, ie, isd, idf, iq, node
!============================
! MPI varibalbes
  integer ne_local ! # of element on each processor
  integer ien_local(ne_local) ! subregion-->wholeregion element index
  integer ie_local ! loop parameter


	mu = 1.0
	la = mu * landa_over_mu
        do ie_local=1,ne_local
	   ie=ien_local(ie_local) 
	   do inl=1,nen
	      do isd=1,nsd
		 x(isd,inl) = xloc(isd,ien(inl,ie))
		 d(isd,inl) = dloc(isd,ien(inl,ie))
	      enddo
	   enddo
	   do iq=1,nquad
	if (nsd == 3) then !3-D case
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




!c	      eft0 = wq(iq)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        CALCULATE k*delta(d)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.........  initialize variables
	      do isd = 1,nsd
		 drx(isd) = 0.0
		 dry(isd) = 0.0
        if (nsd == 3) then ! To be justified for both  2, 3 - D cases
		 drz(isd) = 0.0
	end if
	      enddo
	      
	      do inl=1,nen
!c............... calculate the first derivative
		 do isd=1,nsd
		    drx(isd)=drx(isd)+sh(1,inl)*d(isd,inl)      
		    dry(isd)=dry(isd)+sh(2,inl)*d(isd,inl) 
        if (nsd == 3) then ! To be justified for both  2, 3 - D cases      
		    drz(isd)=drz(isd)+sh(3,inl)*d(isd,inl)    
	end if
		 enddo
	      enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!c.....Elastic Equation (calculate k*(delta(d))=p)
		 p(xsd,node) = p(xsd,node) + &
		      ph(xsd,inl) * ttt + &
		      ph(xsd,inl) * txx + &
		      ph(ysd,inl) * tyx + &
		      ph(zsd,inl) * tzx 
		 p(ysd,node) = p(ysd,node) + &
		      ph(ysd,inl) * ttt + &
		      ph(xsd,inl) * txy + &
		      ph(ysd,inl) * tyy + &
		      ph(zsd,inl) * tzy 
		 p(zsd,node) = p(zsd,node) + &
		      ph(zsd,inl) * ttt + &
		      ph(xsd,inl) * txz + &
		      ph(ysd,inl) * tyz + &
		      ph(zsd,inl) * tzz 
	      enddo
	end if
        if (nsd == 2) then ! 2-D case, plain strain model
           do inl=1,nen
		node=ien(inl,ie)
                p(xsd,node)=p(xsd,node) + &
                        ph(xsd,inl)*(la+2*mu)*drx(xsd) + &
                        ph(ysd,inl)*mu*dry(xsd) + &
                        ph(xsd,inl)*la*dry(ysd) + &
                        ph(ysd,inl)*mu*drx(ysd)

                p(ysd,node)=p(ysd,node) + &
                        ph(ysd,inl)*la*drx(xsd) + &
                        ph(xsd,inl)*mu*dry(xsd) + &
                        ph(ysd,inl)*(la+2*mu)*dry(ysd) + &
                        ph(xsd,inl)*mu*drx(ysd)
          end do
        end if
	   enddo ! Gausian point loop
	enddo ! element loop
      return
      end



