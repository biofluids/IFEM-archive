!
!  energy_solid.f90
!  
!
!  Created by Zhang on 9/3/08.
!  Copyright 2008 __MyCompanyName__. All rights reserved.
!

       subroutine energy_solid(x_solid,nsd,nen,ien,nn_solid,ne_solid,x_solid_init,v_solid)
! x_solid: the potision of the solid domain at current time step
! nsd: dimension of the problem
! nen: number of nodes per element of solid mesh
! ien: connectiviy of the solid mesh
! nn_solid: number of nodes of solid mesh
! ne_solid: number of elements of solid mesh
! e_area: the area of every element
! rhs: (det(F_{ij}))^{-1}-1 right hand side of the equations 

use r_common, only: rc1,rc2,rk, density_solid, xmg, young_mod, Poisson
! constants for Mooney-Rivlin model
! gravity coefficients in all direction
      integer nsd
	  integer nn_solid
	  integer ne_solid
	  integer ien(ne_solid,nen)
      real(8) x_solid(nsd,nn_solid)
	  real(8) x_solid_init(nsd,nn_solid)
!=============================================================
! variable for E_int
      real(8) eint
! internal energy for whole solid domain
      real(8) eint_e(ne_solid)	
! internal energy per element	  
	  real(8) j1
	  real(8) j2
	  real(8) j3
! tensor invariance variables for c_ij=F_im F_mj	  	  
	  integer isd
	  integer jsd
	  integer ksd
! loop variable for nsd	  
	  real(8) cij(3,3)
	  real(8) cij2(3,3)
          real(8) detcij
! C_ij= F_im F_mj	  
! C_ijj2= C_im C_mj
! detcij = det{c_ij}
	  
	  
	  
	  
	  
!=============================================================


!============================================
! variable for EPt & EKt
	   real(8) delta_rho
	   ! density difference between fluid and solid: rho_solid - rho_fluid
	   real(8) x_gauss(nsd)
	   ! the coordinates of the integration point
	   real(8) ept_e(ne_solid)
	   ! the protential energy of solid (the artifical fluid has been substract) per element
	   real(8) ept
	   ! the prootential energy of the whole solid domain
	   real(8) v_solid(nsd,nn_solid)
	   ! grobal solid velocity
	   real(8) v(nsd,nen)
	   ! local solid velocity
	   real(8) v_gauss(nsd)
	   ! velocity at integration point
	   real(8) ekt_e(ne_solid)
	   ! kenectic energy per element (subcrabt the artificial fluid part) 1/2 \Delta \rho * v_s \dot v_s	
	   real(8) ekt
	   ! sum of ekt_e
!==============================================	  

! shape funciton variable       
	  integer iquad
	  integer nquad
	  real(8) xq(nsd,1)
	  real(8) wq
	  integer iq
	  real(8) sq(0:nsd,nen,1)
	  real(8) sh(0:nsd,nen)
	  real(8) xr(nsd,nsd)
	  real(8) cf(nsd,nsd)
	  real(8) sx(3,3)
	  real(8) det
! iquad: number of integration points always 1 here
! nquad: always 1
! sq: the same meaning of in the main program but only 1 in 3nd dimension here
! iq: always 1 here just for convenience
!===============================================================
! variables to calculate rhs & e_area
      real(8) e_area(ne_solid)
	  real(8) rhs(ne_solid)
	  real(8) y(nsd,nen)
	  real(8) rs(nsd)
	  real(8) xj(nsd,nsd)
	  real(8) xji(nsd,nsd)
	  real(8) detJ
	  real(8) toxj(nsd,nsd)
	  real(8) toxji(nsd,nsd)
	  real(8) todet
	  real(8) xto(nsd,nsd)
	  real(8) xot(nsd,nsd)
	  real(8) toc(nsd,nsd)
	  real(8) detf
!------------------------------------------------
          real(8) ge(nsd*2,ne_solid,1)
          real(8) cstr_element(nsd*2)
!------------------------------------------------

! det: determinant of Jacobian at current time step
! detf: det(F_{ij}
!=======================================
! Seting up left and right hand side variables
      integer ine
	  integer inen
	  integer index
	  integer jne
	  integer jnen
	  integer indey
	  integer inn
	  real(8) x(nsd,nen)
!	  real(8) y(nsd,nen)
	  
! crrm: the whole matirx of the left hand side 		
! crra: the whole matrix of the right hand side
! s_area: the sum area of the elements sounding every point
!====================================================
! varialbe to solve the linear system
      integer np
      integer mp
      integer error
      real(8) crra_cg(nsd*nn_solid)
!=============================================================
! First step: calculate the shape funciton in the parent domain, based on nsd & nen      
	  iquad=1
	  iq=1
          xq(:,:)=0.0d0
	  if (nsd==2) then
	     if (nen==3) then
 	!              write(*,*) 'I am here in volcorr 84', xq(:,:)
		     call quad2d3n(iquad,nquad,xq,wq,nsd,1)
         !             write(*,*) 'I am after call quad', xq(:,:)
		 else if (nen==4) then
		     call quad2d4n(iquad,nquad,xq,wq,nsd,1)
		 end if 
		 if(nen==3) then
		  sq(0,1,iq) = xq(1,iq)
		  sq(0,2,iq) = xq(2,iq)
		  sq(0,3,iq) = 1 - xq(1,iq) - xq(2,iq)
        else if (nen==4) then
		  sq(0,1,iq) = (1 - xq(1,iq)) * (1 - xq(2,iq)) / 4
		  sq(0,2,iq) = (1 + xq(1,iq)) * (1 - xq(2,iq)) / 4
		  sq(0,3,iq) = (1 + xq(1,iq)) * (1 + xq(2,iq)) / 4
		  sq(0,4,iq) = (1 - xq(1,iq)) * (1 + xq(2,iq)) / 4
		  
		  sq(1,1,iq) = - (1 - xq(2,iq)) / 4
		  sq(1,2,iq) = + (1 - xq(2,iq)) / 4
		  sq(1,3,iq) = + (1 + xq(2,iq)) / 4
		  sq(1,4,iq) = - (1 + xq(2,iq)) / 4
	
		  sq(2,1,iq) = - (1 - xq(1,iq)) / 4
		  sq(2,2,iq) = - (1 + xq(1,iq)) / 4
		  sq(2,3,iq) = + (1 + xq(1,iq)) / 4
		  sq(2,4,iq) = + (1 - xq(1,iq)) / 4

        endif
	  else if (nsd==3) then
	     if (nen==4) then 
		     call quad3d4n(iquad,nquad,xq,wq,nsd,1)
		 else if  (nen==8) then
		     call quad3d8n(iquad,nquad,xq,wq,nsd,1)
		 end if
	  
		 if(nen.eq.4) then
		  sq(0,1,iq) = xq(1,iq)
		  sq(0,2,iq) = xq(2,iq)
		  sq(0,3,iq) = xq(3,iq)
		  sq(0,4,iq) = 1 - xq(1,iq) - xq(2,iq) - xq(3,iq)
        else
		  sq(0,1,iq) = (1 - xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,2,iq) = (1 + xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,3,iq) = (1 + xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,4,iq) = (1 - xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,5,iq) = (1 - xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,6,iq) = (1 + xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,7,iq) = (1 + xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,8,iq) = (1 - xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,1,iq) = - (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,2,iq) = + (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,3,iq) = + (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,4,iq) = - (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,5,iq) = - (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,6,iq) = + (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,7,iq) = + (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,8,iq) = - (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,1,iq) = - (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,2,iq) = - (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,3,iq) = + (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,4,iq) = + (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,5,iq) = - (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,6,iq) = - (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,7,iq) = + (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,8,iq) = + (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(3,1,iq) = - (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,2,iq) = - (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,3,iq) = - (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,4,iq) = - (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,5,iq) = + (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,6,iq) = + (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,7,iq) = + (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,8,iq) = + (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
        endif
	  end if

!==================================================================
! Second step: Internal energy due to the deformation

! eint_e: internal energy interal per element = c1*(J1-3)+c2*(J2-3)+\kappa /2 *(J3-1)^2
          eint_e(:)=0
          eint=0
          ge(:,:,:)=0
          cstr_element(:)=0
          do ine=1,ne_solid
	     do inen=1,nen
		    y(1:nsd,inen)=x_solid(1:nsd,ien(ine,inen)) ! current position
			x(1:nsd,inen)=x_solid_init(1:nsd,ien(ine,inen)) ! initial position
		 end do
		 rs(1:nsd)=xq(1:nsd,1)
		 call r_element(rs)
		 call r_jacob(y,xj,xji,detJ)
		 call r_jacob(x,toxj,toxji,todet)
		 call r_bdpd_curr(xji)
		 call r_bdpd_init(toxji)
		 call r_stoxc(xto,xot,xj,xji,toxj,toxji,toc,ine)
		 e_area(ine)=wq*detJ ! Jacobbi determinant times weight
		 call determinant(xto,nsd,nsd,detf)
                 call r_sstrain(toc,xto,iq,ine,ge)
               !  write(*,*) ge(:,:,:)
                 call r_spiola_elastic(det,xot,ge,iq,ine,cstr_element)
                 eint_e(ine)= (1.0d0/(2*young_mod))*((cstr_element(1)+cstr_element(2))**2)- &
                (2.0d0*(1+Poisson)/young_mod)*(cstr_element(1)*cstr_element(2)-cstr_element(3)**2)

                 eint_e(ine)=eint_e(ine)*e_area(ine)
!=======================================================================
            ! For hyperelastic model
             !    cij(:,:)=0
             !    cij2(:,:)=0
		! do isd=1,nsd
		! do jsd=1,nsd
		! do ksd=1,nsd
		!	cij(isd,jsd)=cij(isd,jsd)+xto(ksd,isd)*xto(ksd,jsd)
                  !      write(*,*) 'cij', cij(isd,jsd) ! C_ij
	       !  end do
		! end do
		! end do
		! if (nsd==2) then
                !    cij(1,3)=0
                !    cij(2,3)=0
                !    cij(3,1)=0
                !    cij(3,2)=0
                !    cij(3,3)=1
               ! end if
	    ! do isd=1,3
		! do jsd=1,3
		! do ksd=1,3
		!	cij2(isd,jsd)=cij2(isd,jsd)+cij(isd,ksd)*cij(ksd,jsd) ! C_ij ^2
	    ! end do
            ! end do
	    ! end do
		!call determinant(cij,3,3,detcij) 
		! j3=detcij
		! j1=(cij(1,1)+cij(2,2)+cij(3,3))*(j3**(-0.3333))
		! j2=(j1*j1-cij2(1,1)-cij2(2,2)-cij2(3,3))*0.5d0*(j3**(-0.6666))
                ! j3=j3**(0.5)
           !      call determinant(cij,3,3,detcij)
		 
		 ! get J1 J2 J3
		! eint_e(ine)=(rc1*(j1-3)+rc2*(j2-3)+rk*(j3-1)*(j3-1)*0.5d0)*e_area(ine)
		 ! W*weight*D
	  end do
	     eint=sum(eint_e)
             write(*,*) sum(e_area)
       continue


!loop over element


!============================================
!step 3: calculate the protential energy for solid
      
	  ept_e(:)=0
	  ekt_e(:)=0
	  delta_rho=density_solid
          ept=0
          ekt=0
	  do ine=1,ne_solid
      ! assembling from global to local for coordinates and velocity
         do inen=1,nen
		    x(:,inen)=x_solid(1:nsd,ien(ine,inen))
			v(:,inen)=v_solid(1:nsd,ien(ine,inen))
		 end do  
	  ! get N_{j,i} 
      !...  calculate the shape function and the weight at quad point
		if (nsd==2) then
		    if (nen.eq.3) then !calculate shape function at quad point
			   include "sh2d3n.h"
			elseif (nen.eq.4) then
				include "sh2d4n.h"
			endif
		elseif (nsd==3) then
		    if (nen.eq.4) then !calculate shape function at quad point
			   include "sh3d4n.h"
			elseif (nen.eq.8) then
				include "sh3d8n.h"
			endif
		endif
                x_gauss(:)=0
                v_gauss(:)=0
		do isd=1,nsd
                do inen=1,nen
		x_gauss(isd)=x_gauss(isd)+x(isd,inen)*sh(0,inen)
		v_gauss(isd)=v_gauss(isd)+v(isd,inen)*sh(0,inen)
		end do
                end do
		
       do isd=1,nsd
		  ept_e(ine)=ept_e(ine)+xmg(isd)*x_gauss(isd)*delta_rho*e_area(ine)
	   end do
	   
	   do isd=1,nsd
	      ekt_e(ine)=ekt_e(ine)+v_gauss(isd)*v_gauss(isd)*delta_rho*e_area(ine)*0.5d0
	   end do
	   
	  end do
	   ept=sum(ept_e)
	   ekt=sum(ekt_e)
	 open(unit=92, file='EPts.txt', status='unknown')
	 open(unit=93, file='EKts.txt', status='unknown')
	 open(unit=94, file='Eints.txt', status='unknown')
	 write(92,*) ept
	 write(93,*) ekt
	 write(94,*) eint
       return
       end
