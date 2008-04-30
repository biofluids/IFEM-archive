!
!  energy_fluid.f90
!  
!
!  Created by Zhang on 9/3/08.
!  Copyright 2008 __MyCompanyName__. All rights reserved.
!

       subroutine energy_fluid(x_fluid,v_fluid,ien)     
       use fluid_variables, only: vis_liq, den_liq, nn, nsd, nen, ne
	   ! obtain: sq--for shape function in parent domain;
	   !         vis_liq: fluid viscosity
	   !         den_liq: fluid density
	   
	   
	   
	   integer ien(nen,ne)
!========================================================	   
	   
	   ! mesh information for fluid domain
	   real(8) x_fluid(nsd,nn)
	   real(8) v_fluid(nsd,nn)
	   real(8) x(nsd,nen)
	   real(8) v(nsd,nen)
           real(8) v_gauss(nsd)
	   ! local and global coordinates and velocity for fluid
!=========================================================
	   
	   integer iq
	   integer nquad
	   integer iquad
	   real(8) wq
	   real(8) sq(0:nsd,nen,1)
           real(8) xq(nsd,1)
	   real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
	   real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
	   ! variable defined to get shape function
	   real(8) e_area
	   ! determinant of jacobbi times weight
	   real(8) det
	   ! determinant of jacobbi
!=========================================================	   
	   integer ine
	   integer inen
	   integer isd
	   integer jsd
	   ! loop variables
	   real(8) ekf_e(ne)
	   real(8) ekf
	   ! local and global kenactic energy for fluid
	    real(8) dvdx(nsd,nsd)
	   !spacial derivatives of fluid velocity 
	
		real(8) visrate_e(ne)
		real(8) visrate
		! local and global viscous disspation rate
		
	  
	   
	   
	   ! First step: calculate the shape funciton in the parent domain, based on nsd & nen      
	  iquad=1
	  iq=1
          xq(:,1)=0.0d0
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
	   
	   
	   
	   ekf_e(:)=0
	   visrate_e(:)=0
           visrate=0
           ekf=0
	   do ine=1,ne
	      do inen=1,nen
		  x(1:nsd,inen)=x_fluid(1:nsd,ien(inen,ine))
                  v(1:nsd,inen)=v_fluid(1:nsd,ien(inen,ine))
	       end do
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
        e_area=abs(det)*wq
                 v_gauss(:)=0
		 do inen=1,nen
		 v_gauss(1:nsd)=v_gauss(1:nsd)+v(1:nsd,inen)*sh(0,inen)
		 end do
		 do isd=1,nsd
		 ekf_e(ine)=ekf_e(ine)+v_gauss(isd)*v_gauss(isd)*den_liq*e_area*0.5d0
		 end do
		 dvdx(:,:)=0
		 do isd=1,nsd
		 do jsd=1,nsd
		 do inen=1,nen
		 dvdx(isd,jsd)=dvdx(isd,jsd)+v(isd,inen)*sh(jsd,inen)
		 ! isd---direction index for velocity
		 ! jsd---direction index for spacial derivative
		 end do
		 end do
		 end do
		 if (nsd==2) then
		 visrate_e(ine)=((dvdx(1,1)*dvdx(1,1)+dvdx(2,2)*dvdx(2,2))*2.0d0+(dvdx(2,1)+dvdx(1,2))*(dvdx(2,1)+dvdx(1,2)))*vis_liq*e_area
		 end if
	   end do
		 ekf=sum(ekf_e)
		 visrate=sum(visrate_e)
		
       open(unit=95, file='EKtf.txt',status='unknown')
       open(unit=96, file='Evis.txt',status='unknown')
       write(95,*) ekf
       write(96,*) visrate	

       
	   
	   
	   
	   return
       end
