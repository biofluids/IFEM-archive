c
      subroutine extf_loadpar(istep,nstep,iLoadingType,vloadpar)
c
c*** input parameters:
c	istep
c	nstep
c	iLoadingType
c*** output parameter:
c       vloadpar:     between 0 and 1
c      
      implicit none
      integer istep,nstep
      integer iLoadingType
      real*8  vloadpar
c      
      integer istep1,istep2
      
      real*8  r,x,xd,yd,vk,angle,xc,yc,x1,y1,x2
      
      if  (iLoadingType.eq.0) then
      			! constant loading
         vloadpar=1.
         
      elseif (iLoadingType.eq.1) then      
      
	      		! static loading
         vloadpar=real(istep)/real(nstep)
         
      elseif (iLoadingType.eq.2) then
      
			! heaviside function loading
         if ( istep.eq.0 ) then      
            vloadpar=0.
         else
            vloadpar=1.
         endif
         
      elseif (iLoadingType.eq.3) then
         ! use 2 lines to simulate heaviside function loading
         
         istep1=nstep/8

         if (istep.le.istep1) then
            vloadpar=real(istep)/real(istep1)
         else
            vloadpar=1.
         endif            
         
      elseif (iLoadingType.eq.4) then
                        ! use circle to connect 2 lines of heaviside function
         
         
         x=real(istep)/real(nstep)

         r=0.05

		! (xd,yd) is the intersection point of line 1 and line 2
         xd=0.7
         yd=1.
         
         	! the slope of the line 1
         vk=yd/xd
         angle=atan(vk)
         
         
         yc=yd-r
         
         xc=(yc+r/ cos(angle) ) /vk
         
         x1=(xc+vk*yc)/(1.+vk*vk)
         y1=x1*vk
         
         x2=xc
c         
         if ( x.le.x1) then
            vloadpar=x*vk
         elseif ( x.le.x2) then
            vloadpar=yc+sqrt( r*r-(x-xc)**2 )
         else
            vloadpar=1.
         endif
      elseif ( iLoadingType.eq.5) then
         		! impulse impact
         istep1=1
         istep2=nstep/10
         if (istep.le.istep1) then
            vloadpar=real(istep)/real(istep1)
         elseif ( (istep.gt.istep1) 
     &         .and. (istep.le.istep2) ) then
            vloadpar=1.
         else
            vloadpar=0.
         endif               
      else
         write(*,*) 'In extf_loadpar, Wrong number of iLoadingType=',
     &              iLoadingType
         stop
      endif		! endcase of iLoadingType
      
      end	!ends extf_loadpar


      subroutine vel_loadpar(v_istep,nstep,iLoadingType,vloadpar)
c
c*** input parameters:
c	istep
c	nstep
c	iLoadingType
c*** output parameter:
c       vloadpar:     between 0 and 1
      
      implicit none
c      integer istep
      real*8 v_istep
      integer nstep
      integer iLoadingType
      real*8  vloadpar
      

      integer istep1,istep2
      
      real*8  r,x,xd,yd,vk,angle,xc,yc,x1,y1,x2
      
      if  (iLoadingType.eq.0) then
      			! constant loading
         vloadpar=1.
      elseif (iLoadingType.eq.1) then      
      
	      		! static loading
         vloadpar=real(v_istep)/real(nstep)
         
      elseif (iLoadingType.eq.2) then
      
			! heaviside function loading
         if ( v_istep.eq.0 ) then      
            vloadpar=0.
         else
            vloadpar=1.
         endif
         
      elseif (iLoadingType.eq.3) then
            		! use 2 lines to simulate heaviside function loading
         
         istep1=nstep/8

         if (v_istep.le.istep1) then
            vloadpar=real(v_istep)/real(istep1)
         else
            vloadpar=1.
         endif            
         
      elseif (iLoadingType.eq.4) then
                        ! use circle to connect 2 lines of heaviside function
         
         
         x=real(v_istep)/real(nstep)

         r=0.05

		! (xd,yd) is the intersection point of line 1 and line 2
         xd=0.7
         yd=1.
         
         	! the slope of the line 1
         vk=yd/xd
         angle=atan(vk)
         
         
         yc=yd-r
         
         xc=(yc+r/ cos(angle) ) /vk
         
         x1=(xc+vk*yc)/(1.+vk*vk)
         y1=x1*vk
         
         x2=xc
c         
         if ( x.le.x1) then
            vloadpar=x*vk
         elseif ( x.le.x2) then
            vloadpar=yc+sqrt( r*r-(x-xc)**2 )
         else
            vloadpar=1.
         endif
      elseif ( iLoadingType.eq.5) then
         		! impulse impact
         istep1=1
         istep2=nstep/10
         if (v_istep.le.istep1) then
            vloadpar=real(v_istep)/real(istep1)
         elseif ( (v_istep.gt.istep1) 
     &         .and. (v_istep.le.istep2) ) then
            vloadpar=1.
         else
            vloadpar=0.
         endif               
      else
         write(*,*) 'In vel_loadpar:Wrong number of iLoadingType=',
     &              iLoadingType
         stop
      endif		! endcase of iLoadingType
      
      end	!ends vel_loadpar

