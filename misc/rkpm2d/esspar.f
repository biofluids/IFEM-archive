c
      subroutine esspar(v_istep,nstep,dt,iLoadingType,
     &                  dispar,velpar,accpar)
c
c*** input parameters:
c
c    v_istep := istep + theta;
c    nstep   :  the total time step in the computations;
c
c    iLoadingType:  
c
c    0 : constant;
c    1 : Heaviside;
c    2 : linear; (natural proportional)
c    3 : bilinear (set the mark at int(nstep/20)  as default)
c
c
c
c*** output parameter:
c     
c  dispar:  Coefficient to determine the exact displacement on Vbc   
c  velpar:  Coefficeint to determine the exact velocity on Vbc 
c  accpar:  Coefficient to determine the exact acceleration on Vbc
c
c
c--------------------------------------------------
c
      implicit none
      integer nstep,iLoadingType 
      real*8  v_istep,dt
      real*8  dispar,velpar,accpar
      real*8  Tcr1,Tcr2,Treal
c      

      integer istep1
c
c
      if (iLoadingType.eq.0) then
c
c.............Constant velocity .............
c
	 accpar = 0.0
         velpar = 1.0
	 dispar = velpar * v_istep * dt
c
      elseif (iLoadingType .eq. 1) then      
c
c......... Heaviside function ...............
c
         if (v_istep .eq. 0) then      
	    accpar = 0.0
            velpar = 0.0
	    dispar = 0.0
         else
	    accpar = 0.0
            velpar = 1.0
	    dispar = v_istep * dt
         endif
c
      elseif (iLoadingType .eq. 2) then
c
c............. Proportional Loading..........
c      
	 accpar = 1.0/real(nstep*dt)
         velpar = real(v_istep)/real(nstep)
	 dispar = 0.5 * velpar * v_istep * dt
c      
      elseif (iLoadingType .eq. 3) then
c
c............Use 2 lines to simulate heaviside function loading.....
c
c       set the default
c
         istep1= int(nstep/20)

         if (v_istep .le. istep1) then
	    accpar = 1.0/real(istep1*dt)
            velpar = v_istep/real(istep1)
	    dispar = 0.5 * velpar * v_istep * dt
         else
	    accpar = 0.0
            velpar = 1.0
	    dispar = 0.5 * real(istep1) * dt 
     &             + (v_istep - real(istep1)) * dt 
         endif            
      elseif (iLoadingType .eq. 4) then
c
c............kalthoff loading.....
c
c       set the default
c
        Tcr1  = 0.50d-6
	Tcr2  = 47.0d-6
        Treal = v_istep * dt
c
	if(Treal .le. Tcr1) then
	    accpar = 1.0/Tcr1
            velpar = v_istep*dt/Tcr1
	    dispar = 0.5 * velpar * v_istep * dt
	elseif((Treal .ge. Tcr1) 
     &        .and. (Treal .le. Tcr2)) then
c
	    accpar = 0.0
            velpar = 1.0
	    dispar = 0.5 * Tcr1 + (Treal - Tcr1)
c
	elseif(Treal .gt. Tcr2) then
	    accpar = 0.0
            velpar = 0.0
	    dispar = 0.0
        endif
c
      elseif (iLoadingType .eq. 5) then
c
c............kalthoff loading.....
c
c       set the default
c
        Tcr1  = 0.50d-6
	Tcr2  = 42.0d-6
        Treal = v_istep * dt
c
	if(Treal .le. Tcr1) then
	    accpar = 1.0/Tcr1
            velpar = v_istep*dt/Tcr1
	    dispar = 0.5 * velpar * v_istep * dt
	elseif((Treal .ge. Tcr1) 
     &     .and. (Treal .le. Tcr2)) then
c
	    accpar = 0.0
            velpar = 1.0
	    dispar = 0.5 * Tcr1 + (Treal - Tcr1)
c
	elseif(Treal .gt. Tcr2) then
	    accpar = 0.0
            velpar = 0.0
	    dispar = 0.0
        endif
c
      else
         write(*,*) 'In vel_loadpar:Wrong number of iLoadingType=',
     &              iLoadingType
         stop
      endif		! endcase of iLoadingType
c
      return      
      end	!ends vel_loadpar
c
