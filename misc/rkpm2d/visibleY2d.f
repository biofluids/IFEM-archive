c
c****************************************************************
c     subroutine visibleY2d(xbar,n_pt,Lnsp)
c
c     Purpose:
c     =======
c     Using visibility criterion to find the effective nodal
c     points inside the support.
c     
c     The sampling point position: xbar(2): it could be Gauss pt
c                                           or nodal pt. as well
c
c
c     The crack surface position: Cyplane; xCmin,xCmax
c
c     n_pt: total particles in the domain;
c     Lnsp: local connectivity map;
c
c     local parameters:
c
c     pt(2,mnschG): interception point of crack plane and line
c                   segment;
c
c
c     intersect(i) = true;  intersected
c                    flase; not intersected
c
c     testngbr(i) =  0 : exluded
c                    1 : included
c
c     Date: July, 1999
c
c*****************************************************************
c
      subroutine visibleY2d(xbar,n_pt,Lnsp)
      implicit none
      include 'parameter.h'
c
      integer n_pt,Lnsp(mnsch)
      real*8  xbar(2),eps
c
c......Local variables
c
      logical intersect(mnsch)
      integer testngbr(mnsch),n_shift(mnsch)
      integer i,ipt,ncount,ncount1
c
c     declarations for the nodal search 
c     (neighbors for the sampling point)
c
      real*8 dx(mnsch),dy(mnsch),
     &       pt(2,mnsch),dyp,
     &       dypI(mnsch),t(mnsch)
c
      real*8 xm(2,maxNumnp),dxm(2,maxNumnp),
     &       dvm(maxNumnp)
      real*8  Cyplane,Cxmin,Cxmax
      integer iCrack,numnp
c
      common /mesh/xm,dxm,dvm,numnp
      common /crack/Cyplane,Cxmin,Cxmax
      common /Type1/iCrack
c 
c
          eps = 1.0d-10
          do i = 1, n_pt
	     ipt   = Lnsp(i)
             dx(i) = xbar(1) - xm(1,ipt)
             dy(i) = xbar(2) - xm(2,ipt)
c
             testngbr(i) = 1
             n_shift(i)  = 0
             t(i)        = 0.0
             pt(1,i)     = 0.0
             pt(2,i)     = 0.0
          enddo ! n_pt
c
c......Check nodes for the visibility criterion............
c
         dyp  = Cyplane - xbar(2)    ! vertical distance s-->c
         if (abs(dyp) .le. eps) then !P-bar is parallel with crack
	    return
	 endif
c
         do i = 1,n_pt
	    ipt     = Lnsp(i)
            dypI(i) = Cyplane - xm(2,ipt)  ! vertical distance i-->c
c
c     determine if the line joining the node to 
c     the sampling point intersects the crack plane
c
	    if (abs(dypI(i)) .le. eps ) then ! nodal is parallel with crack
	       intersect(i) = .false.
            elseif (sign(1.0,dyp) .eq. sign(1.0,dypI(i))) then
               intersect(i) = .false.
            else
               intersect(i) = .true.
               t(i)    = dypI(i)/dy(i)
               pt(1,i) = dx(i)*t(i) + xm(1,ipt)
               pt(2,i) = Cyplane
            endif
          enddo !i
c
	 ncount = 0
         do i = 1,n_pt
            if (intersect(i)) then
c
               if ((dabs(pt(1,i)) .ge. Cxmin) .and.
     &             (dabs(pt(1,i)) .le. Cxmax)) then 
                     testngbr(i) = 0   ! the point should be excluded
		     ncount      = ncount + 1 
               else
		  ncount = ncount
               endif
c                                  
            endif ! intersect = .true.
	    n_shift(i) = ncount
         enddo  ! n_pt
c
c........Modify the connectivity map..............
c
	 ncount1 = 0
	 do i = 1, n_pt
	    if (testngbr(i) .eq. 1) then ! shift the useful pt.  
               ncount1         = n_shift(i)
               Lnsp(i-ncount1) = Lnsp(i)
	    endif
	 enddo
c
c....Update n_pt
c
	n_pt = n_pt - ncount 
c
      return
      end
c
