      subroutine gauss2D(gp_loc2D,gp_weight2D,nintElem,nnode)
c** this subroutine to get the loc coordinates and weights for the 2D elements
c** input parameter:
c   nintElem: the total number of integration points per element
c   nnode:    the node number of each element
c             nnode = 3: triangle element
c                     4: quadrilateral element

      implicit none
      
      include 'parameter.h'
      
      integer nintElem,nnode
      real*8 gp_loc2D(2,maxIntElem),gp_weight2D(maxIntElem)
      
        !** loc vars
      real*8  gp_loc1D(maxInt),gp_weight1D(maxInt)
      integer nint1D,iLocInt
      integer ig,jg
      real*8  xsi,eta

      if (nnode.eq.4) then
         nint1D=sqrt(real(nintElem))
         call gauss1D(gp_loc1D,gp_weight1D,nint1D)
         
         iLocInt=0
         do ig=1,nint1D
            xsi=gp_loc1D(ig)
            do jg=1,nint1D
               iLocInt=iLocInt+1
               
               eta=gp_loc1D(jg)
               
               gp_loc2D(1,iLocInt)=xsi
               gp_loc2D(2,iLocInt)=eta
               gp_weight2D(iLocInt)=gp_weight1D(ig)*gp_weight1D(jg)
            enddo
         enddo
      elseif (nnode.eq.3) then
         if (nintElem.eq.1) then
            gp_loc2D(1,1)=1./3.
            gp_loc2D(2,1)=1./3.
            gp_weight2D(1)=1.
         elseif (nintElem.eq.3) then
            gp_loc2D(1,1)=1./2.
            gp_loc2D(2,1)=1./2.
            gp_weight2D(1)=1./3.  

            gp_loc2D(1,2)=0.
            gp_loc2D(2,2)=1./2.
            gp_weight2D(2)=1./3.  

            gp_loc2D(1,3)=1./2.
            gp_loc2D(2,3)=0.
            gp_weight2D(3)=1./3.
         elseif (nintElem.eq.4) then
            gp_loc2D(1,1)=1./3.
            gp_loc2D(2,1)=1./3.
            gp_weight2D(1)=-27./48.  

            gp_loc2D(1,2)=0.6
            gp_loc2D(2,2)=0.2
            gp_weight2D(2)=25./48.

            gp_loc2D(1,3)=0.2
            gp_loc2D(2,3)=0.6
            gp_weight2D(3)=25./48.
            
            gp_loc2D(1,4)=0.2
            gp_loc2D(2,4)=0.2
            gp_weight2D(4)=25./48.

            
         else
            write(*,*) 'In sub gauss2D,for nnode=3,',
     &                 'unknown nintElem=',nintElem
            stop
         endif
      else
         write(*,*) 'In sub gauss2D, unknown nnode=',nnode
         stop
      endif
      
      return
      end ! ends gauss2D
      
      
      subroutine gauss1D(s,w,nint)
c
c     subroutine to give gaussian pts (up to 10) of 1D
c     for intergration over -1 to 1 !!!!!!
c
      implicit double precision (a-h,o-z)
c
      dimension s(nint),w(nint)

      if ( nint.gt.10) then
         write(*,*) 'nint > 10 in subroutine of gaussian. STOP'
         write(*,*) 'nint=',nint
         stop
      endif

      go to (19,29,39,49,59,69,79,89,99,109) nint
      
      		! ==1
 19   s(1)=0.
      w(1)=2.
      return
      
      		! ==2
 29   s(1)=-0.577350269189626d0
      s(2)=-s(1)
      w(1)=1.d0
      w(2)=1.d0
      return
      
      		! ==3
 39   s(1)=-0.774596669241483d0
      s(2)=0.
      s(3)=-s(1)
      w(1)=0.555555555555556d0
      w(2)=0.888888888888889d0
      w(3)=w(1)
      return
      
      		! ==4
 49   s(1)=-0.861136311594053d0
      s(2)=-0.339981043584856d0
      s(3)=-s(2)
      s(4)=-s(1)
      w(1)=0.347854845137454d0
      w(2)=0.652145154862546d0
      w(3)=w(2)
      w(4)=w(1)
      return
      
      		! ==5
 59   s(1)=-0.906179845938664d0
      s(2)=-0.538469310105683d0
      s(3)=0.
      s(4)=-s(2)
      s(5)=-s(1)
      
      w(1)=0.236936885056189d0
      w(2)=0.478638670499366d0
      w(3)=0.568888888888889d0
      w(4)=w(2)
      w(5)=w(1)
      return
      
      		! ==6
 69   s(1)=-.932469514203252d0
      s(2)=-.661209386466265d0
      s(3)=-.238619186083197d0
      s(4)=-s(3)
      s(5)=-s(2)
      s(6)=-s(1)
      w(1)=0.171324492379170d0
      w(2)=0.360761573048139d0
      w(3)=0.467913934572691d0
      w(4)=w(3)
      w(5)=w(2)
      w(6)=w(1)
      return
      
      		! ==7
 79   s(1)=-0.949107912342759d0
      s(2)=-0.741531185599394d0
      s(3)=-0.405845151377397d0
      s(4)=0.
      s(5)=-s(3)
      s(6)=-s(2)
      s(7)=-s(1)
      w(1)=0.129484966168870d0
      w(2)=0.279705391489277d0
      w(3)=0.381830050505119d0
      w(4)=0.417959183673469d0
      w(5)=w(3)
      w(6)=w(2)
      w(7)=w(1)
      return
      
      		! ==8
 89   s(1)=-0.960289856497536d0
      s(2)=-0.796666477413627d0
      s(3)=-0.525532409916329d0
      s(4)=-0.183434642495650d0
      s(5)=-s(4)
      s(6)=-s(3)
      s(7)=-s(2)
      s(8)=-s(1)                  
      w(1)=0.101228536290376d0
      w(2)=0.222381034453374d0
      w(3)=0.313706645877887d0
      w(4)=0.362683783378362d0
      w(5)=w(4)
      w(6)=w(3)
      w(7)=w(2)
      w(8)=w(1)                  
      return
      
      		! ==9
 99   s(1)=-0.968160239507626d0
      s(2)=-0.836031107326636d0
      s(4)=-0.613371432700590d0
      s(4)=-0.324253423403809d0
      s(5)=0.
      s(6)=-s(4)
      s(7)=-s(3)
      s(8)=-s(2)
      s(9)=-s(1)                  
      w(1)=0.081274388361574d0
      w(2)=0.180648160694857d0
      w(3)=0.260610696402935d0
      w(4)=0.312347077040003d0
      w(5)=0.330239355001260d0
      w(6)=w(4)
      w(7)=w(3)
      w(8)=w(2)
      w(9)=w(1)                  
      return
      
 109  s(1)=-0.973906528517172d0
      s(2)=-0.865063366688985d0
      s(3)=-0.679409568299024d0
      s(4)=-0.433395394129247d0
      s(5)=-0.148874338981631d0
      s(6)=-s(5)
      s(7)=-s(4)
      s(8)=-s(3)
      s(9)=-s(2)
      s(10)=-s(1)
      w(1)=0.066671344308688d0
      w(2)=0.149451349150581d0
      w(3)=0.219086362515982d0
      w(4)=0.269266719309996d0
      w(5)=0.295524224714753d0
      w(6)=w(5)
      w(7)=w(4)
      w(8)=w(3)
      w(9)=w(2)
      w(10)=w(1)                  
      
      return
      end		!ends gauss1D
