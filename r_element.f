c     
c     isoparametric implementation
c     element calculation
c
      subroutine r_element(rs)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension rs(2)
      r=rs(1)
      s=rs(2)
      h(9)=0.0d0
      if (nis .eq. 9) then 
         h(9)=(1.0d0-r**2)*(1.0d0-s**2)
      endif
      if (nis .ne. 4) then
         h(5)=0.5d0*(1.0d0-r**2)*(1.0d0+s)-0.5d0*h(9)
         h(6)=0.5d0*(1.0d0-s**2)*(1.0d0-r)-0.5d0*h(9)
         h(7)=0.5d0*(1.0d0-r**2)*(1.0d0-s)-0.5d0*h(9)
         h(8)=0.5d0*(1.0d0-s**2)*(1.0d0+r)-0.5d0*h(9)
         h(1)=0.25d0*(1.0d0+r)*(1.0d0+s)-0.5d0*h(5)-
     $        0.5d0*h(8)-0.25d0*h(9)
         h(2)=0.25d0*(1.0d0-r)*(1.0d0+s)-0.5d0*h(5)-
     $        0.5d0*h(6)-0.25d0*h(9)
         h(3)=0.25d0*(1.0d0-r)*(1.0d0-s)-0.5d0*h(6)-
     $        0.5d0*h(7)-0.25d0*h(9)
         h(4)=0.25d0*(1.0d0+r)*(1.0d0-s)-0.5d0*h(7)-
     $        0.5d0*h(8)-0.25d0*h(9)
      endif
      if (nis .eq. 4) then 
         h(1)=0.25d0*(1.0d0+r)*(1.0d0+s)
         h(2)=0.25d0*(1.0d0-r)*(1.0d0+s)
         h(3)=0.25d0*(1.0d0-r)*(1.0d0-s)
         h(4)=0.25d0*(1.0d0+r)*(1.0d0-s)
      endif
c     
c     derivative of interpolation function
c     
c     first derivative with respect to r
c     
      r_p(1,9)=0.0d0
      if (nis .eq. 9) then 
         r_p(1,9)=-2.0d0*r*(1.0d0-s**2)      
      endif
      if (nis .ne. 4) then      
c     element node 5
         r_p(1,5)=-r*(1.0d0+s)-0.5d0*r_p(1,9)
c     element node 6
         r_p(1,6)=-0.5d0*(1.0d0-s**2)-0.5d0*r_p(1,9)
c     element node 7
         r_p(1,7)=-r*(1.0d0-s)-0.5d0*r_p(1,9)
c     element node 8 
         r_p(1,8)=0.5d0*(1.0d0-s**2)-0.5d0*r_p(1,9)
c     element node 1
         r_p(1,1)=0.25d0*(1.0d0+s)-0.5d0*r_p(1,5)-
     $        0.5d0*r_p(1,8)-0.25d0*r_p(1,9)
c     element node 2
         r_p(1,2)=-0.25d0*(1.0d0+s)-0.5d0*r_p(1,5)-
     $        0.5d0*r_p(1,6)-0.25d0*r_p(1,9)
c     element node 3
         r_p(1,3)=-0.25d0*(1.0d0-s)-0.5d0*r_p(1,6)-
     $        0.5d0*r_p(1,7)-0.25d0*r_p(1,9)
c     element node 4
         r_p(1,4)=0.25d0*(1.0d0-s)-0.5d0*r_p(1,7)-
     $        0.5d0*r_p(1,8)-0.25d0*r_p(1,9)
      endif
      if (nis .eq. 4) then
         r_p(1,1)=0.25d0*(1.0d0+s)
         r_p(1,2)=-0.25d0*(1.0d0+s)
         r_p(1,3)=-0.25d0*(1.0d0-s)
         r_p(1,4)=0.25d0*(1.0d0-s)
      endif
c     
c     first derivative with respect to s
c     
      r_p(2,9)=0.0d0
      if (nis .eq. 9) then 
         r_p(2,9)=-2.0d0*s*(1.0d0-r**2)         
      endif
      if (nis .ne. 4) then
c     element node 5      
         r_p(2,5)=0.5d0*(1.0d0-r**2)-0.5d0*r_p(2,9)
c     element node 6
         r_p(2,6)=-s*(1.0d0-r)-0.5d0*r_p(2,9)
c     element node 7
         r_p(2,7)=-0.5d0*(1.0d0-r**2)-0.5d0*r_p(2,9)
c     element node 8
         r_p(2,8)=-s*(1.0d0+r)-0.5d0*r_p(2,9)
c     element node 1
         r_p(2,1)=0.25d0*(1.0d0+r)-0.5d0*r_p(2,5)-
     $        0.5d0*r_p(2,8)-0.25d0*r_p(2,9)
c     element node 2
         r_p(2,2)=0.25d0*(1.0d0-r)-0.5d0*r_p(2,5)-
     $        0.5d0*r_p(2,6)-0.25d0*r_p(2,9)
c     element node 3
         r_p(2,3)=-0.25d0*(1.0d0-r)-0.5d0*r_p(2,6)-
     $        0.5d0*r_p(2,7)-0.25d0*r_p(2,9)
c     element node 4
         r_p(2,4)=-0.25d0*(1.0d0+r)-0.5d0*r_p(2,7)-
     $        0.5d0*r_p(2,8)-0.25d0*r_p(2,9)
      endif
      if (nis .eq. 4) then
         r_p(2,1)=0.25d0*(1.0d0+r)
         r_p(2,2)=0.25d0*(1.0d0-r)
         r_p(2,3)=-0.25d0*(1.0d0-r)
         r_p(2,4)=-0.25d0*(1.0d0+r)         
      endif
      return
      end
