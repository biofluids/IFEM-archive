c     
c     isoparametric implementation
c     element calculation
c
      subroutine r_element(rs)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension rs(3)
      r=rs(1)
      s=rs(2)
      t=rs(3)
*
      if (nis .eq. 8) then 
         h(1)=0.125d0*(1.0d0+r)*(1.0d0+s)*(1.0d0+t)
         h(2)=0.125d0*(1.0d0-r)*(1.0d0+s)*(1.0d0+t)
         h(3)=0.125d0*(1.0d0-r)*(1.0d0-s)*(1.0d0+t)
         h(4)=0.125d0*(1.0d0+r)*(1.0d0-s)*(1.0d0+t)
         h(5)=0.125d0*(1.0d0+r)*(1.0d0+s)*(1.0d0-t)
         h(6)=0.125d0*(1.0d0-r)*(1.0d0+s)*(1.0d0-t)
         h(7)=0.125d0*(1.0d0-r)*(1.0d0-s)*(1.0d0-t)
         h(8)=0.125d0*(1.0d0+r)*(1.0d0-s)*(1.0d0-t)
      endif
c     
c     derivative of interpolation function
c     
c     first derivative with respect to r
c     
      if (nis .eq. 8) then
         r_p(1,1)= 0.125d0*(1.0d0+s)*(1.0d0+t)
         r_p(1,2)=-0.125d0*(1.0d0+s)*(1.0d0+t)
         r_p(1,3)=-0.125d0*(1.0d0-s)*(1.0d0+t)
         r_p(1,4)= 0.125d0*(1.0d0-s)*(1.0d0+t)
         r_p(1,5)= 0.125d0*(1.0d0+s)*(1.0d0-t)
         r_p(1,6)=-0.125d0*(1.0d0+s)*(1.0d0-t)
         r_p(1,7)=-0.125d0*(1.0d0-s)*(1.0d0-t)
         r_p(1,8)= 0.125d0*(1.0d0-s)*(1.0d0-t)
      endif
c     
c     first derivative with resr_pect to s
c     
      if (nis .eq. 8) then
         r_p(2,1)= 0.125d0*(1.0d0+r)*(1.0d0+t)
         r_p(2,2)= 0.125d0*(1.0d0-r)*(1.0d0+t)
         r_p(2,3)=-0.125d0*(1.0d0-r)*(1.0d0+t)
         r_p(2,4)=-0.125d0*(1.0d0+r)*(1.0d0+t)         
         r_p(2,5)= 0.125d0*(1.0d0+r)*(1.0d0-t)
         r_p(2,6)= 0.125d0*(1.0d0-r)*(1.0d0-t)
         r_p(2,7)=-0.125d0*(1.0d0-r)*(1.0d0-t)
         r_p(2,8)=-0.125d0*(1.0d0+r)*(1.0d0-t)      
	endif   
c     
c     first derivative with resr_pect to t
c     
      if (nis .eq. 8) then
         r_p(3,1)= 0.125d0*(1.0d0+r)*(1.0d0+s)
         r_p(3,2)= 0.125d0*(1.0d0-r)*(1.0d0+s)
         r_p(3,3)= 0.125d0*(1.0d0-r)*(1.0d0-s)
         r_p(3,4)= 0.125d0*(1.0d0+r)*(1.0d0-s)         
         r_p(3,5)=-0.125d0*(1.0d0+r)*(1.0d0+s)
         r_p(3,6)=-0.125d0*(1.0d0-r)*(1.0d0+s)
         r_p(3,7)=-0.125d0*(1.0d0-r)*(1.0d0-s)
         r_p(3,8)=-0.125d0*(1.0d0+r)*(1.0d0-s)         
      endif
      return
      end