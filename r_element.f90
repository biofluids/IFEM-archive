!     
!     isoparametric implementation
!     element calculation
!
subroutine r_element(rs)
  use solid_variables, only: nis
  use r_common, only: h,r_p
  implicit none

  real*8 :: rs(3)

  real*8 :: r,s,t

  r=rs(1)
  s=rs(2)
  t=rs(3)

  cube: if (nis .eq. 8) then 
     h(1)=0.125d0*(1.0d0+r)*(1.0d0+s)*(1.0d0+t)
     h(2)=0.125d0*(1.0d0-r)*(1.0d0+s)*(1.0d0+t)
     h(3)=0.125d0*(1.0d0-r)*(1.0d0-s)*(1.0d0+t)
     h(4)=0.125d0*(1.0d0+r)*(1.0d0-s)*(1.0d0+t)
     h(5)=0.125d0*(1.0d0+r)*(1.0d0+s)*(1.0d0-t)
     h(6)=0.125d0*(1.0d0-r)*(1.0d0+s)*(1.0d0-t)
     h(7)=0.125d0*(1.0d0-r)*(1.0d0-s)*(1.0d0-t)
     h(8)=0.125d0*(1.0d0+r)*(1.0d0-s)*(1.0d0-t)

!     
!     derivative of interpolation function
!     
!     first derivative with respect to r
!     

     r_p(1,1)= 0.125d0*(1.0d0+s)*(1.0d0+t)
     r_p(1,2)=-0.125d0*(1.0d0+s)*(1.0d0+t)
     r_p(1,3)=-0.125d0*(1.0d0-s)*(1.0d0+t)
     r_p(1,4)= 0.125d0*(1.0d0-s)*(1.0d0+t)
     r_p(1,5)= 0.125d0*(1.0d0+s)*(1.0d0-t)
     r_p(1,6)=-0.125d0*(1.0d0+s)*(1.0d0-t)
     r_p(1,7)=-0.125d0*(1.0d0-s)*(1.0d0-t)
     r_p(1,8)= 0.125d0*(1.0d0-s)*(1.0d0-t)

!     
!     first derivative with resr_pect to s
!     

     r_p(2,1)= 0.125d0*(1.0d0+r)*(1.0d0+t)
     r_p(2,2)= 0.125d0*(1.0d0-r)*(1.0d0+t)
     r_p(2,3)=-0.125d0*(1.0d0-r)*(1.0d0+t)
     r_p(2,4)=-0.125d0*(1.0d0+r)*(1.0d0+t)         
     r_p(2,5)= 0.125d0*(1.0d0+r)*(1.0d0-t)
     r_p(2,6)= 0.125d0*(1.0d0-r)*(1.0d0-t)
     r_p(2,7)=-0.125d0*(1.0d0-r)*(1.0d0-t)
     r_p(2,8)=-0.125d0*(1.0d0+r)*(1.0d0-t)      
  
!     
!     first derivative with resr_pect to t
!     

     r_p(3,1)= 0.125d0*(1.0d0+r)*(1.0d0+s)
     r_p(3,2)= 0.125d0*(1.0d0-r)*(1.0d0+s)
     r_p(3,3)= 0.125d0*(1.0d0-r)*(1.0d0-s)
     r_p(3,4)= 0.125d0*(1.0d0+r)*(1.0d0-s)         
     r_p(3,5)=-0.125d0*(1.0d0+r)*(1.0d0+s)
     r_p(3,6)=-0.125d0*(1.0d0-r)*(1.0d0+s)
     r_p(3,7)=-0.125d0*(1.0d0-r)*(1.0d0-s)
     r_p(3,8)=-0.125d0*(1.0d0+r)*(1.0d0-s)         
  endif cube

  tetr: if (nis .eq. 4) then 
     h(1)=r
     h(2)=s
     h(3)=t
     h(4)=1-r-s-t
!     
!     derivative of interpolation function
!     
!     first derivative with respect to r
!     
     r_p(1,1)= 1.0d0
     r_p(1,2)= 0.0d0
     r_p(1,3)= 0.0d0
     r_p(1,4)=-1.0d0
! 
!     first derivative with respect to s
!     
     r_p(2,1)= 0.0d0
     r_p(2,2)= 1.0d0
     r_p(2,3)= 0.0d0
     r_p(2,4)=-1.0d0         
!  
!     first derivative with respect to t
!     
     r_p(3,1)= 0.0d0
     r_p(3,2)= 0.0d0
     r_p(3,3)= 1.0d0
     r_p(3,4)=-1.0d0         
  endif tetr

  return
end subroutine r_element
