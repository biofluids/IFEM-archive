!     
!     isoparametric implementation
!     element calculation
!
subroutine r_element(rs)
  use solid_variables, only: nsd_solid,nen_solid
  use r_common, only: h,r_p
  implicit none

  real(8) :: rs(1:nsd_solid)
  real(8) :: r,s,t

  cube: if (nen_solid .eq. 8) then 
		
	 r=rs(1)
	 s=rs(2)
	 t=rs(3)

     h(1)=0.125d0*(1.0d0+r)*(1.0d0+s)*(1.0d0+t)
     h(2)=0.125d0*(1.0d0-r)*(1.0d0+s)*(1.0d0+t)
     h(3)=0.125d0*(1.0d0-r)*(1.0d0-s)*(1.0d0+t)
     h(4)=0.125d0*(1.0d0+r)*(1.0d0-s)*(1.0d0+t)
     h(5)=0.125d0*(1.0d0+r)*(1.0d0+s)*(1.0d0-t)
     h(6)=0.125d0*(1.0d0-r)*(1.0d0+s)*(1.0d0-t)
     h(7)=0.125d0*(1.0d0-r)*(1.0d0-s)*(1.0d0-t)
     h(8)=0.125d0*(1.0d0+r)*(1.0d0-s)*(1.0d0-t)

!     derivative of interpolation function
!     first derivative with respect to r
     r_p(1,1)= 0.125d0*(1.0d0+s)*(1.0d0+t)
     r_p(1,2)=-0.125d0*(1.0d0+s)*(1.0d0+t)
     r_p(1,3)=-0.125d0*(1.0d0-s)*(1.0d0+t)
     r_p(1,4)= 0.125d0*(1.0d0-s)*(1.0d0+t)
     r_p(1,5)= 0.125d0*(1.0d0+s)*(1.0d0-t)
     r_p(1,6)=-0.125d0*(1.0d0+s)*(1.0d0-t)
     r_p(1,7)=-0.125d0*(1.0d0-s)*(1.0d0-t)
     r_p(1,8)= 0.125d0*(1.0d0-s)*(1.0d0-t)

!     first derivative with resr_pect to s
     r_p(2,1)= 0.125d0*(1.0d0+r)*(1.0d0+t)
     r_p(2,2)= 0.125d0*(1.0d0-r)*(1.0d0+t)
     r_p(2,3)=-0.125d0*(1.0d0-r)*(1.0d0+t)
     r_p(2,4)=-0.125d0*(1.0d0+r)*(1.0d0+t)         
     r_p(2,5)= 0.125d0*(1.0d0+r)*(1.0d0-t)
     r_p(2,6)= 0.125d0*(1.0d0-r)*(1.0d0-t)
     r_p(2,7)=-0.125d0*(1.0d0-r)*(1.0d0-t)
     r_p(2,8)=-0.125d0*(1.0d0+r)*(1.0d0-t)      
  
!     first derivative with resr_pect to t
     r_p(3,1)= 0.125d0*(1.0d0+r)*(1.0d0+s)
     r_p(3,2)= 0.125d0*(1.0d0-r)*(1.0d0+s)
     r_p(3,3)= 0.125d0*(1.0d0-r)*(1.0d0-s)
     r_p(3,4)= 0.125d0*(1.0d0+r)*(1.0d0-s)         
     r_p(3,5)=-0.125d0*(1.0d0+r)*(1.0d0+s)
     r_p(3,6)=-0.125d0*(1.0d0-r)*(1.0d0+s)
     r_p(3,7)=-0.125d0*(1.0d0-r)*(1.0d0-s)
     r_p(3,8)=-0.125d0*(1.0d0+r)*(1.0d0-s)         
  endif cube

  tetr: if ((nen_solid .eq. 4).AND.(nsd_solid.eq.3)) then 

  	 r=rs(1)
	 s=rs(2)
	 t=rs(3)

     h(1)=r
     h(2)=s
     h(3)=t
     h(4)=1-r-s-t
!     
!     derivative of interpolation function
!     first derivative with respect to r
     r_p(1,1)= 1.0d0
     r_p(1,2)= 0.0d0
     r_p(1,3)= 0.0d0
     r_p(1,4)=-1.0d0
! 
!     first derivative with respect to s
     r_p(2,1)= 0.0d0
     r_p(2,2)= 1.0d0
     r_p(2,3)= 0.0d0
     r_p(2,4)=-1.0d0         
!  
!     first derivative with respect to t
     r_p(3,1)= 0.0d0
     r_p(3,2)= 0.0d0
     r_p(3,3)= 1.0d0
     r_p(3,4)=-1.0d0         
  endif tetr

!	Yaling Liu added triaguler element
!	h(nen_solid), r_p(nsd_solid,nen_solid)
  tria: if (nen_solid .eq. 3) then 

  	 r=rs(1)
	 s=rs(2)

     h(1)=r
     h(2)=s
     h(3)=1-r-s
!     
!     derivative of interpolation function
!     first derivative with respect to r
     r_p(1,1)= 1.0d0
     r_p(1,2)= 0.0d0
     r_p(1,3)=-1.0d0
! 
!     first derivative with respect to s
     r_p(2,1)= 0.0d0
     r_p(2,2)= 1.0d0
     r_p(2,3)=-1.0d0         
	        
  endif tria

  quad: if ((nen_solid .eq. 4).and.(nsd_solid.eq.2)) then 

  	 r=rs(1)
	 s=rs(2)

     h(1)=0.25d0*(1.0d0-r)*(1.0d0-s)
     h(2)=0.25d0*(1.0d0+r)*(1.0d0-s)
     h(3)=0.25d0*(1.0d0+r)*(1.0d0+s)
     h(4)=0.25d0*(1.0d0-r)*(1.0d0+s)

!     derivative of interpolation function
!     first derivative with respect to r
     r_p(1,1)=-0.25d0*(1.0d0+s)
     r_p(1,2)= 0.25d0*(1.0d0+s)
     r_p(1,3)= 0.25d0*(1.0d0-s)
     r_p(1,4)=-0.25d0*(1.0d0-s)
!     
!     first derivative with resr_pect to s
     r_p(2,1)=-0.25d0*(1.0d0+r)
     r_p(2,2)=-0.25d0*(1.0d0-r)
     r_p(2,3)=+0.25d0*(1.0d0-r)
     r_p(2,4)=+0.25d0*(1.0d0+r)        
	        
  endif quad

  return
end subroutine r_element
