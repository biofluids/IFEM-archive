!     
!     jacobian calculation
!     
subroutine r_jacob(x,xj,xji,det)
  use solid_variables, only:nsd_solid,nen_solid
  use r_common
  implicit none

  real(8) :: x(nsd_solid,nen_solid),xj(nsd_solid,nsd_solid),xji(nsd_solid,nsd_solid),cf(nsd_solid,nsd_solid)
  real(8),intent(out) :: det


  real(8) :: dum
  integer :: i,j,k

 !...compute the jacobians
  do i=1,nsd_solid
     do j=1,nsd_solid
        dum = 0.0d0
        do k=1,nen_solid
           dum = dum + r_p(j,k)*x(i,k)
        enddo
        xj(i,j) = dum
     enddo
  enddo
 !...compute the determinant of the jacobian matrix

  threedim: if (nsd_solid .eq. 3) then 
 !...3-D determinant
  cf(1,1) = + (xj(2,2)*xj(3,3) - xj(2,3)*xj(3,2))
  cf(1,2) = - (xj(2,1)*xj(3,3) - xj(2,3)*xj(3,1))
  cf(1,3) = + (xj(2,1)*xj(3,2) - xj(2,2)*xj(3,1))
  cf(2,1) = - (xj(1,2)*xj(3,3) - xj(1,3)*xj(3,2))
  cf(2,2) = + (xj(1,1)*xj(3,3) - xj(1,3)*xj(3,1))
  cf(2,3) = - (xj(1,1)*xj(3,2) - xj(1,2)*xj(3,1))
  cf(3,1) = + (xj(1,2)*xj(2,3) - xj(1,3)*xj(2,2))
  cf(3,2) = - (xj(1,1)*xj(2,3) - xj(1,3)*xj(2,1))
  cf(3,3) = + (xj(1,1)*xj(2,2) - xj(1,2)*xj(2,1))

  det =( xj(1,1) * cf(1,1) + &
         xj(1,2) * cf(1,2) + &
         xj(1,3) * cf(1,3) )

  if (det .lt. 1.0d-15) then
     write(*,100) 
     stop
  endif
100 format(6x, 'error, zero or negative jacobian determinant')


!...compute the inverse of the jacobian matrix

  xji(1,1) = cf(1,1)/det
  xji(1,2) = cf(2,1)/det
  xji(1,3) = cf(3,1)/det
  xji(2,1) = cf(1,2)/det
  xji(2,2) = cf(2,2)/det
  xji(2,3) = cf(3,2)/det
  xji(3,1) = cf(1,3)/det
  xji(3,2) = cf(2,3)/det
  xji(3,3) = cf(3,3)/det

  endif  threedim

  twodim: if (nsd_solid .eq. 2) then 
 !...2-D determinant
  det = ( xj(1,1) * xj(2,2) - &
         xj(1,2) * xj(2,1) )

  if (det .lt. 1.0d-15) then
     write(*,100) 
     stop
  endif
	
!...compute the inverse of the jacobian matrix

  xji(1,1) = xj(2,2)/det
  xji(1,2) = -xj(1,2)/det
  xji(2,1) = -xj(2,1)/det
  xji(2,2) = xj(1,1)/det

	endif twodim


  return
end subroutine r_jacob
