c     
c     jacobian calculation
c     
      subroutine r_jacob(x,xj,xji,det)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension x(3,9),xj(3,3),xji(3,3),cf(3,3)
c     compute the jacobians
      do i=1,3
         do j=1,3
            dum=0.0d0
            do k=1,nis
               dum=dum+r_p(i,k)*x(j,k)
		  enddo
            xj(i,j)=dum
	   enddo
	enddo
c     compute the determinant of the jacobian matrix
c   2-D determinant
c      det=xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
c   3-D determinant
	cf(1,1) = + (xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3))
      cf(1,2) = - (xj(1,2)*xj(3,3) - xj(3,2)*xj(1,3))
      cf(1,3) = + (xj(1,2)*xj(2,3) - xj(2,2)*xj(1,3))
      cf(2,1) = - (xj(2,1)*xj(3,3) - xj(3,1)*xj(2,3))
      cf(2,2) = + (xj(1,1)*xj(3,3) - xj(3,1)*xj(1,3))
      cf(2,3) = - (xj(1,1)*xj(2,3) - xj(2,1)*xj(1,3))
      cf(3,1) = + (xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2))
      cf(3,2) = - (xj(1,1)*xj(3,2) - xj(3,1)*xj(1,2))
      cf(3,3) = + (xj(1,1)*xj(2,2) - xj(2,1)*xj(1,2))

      det = ( xj(1,1) * cf(1,1)
     &     + xj(2,1) * cf(1,2)
     &     + xj(3,1) * cf(1,3) )

c	det=xj(1,1)*(xj(2,2)*xj(3,3)-xj(2,3)*xj(3,2))
c     +   -xj(1,2)*(xj(2,1)*xj(3,3)-xj(3,1)*xj(2,3))
c     +   +xj(1,3)*(xj(2,1)*xj(3,2)-xj(3,1)*xj(2,2))
c      det=xj(1,1)*xj(2,2)*xj(3,3)+
c     $     xj(1,2)*xj(2,3)*xj(3,1)+
c     $     xj(2,1)*xj(3,2)*xj(1,3)-
c     $     xj(1,3)*xj(2,2)*xj(3,1)-
c     $     xj(1,2)*xj(2,1)*xj(3,3)-
c     $     xj(1,1)*xj(2,3)*xj(3,2)
      if (det .lt. 1.0d-15) then
         write(*,100) 
         stop
      endif
  100 format(6x, 'error, zero or negative jacobian determinant')
c     compute the inverse of the jacobian matrix
c   2-D inverse
c      xji(1,1)=xj(2,2)/det
c      xji(1,2)=-xj(1,2)/det
c      xji(2,1)=-xj(2,1)/det
c      xji(2,2)=xj(1,1)/det
c   3-D inverse
c	xji(1,1)=xj(2,2)*xj(3,3)-xj(3,2)*xj(2,3)
c	xji(1,2)=xj(2,3)*xj(3,1)-xj(2,1)*xj(3,3)
c	xji(1,3)=xj(2,1)*xj(3,2)-xj(2,2)*xj(3,1)
c	xji(2,1)=xj(1,3)*xj(3,2)-xj(1,2)*xj(3,3)
c	xji(2,2)=xj(1,1)*xj(3,3)-xj(1,3)*xj(3,1)
c     xji(2,3)=xj(3,1)*xj(1,2)-xj(1,1)*xj(3,2)
c	xji(3,1)=xj(1,2)*xj(2,3)-xj(1,3)*xj(2,2)
c	xji(3,2)=xj(1,3)*xj(2,1)-xj(1,1)*xj(2,3)
c	xji(3,3)=xj(1,1)*xj(2,2)-xj(1,2)*xj(2,1)	
c	xji(:,:)=xji(:,:)/det
      xji(1,1) = cf(1,1)/det
      xji(1,2) = cf(1,2)/det
      xji(1,3) = cf(1,3)/det
      xji(2,1) = cf(2,1)/det
      xji(2,2) = cf(2,2)/det
      xji(2,3) = cf(2,3)/det
      xji(3,1) = cf(3,1)/det
      xji(3,2) = cf(3,2)/det
      xji(3,3) = cf(3,3)/det

	return
      end
