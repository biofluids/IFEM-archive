c     
c     jacobian calculation
c     
      subroutine r_jacob(x,xj,xji,det)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension x(2,9),xj(2,2),xji(2,2)
c     compute the jacobians
      do 10 i=1,2
         do 11 j=1,2
            dum=0.0d0
            do 12 k=1,nis
               dum=dum+r_p(i,k)*x(j,k)
   12       continue
            xj(i,j)=dum
   11    continue
   10 continue
c     compute the determinant of the jacobian matrix
      det=xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
      if (det .lt. 1.0d-15) then
         write(*,100) 
         stop
      endif
  100 format(6x, 'error, zero or negative jacobian determinant')
c     compute the inverse of the jacobian matrix
      xji(1,1)=xj(2,2)/det
      xji(1,2)=-xj(1,2)/det
      xji(2,1)=-xj(2,1)/det
      xji(2,2)=xj(1,1)/det
      return
      end
