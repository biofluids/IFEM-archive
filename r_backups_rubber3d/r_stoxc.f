c     
c     calculate deformation gradient
c     
      subroutine r_stoxc(xto,xot,xj,xji,toxj,toxji,toc)
      implicit real*8 (a-h,o-z) 
      dimension xto(3,3),xot(3,3),xj(3,3),xji(3,3),toc(3,3),
     $     toxj(3,3),toxji(3,3)
      do i=1,3
         do j=1,3
            xto(i,j)=0.0d0
            xot(i,j)=0.0d0
            do m=1,3
               xto(i,j)=xto(i,j)+xj(m,i)*toxji(j,m)
               xot(i,j)=xot(i,j)+toxj(m,i)*xji(j,m)
		  enddo
	   enddo
	enddo

      do i=1,3
         do j=1,3
            toc(i,j)=0.0d0
            do m=1,3
               toc(i,j)=toc(i,j)+xto(m,i)*xto(m,j)
		  enddo
	   enddo
	enddo

      return
      end