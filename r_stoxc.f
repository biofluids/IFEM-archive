c     
c     calculate deformation gradient
c     
      subroutine r_stoxc(xto,xot,xj,xji,toxj,toxji,toc)
      implicit real*8 (a-h,o-z) 
      dimension xto(2,2),xot(2,2),xj(2,2),xji(2,2),toc(3,3),
     $     toxj(2,2),toxji(2,2)
      do 10 i=1,2
         do 11 j=1,2
            xto(i,j)=0.0d0
            xot(i,j)=0.0d0
            do 12 m=1,2
               xto(i,j)=xto(i,j)+xj(m,i)*toxji(j,m)
               xot(i,j)=xot(i,j)+toxj(m,i)*xji(j,m)
   12       continue
   11    continue
   10 continue
      do 13 i=1,2
         do 14 j=1,2
            toc(i,j)=0.0d0
            do 15 m=1,2
               toc(i,j)=toc(i,j)+xto(m,i)*xto(m,j)
   15       continue
   14    continue
   13 continue
      return
      end
