c     
c     jacobian,bd,pd calculation
c     
      subroutine r_bdpd(xji)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      dimension xji(2,2)
c     compute the bd matrix
      do 10 i=1,nis
         do 11 k=1,2
            dumcd=0.0d0
            do 12 j=1,2
               dumcd=dumcd+xji(k,j)*r_p(j,i)
   12       continue
            bd(k,i)=dumcd
   11    continue
   10 continue
      return
      end
