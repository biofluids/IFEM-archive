c     
c     jacobian,bd,pd calculation
c     
      subroutine r_bdpd(xji)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      dimension xji(3,3)
c     compute the bd matrix-> strain displacement matrix
      do i=1,nis
         do k=1,3
            dumcd=0.0d0
            do j=1,3
               dumcd=dumcd+xji(j,k)*r_p(j,i)
		  enddo
            bd(k,i)=dumcd
	   enddo
	enddo

      return
      end
