c     
c     caculation of kup
c     
      subroutine r_scalkup(fkup,ocup,i,k,ni)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension ocup(6)
      fkup=0.0d0
      do m=1,6
         fkup=fkup+ocup(m)*dge(m,i,ni)*hp(k)
	enddo
      return
      end


