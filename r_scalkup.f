c     
c     caculation of kup
c     
      subroutine r_scalkup(fkup,ocup,i,k,ni)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension ocup(6)
      fkup=0.0d0
      do 10 m=1,3
         fkup=fkup+ocup(m)*dge(m,i,ni)*hp(k)
   10 continue
      return
      end


