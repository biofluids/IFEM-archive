c     
c     caculation of fu
c     
      subroutine r_scalfu(fu,i,ni)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      fu=0.0d0
      do 10 m=1,3
         fu=fu+tos(m)*dge(m,i,ni)
   10 continue
      return
      end
