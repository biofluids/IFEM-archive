c     
c     pressure interpolation
c
      subroutine r_spress(rs,ne)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension rs(2)
      r=rs(1)
      s=rs(2)
      hp(1)=1.0d0
      hp(2)=r
      hp(3)=s
      hp(4)=r*s
      hp(5)=r**2
      hp(6)=s**2
      cpre=0.0d0
      do 10 i=1,nump
         ntt=(ne-1)*nump+i
         cpre=cpre+prec(ntt)*hp(i)
 10   continue
      return
      end

