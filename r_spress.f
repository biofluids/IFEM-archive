c     
c     pressure interpolation
c
      subroutine r_spress(rs,ne)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension rs(3)
	   r=rs(1)
         s=rs(2)
         t=rs(3)
         hp(1)=1.0d0
         hp(2)=r
         hp(3)=s
         hp(4)=t
      cpre=0.0d0
      do i=1,nump
         ntt=(ne-1)*nump+i
         cpre=cpre+prec(ntt)*hp(i)
	enddo
      return
      end

