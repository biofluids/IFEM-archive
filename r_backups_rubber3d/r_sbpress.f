c     
c     bar pressure and derivative 
c
      subroutine r_sbpress(dxmj,ddxmj,xmj)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension xmj(3),dxmj(3,6),ddxmj(3,6,6)
ccccccccccccccccccccccccccccccccccccccc
      bpre=-rk*(xmj(3)-1.0d0)
      do i=1,6
         dbpre(i)=-rk*dxmj(3,i)
         do j=1,6
            ddbpre(i,j)=-rk*ddxmj(3,i,j)
	   enddo
	enddo
      return
      end
