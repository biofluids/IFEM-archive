c     
c     bar 2nd piola-kichhoff
c     
      subroutine r_spiola(ocpp,xmj,dxmj)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension xmj(3),dxmj(3,6)
      do 10 i=1,4
         btos(i)=rc1*dxmj(1,i)+rc2*dxmj(2,i)+
     $        rk*(xmj(3)-1.0d0)*dxmj(3,i)
         tos(i)=btos(i)+ocpp*(bpre-cpre)*dbpre(i)

   10 continue
      return
      end
