c     
c     material constant c,cpp,cuu,cup
c     
      subroutine r_sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension obc(6,6),xmj(3),dxmj(3,6),ddxmj(3,6,6),
     $     ocuu(6,6),ocup(6)
cccccccccccccccccccccccccc
c     oCUPP, oCUP, oCPP
cccccccccccccccccccccccccc
      ocpp=-1.0d0/rk
      do 10 i=1,3
         ocup(i)=-ocpp*dbpre(i)
         do 11 j=1,3
            obc(i,j)=rc1*ddxmj(1,i,j)+rc2*ddxmj(2,i,j)+
     $           rk*dxmj(3,i)*dxmj(3,j)+
     $           rk*ddxmj(3,i,j)*(xmj(3)-1.0d0)
            ocuu(i,j)=obc(i,j)+ocpp*dbpre(i)*dbpre(j)+
     $           ocpp*(bpre-cpre)*ddbpre(i,j)
   11    continue
   10 continue
      return
      end
