!     
!     material constant c,cpp,cuu,cup
!     

!cccccccccccccccccccccccccc
!     oCUPP, oCUP, oCPP
!cccccccccccccccccccccccccc

subroutine r_sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
  use r_common, only: rc1,rc2,rk,cpre,bpre,dbpre,ddbpre
  implicit none 

  real(8) :: ocpp
  real(8) :: obc(6,6),xmj(3),dxmj(3,6),ddxmj(3,6,6),ocuu(6,6),ocup(6)

  integer :: i,j

  ocpp=-1.0d0/rk
  do i=1,6
     ocup(i)=-ocpp*dbpre(i)
     do j=1,6
        obc(i,j)  = rc1*ddxmj(1,i,j) + rc2*ddxmj(2,i,j) + rk*dxmj(3,i)*dxmj(3,j) + rk*ddxmj(3,i,j)*(xmj(3)-1.0d0)
        ocuu(i,j) = obc(i,j) + ocpp*dbpre(i)*dbpre(j) + ocpp*(bpre-cpre)*ddbpre(i,j)
     enddo
  enddo

  return
end subroutine r_sboc
