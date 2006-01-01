!     
!     material constant c,cpp,cuu,cup
!     

!cccccccccccccccccccccccccc
!     oCUPP, oCUP, oCPP
!cccccccccccccccccccccccccc

subroutine r_sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
  use r_common, only: rc1,rc2,rk,cpre,bpre,dbpre,ddbpre
  use solid_variables, only: nsd_solid
  
  implicit none 

  real(8) :: ocpp
  real(8) :: obc(2*nsd_solid,2*nsd_solid),xmj(3),dxmj(3,2*nsd_solid)
  real(8) :: ddxmj(3,2*nsd_solid,2*nsd_solid),ocuu(2*nsd_solid,2*nsd_solid),ocup(2*nsd_solid)

  integer :: i,j

  ocpp=-1.0d0/rk
  do i=1,2*nsd_solid
     ocup(i)=-ocpp*dbpre(i)
     do j=1,2*nsd_solid
        obc(i,j)  = rc1*ddxmj(1,i,j) + rc2*ddxmj(2,i,j) + rk*dxmj(3,i)*dxmj(3,j) + rk*ddxmj(3,i,j)*(xmj(3)-1.0d0)
        ocuu(i,j) = obc(i,j) + ocpp*dbpre(i)*dbpre(j) + ocpp*(bpre-cpre)*ddbpre(i,j)
     enddo
  enddo

  return
end subroutine r_sboc
