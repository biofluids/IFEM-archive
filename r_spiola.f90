!     
!     bar 2nd piola-kichhoff
!     
subroutine r_spiola(ocpp,xmj,dxmj,xto)
  use r_common
  implicit none

  real*8 :: ocpp
  real*8 :: xto(3,3)

  integer :: i

  real*8 :: xmj(3),dxmj(3,6)
  do i=1,6
     bPK2str(i) = rc1*dxmj(1,i) + rc2*dxmj(2,i) + rk*(xmj(3)-1.0d0)*dxmj(3,i)
      PK2str(i) = bPK2str(i) !+ ocpp*(bpre-cpre)*dbpre(i)
  enddo

  PK1str(1) = PK2str(1)*xto(1,1) + PK2str(6)*xto(1,2) + PK2str(5)*xto(1,3)
  PK1str(2) = PK2str(6)*xto(2,1) + PK2str(2)*xto(2,2) + PK2str(4)*xto(2,3)
  PK1str(3) = PK2str(5)*xto(3,1) + PK2str(4)*xto(3,2) + PK2str(3)*xto(3,3)
  PK1str(4) = PK2str(6)*xto(3,1) + PK2str(2)*xto(3,2) + PK2str(4)*xto(3,3)
  PK1str(5) = PK2str(1)*xto(3,1) + PK2str(6)*xto(3,2) + PK2str(5)*xto(3,3)
  PK1str(6) = PK2str(1)*xto(2,1) + PK2str(6)*xto(2,2) + PK2str(5)*xto(2,3)


  return
end subroutine r_spiola
