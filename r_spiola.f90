!     
!     bar 2nd piola-kichhoff
!     
subroutine r_spiola(ocpp,xmj,dxmj)
  use r_common
  implicit none

  real*8 :: ocpp

  integer :: i

  real*8 :: xmj(3),dxmj(3,6)
  do i=1,6
     btos(i) = rc1*dxmj(1,i) + rc2*dxmj(2,i) + rk*(xmj(3)-1.0d0)*dxmj(3,i)
      tos(i) = btos(i) + ocpp*(bpre-cpre)*dbpre(i)
  enddo

  return
end subroutine r_spiola
