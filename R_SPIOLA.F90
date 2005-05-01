!     
!     bar 2nd piola-kichhoff
!     
subroutine r_spiola(xmj,dxmj,xto)
  use r_common
  implicit none

  real(8) :: xto(3,3)

  real(8) :: PK2str_tens(3,3)

  integer :: i,isd,jsd,ksd

  real(8) :: xmj(3),dxmj(3,6)
  do i=1,6
     bPK2str(i) = rc1*dxmj(1,i) + rc2*dxmj(2,i) + rk*(xmj(3)-1.0d0)*dxmj(3,i)
      PK2str(i) = bPK2str(i) !+ ocpp*(bpre-cpre)*dbpre(i)
  enddo
  
  PK2str_tens(1,1) = PK2str(1)
  PK2str_tens(1,2) = PK2str(6)
  PK2str_tens(1,3) = PK2str(5)
  PK2str_tens(2,1) = PK2str(6)
  PK2str_tens(2,2) = PK2str(2)
  PK2str_tens(2,3) = PK2str(4)
  PK2str_tens(3,1) = PK2str(5)
  PK2str_tens(3,2) = PK2str(4)
  PK2str_tens(3,3) = PK2str(3)


  do isd = 1,3
     do jsd = 1,3
        PK1str_tens(isd,jsd) = 0.0d0
        do ksd = 1,3
           PK1str_tens(isd,jsd) = PK1str_tens(isd,jsd) + PK2str_tens(isd,ksd)*xto(jsd,ksd)
        
		enddo
     enddo
  enddo

  return
end subroutine r_spiola
