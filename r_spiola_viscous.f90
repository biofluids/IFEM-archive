!     
!     bar 2nd piola-kichhoff
!     
subroutine r_spiola_viscous(xot,vel)
  use r_common
  use fluid_variables, only: mu => vis_liq
  use solid_variables, only: nen_solid,nsd => nsd_solid
  implicit none

  !real(8) :: ocpp
  !real(8) :: xmj(3),dxmj(3,6)
  real(8),intent(in) :: xot(1:nsd,1:nsd)
  real(8),intent(in) :: vel(1:nsd,1:nen_solid)

  integer :: nos,isd,jsd,ksd
  real(8) :: sigma(1:nsd,1:nsd),ui_j(1:nsd,1:nsd) !,press


  sigma(:,:) = 0.0
  ui_j(:,:) = 0.0
  !press = 0.0
  !mu = vis_liq
  do nos = 1,nen_solid
     do isd = 1,nsd
        do jsd = 1,nsd
           ui_j(isd,jsd) = ui_j(isd,jsd) + bd_curr(jsd,nos)*vel(isd,nos)
        end do
     end do
     !press = press + bd_curr(0,nos)*d(pdf,ien(nos,ie))
  end do
  do isd = 1,nsd
     do jsd = 1,nsd
        sigma(isd,jsd) = mu*(ui_j(isd,jsd)+ui_j(jsd,isd))
     end do
     !sigma(isd,isd) = sigma(isd,isd) - press
  end do

  do isd = 1,nsd
     do jsd = 1,nsd
        do ksd = 1,nsd
           PK1str_tens(isd,jsd) = PK1str_tens(isd,jsd) - xot(isd,ksd)*sigma(ksd,jsd)
        enddo
     enddo
  enddo

end subroutine r_spiola_viscous
