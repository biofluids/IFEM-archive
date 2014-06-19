!==================================================================
! J. Yang
! Rensselaer Polytechnic Institute
! This subroutine calculates the PML auxiliary variable "qv"
!------------------------------------------------------------------
! Latest update: Jack, 05/21/2014
!==================================================================
subroutine evaluateQvPML(qv,dUnkn)
  use global_constants
  use fluid_variables, only: ndf,nsd,udf,vdf,wdf,nn
  use run_variables, only: dt
  implicit none

  real(8) :: qv(ndf,nn), dUnkn(ndf,nn)
  real(8) :: den_local(nn)
  integer :: idof,inn

!------------------------------------------------------------------
! Jack, 05/23/2014
! only dealing with real fluid for now
! \partial q^u_i / \partial t = \rho * u_i - \rho_0 * u_i^0
! \partial q^p / \partial t = p - p_0 = \hat(p)
! \partial q^\rho / \partial t = \rho - \rho_0 = \hat(\rho) ! not calculated
! q^\rho = q^p / c_s^2                                      ! not calculated

    do inn=1,nn
        den_local(inn)=( 1.0 + dUnkn(ndf,inn)/(ZC*P0) )*dens0
        do idof=1,nsd
            qv(idof,inn)=qv(idof,inn)+den_local(inn)*dUnkn(idof,inn)*dt
        enddo
        qv(ndf,inn)=qv(ndf,inn)+dUnkn(ndf,inn)*dt
    enddo
!-------------------------------------------------------------------

return
end subroutine evaluateQvPML







