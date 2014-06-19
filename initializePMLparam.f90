!==================================================================
! J. Yang
! Rensselaer Polytechnic Institute
! This subroutine initializes the parameters for Perfectly Matched Layer
!------------------------------------------------------------------
! Latest update: Jack, 05/23/2014
!==================================================================
subroutine initializePMLparam(xloc,sigmaPML)
  use global_constants
  use fluid_variables
  implicit none

  real(8) :: sigmaPML(nsd,nn)
  real(8) :: xloc(nsd,nn)
  real(8) :: DPML, distPML
  integer :: irng, inn, isd

    sigmaPML(:,:) = 0.0
    DPML=nDPML*hmax

    if (sumNbcPML .ne. 0) then
        do irng=1,nrng
            isd=flagPML(irng)
            if (isd .ne. 0) then
                do inn=1,nn
                    distPML=abs( xloc(isd,inn)-xyzcPML(irng) )
                    if ( distPML < DPML ) then
                        ! alphaPML = 1.0 + 2*(1.0 - distPML/DPML)**2
                        ! sigmaPML(isd,inn) = sigmaMaxPML(isd) * (1.0 - distPML/DPML)**alphaPML
                        sigmaPML(isd,inn) = sigmaMaxPML(isd) * (1.0 - distPML/DPML)**4.0
                    endif
                enddo
            endif
        enddo
    endif

return
end subroutine initializePMLparam








