module global_constants
    implicit none
    save

    real(8),parameter :: pi = 3.14159265d0
    integer,parameter :: xsd=1,ysd=2,zsd=3
    real(8),parameter :: TC=273.15, ZC=1.4, RC=2.87058e6
    real(8),parameter :: P0=1.01325e6, dens0=1.2922e-3
    real(8),parameter :: nDPML=41.0     ! the width of PML layers are determined by nDPML*(max size of elements)
    real(8),parameter,dimension(3) :: sigmaMaxPML=(/ 2.4e6, 2.4e6, 2.4e6 /)  ! \sigma_{max} for x/y/z-directions
    real(8),parameter,dimension(4) :: uvwpAvgPML=(/ 0.0, 0.0, 0.0, 0.0 /)    ! pseudo mean flow velocities

end module global_constants
