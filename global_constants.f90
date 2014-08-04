module global_constants
  implicit none
  save

  real(8),parameter :: pi = 3.14159265d0
  integer,parameter :: xsd=1,ysd=2,zsd=3
  real(8),parameter :: TC=273.15, ZC=1.4, RC=2.87058e6
  real(8),parameter :: P0=1.01325e6, dens0=1.2922e-3
  real(8),parameter :: nDPML=10.0
  real(8),parameter,dimension(3) :: sigmaMaxPML=(/ 2.4e2, 2.4e2, 4.0 /)

end module global_constants