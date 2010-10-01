module run_variables
  implicit none
  save

  integer :: its,nts,nts_start
  real(8) :: tt,dt,T
  integer :: ntsbout,N
  integer :: restart_freq,restart_klok,restart_unit
  integer,parameter :: restart_u1 = 9611, restart_u2 = 9612
  real*8 com_kappa ! magic number for PS
end module run_variables
