subroutine initialize
  use run_variables
  use restart_lib, only: restart
  use fluid_variables
  implicit none

  call facemap

  alpha = 1.0
  tt = 0.0
  dt = 1.0
  t_start = 0.0
  ref_lgt = 1.0
  ref_vel = 1.0
  ref_den = 1.0
  vis_liq = 1.0
  den_liq = 1.0
  vis_gas = 1.0
  den_gas = 1.0

  gravity(1:nsd) = 0.0d0
  interface(1:nsd) = -999

  turb_kappa = 0.41

  hydro = 0 

  surf(0:maxnsurf) = 0

  ndf = 4
  nsd = 3
  nn = 0
  ne = 0
  nq = 0
  nrng = 0

  iquad = 2
  nts = 1
  nit = 1
  ntsbout = 1
  idisk = 0

  !nfsurf = 0
  !fsurf(:) = 0
  calcforce = .false.

  inner = 1
  kinner = 1
  outer = 1
  kouter = 1
  iscaling = 1

  restart  = 0
  stokes  = .false.
  steady  = .false.
  hg_vol  = .false.
  static  = .false.
  taudt  = .false.
  conserve = .false.


  bc(1:ndfpad,1:maxnsurf) = 0
  bv(1:ndfpad,1:maxnsurf) = -999.0

  bcd(1:nsdpad,1:maxnsurf) = 0
  bvd(1:nsdpad,1:maxnsurf) = -999.0

  ic(1:ndfpad) = -999.0

  delta(0:21) = 0.0

  !mass = 1.0
 
end subroutine initialize


