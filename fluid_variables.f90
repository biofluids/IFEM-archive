module fluid_variables
  implicit none
  save

  integer,parameter :: ndfpad=5,nsdpad=3,nenpad=8,nquadpad=8
  integer :: iquad, nquad, nquad2d      
  real(8) :: sq(0:nsdpad,nenpad,nquadpad)
  real(8) :: xq(nsdpad,nquadpad),wq(nquadpad)
  real(8) :: sq2d(0:nsdpad,nenpad,nquadpad*6)
  real(8) :: xq2d(nsdpad,nquadpad*6),wq2d(nquadpad*6)
  integer :: maxconn
  real(8) :: t_start,alpha,res_g,del_g,res_l,del_l,turb_kappa
  real(8) :: ref_lgt,ref_vel,ref_den,vis_liq,den_liq
  real(8) :: liq,vmin,vmax,hmin,hmax
  integer,parameter :: maxnsurf = 21
  integer :: bc(ndfpad,maxnsurf) 
  real(8) :: bv(ndfpad,maxnsurf),ic(ndfpad)  
  integer :: bcd(nsdpad,maxnsurf)
  real(8) :: bvd(nsdpad,maxnsurf),icd(nsdpad)
  real(8) :: landa_over_mu
  integer :: mapping(6,8,8)
  real(8) :: interface(1:nsdpad),gravity(1:nsdpad),delta(0:21)
  integer :: etype,inner,outer,iscaling,kinner,kouter
  logical :: hg_vol,static,taudt,stokes,steady,conserve
  integer :: restart
  logical :: twod
  logical :: calcforce
  integer :: nn,ne,nq,nen,ndf,nsd,nrng,neface,nnface
  integer :: iit,nit,idisk
  integer,parameter :: tri=1,qud=2,tet=3,hex=4,tris=5,quds=6,tets=7,hexs=8
  integer :: udf,vdf,wdf,pdf
  real(8),parameter :: epsilon = 1.0e-12
  integer :: ne_inflow, edge_inflow ! Neumann "pressure" boundary condition
  real(8) pin                       ! inflow pressure
  integer ptotflag ! flag for use total pressure BC or not 1---> yes, 0 ---> no
  integer outedge
  integer :: flagPML(maxnsurf)      ! flag if PML is defined on this boundary
  integer :: sumNbcPML
  real(8) :: xyzcPML(maxnsurf)      ! coordinate of the maximum coordinate
end module fluid_variables
