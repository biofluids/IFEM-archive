module fluid_variables
  implicit none
  save

  integer,parameter :: ndfpad=5,nsdpad=3,nenpad=8,nquadpad=8

  integer :: iquad, nquad, nquad2d      
  real(8) :: sq(0:nsdpad,nenpad,nquadpad)
  real(8) :: xq(nsdpad,nquadpad),wq(nquadpad)
  real(8) :: sq2d(0:nsdpad,nenpad,nquadpad*6)
  real(8) :: xq2d(nsdpad,nquadpad*6),wq2d(nquadpad*6)


  real(8) :: alpha,res_g,del_g,res_l,del_l,turb_kappa

  real(8) :: ref_lgt,ref_vel,ref_den,vis_liq,vis_gas,den_liq,den_gas

  real(8) :: gas,liq,vmin,vmax,hmin,hmax

  integer,parameter :: maxnsurf = 21

  integer :: bc(ndfpad,maxnsurf) !,bcf(maxnsurf)
  real(8) :: bv(ndfpad,maxnsurf),ic(ndfpad)  !,bvf(maxnsurf),icf


  integer :: bcd(nsdpad,maxnsurf)
  real(8) :: bvd(nsdpad,maxnsurf),icd(nsdpad)


  real(8) :: landa_over_mu



  integer :: surf(0:maxnsurf),mapping(6,8,8)
  real(8) :: interface(1:nsdpad),gravity(1:nsdpad),delta(0:21)

  integer :: hydro,etype,inner,outer,iscaling,kinner,kouter

  logical :: hg_vol,static,taudt,stokes,steady,conserve

  logical :: twod

  logical :: calcforce
  !integer fsurf(0:21),nfsurf

  integer :: nn,ne,nq,nen,ndf,nsd,nrng,neface,nnface
  integer :: ne_surf
  integer :: nen_surf

  integer :: iit,nit,idisk

  integer,parameter :: tri=1,qud=2,tet=3,hex=4,tris=5,quds=6,tets=7,hexs=8
  integer,parameter :: udf=1,vdf=2,wdf=3,pdf=4

  integer :: ale_mesh_update

  real(8),parameter :: epsilon = 1.0e-13
  !real(8),parameter :: epsilon = 1.0e-6

  integer,parameter :: fluid_mesh_filename_length = 15
  character(len = fluid_mesh_filename_length),parameter :: fluid_mesh_filename = "input_fluid.msh"

end module fluid_variables
