module fluid_variables
  implicit none
  save

  integer nec,nnc,nqdc,nqdf,maxnec,maxnnc,maxnqdc,numproc,myid

  integer,parameter :: ndfpad=5,nsdpad=3,nenpad=8,nquadpad=12

  integer iquad, nquad, nquad2d      
  real* 8 sq(0:nsdpad,nenpad,nquadpad)
  real* 8 xq(nsdpad,nquadpad),wq(nquadpad)
  real* 8 sq2d(0:nsdpad,nenpad,nquadpad*6)
  real* 8 xq2d(nsdpad,nquadpad*6),wq2d(nquadpad*6)


  real* 8 t_start,alpha,res_g,del_g,res_l,del_l,turb_kappa       

  real* 8 ref_lgt,ref_vel,ref_den,vis_liq,vis_gas,den_liq,den_gas

  real* 8 gas,liq,vmin,vmax,hmin,hmax

  integer bc(ndfpad,21), bcf(21)
  real*8  bv(ndfpad,21),ic(ndfpad),bvf(21),icf

  integer surf(0:21),mapping(6,8,8)
  real* 8 interface(3),gravity(3),delta(0:21)

  integer hydro,etype,inner,outer,iscaling

  logical hg_vol,static,taudt,restart,stokes,steady,conserve

  logical twod

  logical calcforce
  integer fsurf(0:21),nfsurf

  integer nn,ne,nqd,nq,nqf,nen,ndf,nsd,nrng,neface,nnface

  integer iit,nit,idisk

  integer,parameter :: tri=1,qud=2,tet=3,hex=4,tris=5,quds=6,tets=7,hexs=8
  integer,parameter :: udf=1,vdf=2,wdf=3,pdf=4

end module fluid_variables
