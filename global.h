      include "mpif.h"
      integer nec,nnc,nqdc,maxnec,maxnnc,maxnqdc,numproc,myid,nqdf
      common  nec,nnc,nqdc,maxnec,maxnnc,maxnqdc,numproc,myid,nqdf

      integer    ndfpad,nsdpad,nenpad,nquadpad
      parameter (ndfpad=5,nsdpad=3,nenpad=8,nquadpad=8)

      integer iquad, nquad      
      real* 8 sq(0:nsdpad,nenpad,nquadpad)
      real* 8 xq(nsdpad,nquadpad),wq(nquadpad)
      common  iquad,nquad,sq,xq,wq
      
      real* 8 tt,dt,t_start,alpha,res_g,del_g,res_l,del_l,turb_kappa       
      common  tt,dt,t_start,alpha,res_g,del_g,res_l,del_l,turb_kappa
      
      real* 8 ref_lgt,ref_vel,ref_den,vis_liq,vis_gas,den_liq,den_gas
      common  ref_lgt,ref_vel,ref_den,vis_liq,vis_gas,den_liq,den_gas
      
      real* 8 gas,liq,vmin,vmax,hmin,hmax
      common  gas,liq,vmin,vmax,hmin,hmax
      
      integer bc(ndfpad,21), bcf(21)
      real*8  bv(ndfpad,21),ic(ndfpad),bvf(21),icf
      common  bc,bcf,bv,ic,bvf,icf
      
      integer bcd(ndfpad,21)
      real*8  bvd(ndfpad,21),icd(ndfpad)
      common  bcd,bvd,icd
      
      real* 8 landa_over_mu
      common landa_over_mu

      integer surf(0:21),map(6,8,8)
      real* 8 interface(3),gravity(3),delta(0:21)
      common  surf,map,interface,gravity,delta

      integer hydro,etype,inner,outer,iscaling,kinner,kouter
      common  hydro,etype,inner,outer,iscaling,kinner,kouter
      
      logical hg_vol,static,taudt,restart,stokes,steady,conserve
      common  hg_vol,static,taudt,restart,stokes,steady,conserve
      
      logical twod
      common /twodimension/ twod

      logical calcforce
      integer fsurf(0:21),nfsurf
      common  calcforce,fsurf,nfsurf
      
      integer nn,ne,nqd,nq,nqf,nen,ndf,nsd,nrng,neface,nnface,maxconn
      common  nn,ne,nqd,nq,nqf,nen,ndf,nsd,nrng,neface,nnface,maxconn
      
      real*8 amp,spring,damper,mass
      common amp,spring,damper,mass

      integer its,iit,nts,nit,ntsbout,idisk,idiskmsh
      common  its,iit,nts,nit,ntsbout,idisk,idiskmsh
      
      integer tri,tet,qud,hex,tris,tets,quds,hexs
      parameter (tri=1,qud=2,tet=3,hex=4,tris=5,quds=6,tets=7,hexs=8)
            
      integer    xsd,ysd,zsd,udf,vdf,wdf,pdf
      parameter (xsd=1,ysd=2,zsd=3,udf=1,vdf=2,wdf=3,pdf=4)
      
      real* 8   pi,epsilon
      parameter (pi=3.141592654,epsilon=1.0e-6)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer nn_loc,nn_glb,nnon
      integer node_loc(1),node_glb(1)
	pointer (node_locptr,node_loc),(node_glbptr,node_glb)
      common  nn_loc,nn_glb,nnon,node_locptr,node_glbptr

      integer loc(1),locfrom(1),glb(1),glbfrom(1),des(1),src(1)
      pointer (locptr,loc),(locfromptr,locfrom)
      pointer (glbptr,glb),(glbfromptr,glbfrom)
      pointer (desptr,des),(srcptr,src)
      common  locptr,locfromptr,glbptr,glbfromptr,desptr,srcptr

      real*8  bufloc(1),bufglb(1)
      pointer (buflocptr,bufloc),(bufglbptr,bufglb)
      common  buflocptr,bufglbptr


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nn_on,nn_al,mmon
	integer node_on(1),node_al(1)
	pointer (node_onptr,node_on),(node_alptr,node_al)
	common  nn_on,nn_al,mmon,node_onptr,node_alptr

	integer on(1),onfrom(1),al(1),alfrom(1)
	integer destin(1),source(1)
	pointer (onptr,on),(onfromptr,onfrom)
	pointer (alptr,al),(alfromptr,alfrom)
	pointer (destinptr,destin),(sourceptr,source)
	common  onptr,onfromptr,alptr,alfromptr,destinptr,sourceptr

	real*8  bufon(1),bufal(1)
	pointer (bufonptr,bufon),(bufalptr,bufal)
	common  bufonptr,bufalptr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



