!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains all variables related to the FEM rubber model (Xiaodong Wang, PU
! 
! 22.April.2003,  Axel Gerstenberger

module r_common
  implicit none
  save



  integer,parameter :: nup=4,nnis=9,nnumr=4
  integer,parameter :: mnos=0,nels=0
  integer,parameter :: mno=23000,mno2=3*mno,nel=15000
  integer,parameter :: n_solid_max=3
!ccccccccccccccccccccccccccccccccccccccccccccccc
  integer :: nxt(3,mno)
  real*8  :: xk,xtedis,xvisc,xviss,xstretch
  real*8  :: xg(4,4),wgt(4,4)
  real*8  :: xg_tetra(4,4),wgt_tetra(4,4)
  integer :: nis,nump,nrigid
  integer :: nbouc,numgb,numeb,numfn
  integer :: nint
!  integer :: iti
  real*8  :: tfun(10)
  real*8  :: xkup(3*nnis,nup,nel)
  real*8  :: xkpp(nup,nup,nel),xfp(nup,nel)
  real*8  :: shift(3,n_solid_max)
  real*8  :: ge(6,nel,3,3,3),cstr(6,nel,3,3,3)  !...Cauchy stress
   real*8 :: tstress(6,mno),tstrain(6,mno)
  real*8  :: tos(6),btos(6)
  real*8  :: cpre,bpre
  real*8  :: h(9),hp(6),dbpre(6),ddbpre(6,6),r_p(3,9),bd(3,9)
  !integer :: ifc(30)
  real*8  :: fnod(mno,3),boup(mno,3)
  integer :: nodefn(mno),ndirfn(mno)
  integer :: nbe(nel),nbn(nel,3),nface(nel),numdir(6),ndirgb(6),nodegb(6,mno)
!  real*8  :: diso(3,mno)
  real*8  :: du(3,mno)!,dvel(3,mno)
  !real*8  :: acm(3,mno),acm1(3,mno)
  real*8  :: drf(mno2),predrf(mno2)
  real*8  :: drf2(mno2),predrf2(mno2)
  real*8  :: y(3,9)
!  real*8  :: coor(mno,3)
  integer :: nea(nel,9)!,neap(nel,4)
  real*8  :: dge(6,3,9),ddge(6,3,3,9,9)
  real*8  :: pre(nup,nel)
!  integer :: n1link(mno),n2link(mno)
!  integer :: nod2link(4,mno)!,nod1link(4,mno)
!  real*8  :: coe1link(4,mno)!,coe2link(4,mno)
  real*8  :: rc1,rc2,rk
  
!ccccccccccccccccccccccccc
  real*8  :: fnodo(mno,3),boupo(mno,3)
  integer :: ntfun
  integer :: ninit,npr,npdis,ntprint,nprint(10)
  integer :: nfuns(10)
!  real*8  :: xome(10),fbacc(3)      
!  integer :: imasdir(mno)!,islavdir(mno)
!  integer :: numint(mno),ninsk(mno)
!  real*8  :: xang(mno)
!  integer :: nodesl(mno),nodema(mno)
!  real*8  :: amct(mno)
  integer :: nchkread,nndtem
  real*8  :: prec(nup*nel),density_solid
!  real*8  :: raf(mno,3)   !,xindis(3,mno)
  integer :: nreact,initdir,nrtp,nraf(mno2),ndraf(mno2)
  integer :: nina,ndprint(mno)
  real*8  :: x13,x23,x43,x53,x73,x83,x49,x109
!  real*8  :: jump_frame,jump_fd,jump_mk,mk_start,mk_finish
!  integer :: npoint_fiber(0:40000),nfiber_file
  real*8  :: vnorm,fnorm,vtol,ftol,alpha_solid,beta_solid
!  integer :: nmove,nbou
  real*8  :: xmg1,xmg2,xmg3

end module r_common

