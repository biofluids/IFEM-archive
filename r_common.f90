!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains all variables related to the FEM rubber model (Xiaodong Wang, PU)
! 
! 22.April.2003,  Axel Gerstenberger

module r_common
  implicit none
  save


  integer,parameter :: nup=4,nnumr=4
  integer,parameter :: mnos=0,nels=0
  integer,parameter :: mno=100000,mno2=3*mno,nel=106454
  real(8),parameter :: x13  = 1.0d0/3.0d0
  real(8),parameter :: x23  = 2.0d0/3.0d0
  real(8),parameter :: x43  = 4.0d0/3.0d0
  real(8),parameter :: x53  = 5.0d0/3.0d0
  real(8),parameter :: x73  = 7.0d0/3.0d0
  real(8),parameter :: x83  = 8.0d0/3.0d0
  real(8),parameter :: x49  = 4.0d0/9.0d0
  real(8),parameter :: x109 = 10.0d0/9.0d0
  real(8) :: xk,xtedis,xvisc,xviss,xstretch
  integer :: nump,nrigid
  integer :: nbouc,numgb,numeb,numfn
  real(8) :: tfun(10)
  real(8) :: PK2str(6),bPK2str(6)
  real(8) :: PK1str_tens(3,3)
  real(8) :: cpre,bpre
  real(8) :: h(9),hp(6),dbpre(6),ddbpre(6,6),r_p(3,9)
  real(8) :: bd(3,9),bd_curr(3,9)
  real(8) :: fnod(mno,3),boup(mno,3)
  integer :: nodefn(mno),ndirfn(mno)
  integer :: nface(nel)
  real(8) :: du(3,mno)
  real(8) :: predrf(mno2)
  real(8) :: fnodo(mno,3),boupo(mno,3)
  integer :: ntfun
  integer :: ninit,npr,npdis,ntprint,nprestress
  integer :: nfuns(10)
  integer :: nchkread
  real(8) :: prec(nup*nel),density_solid
  integer :: nreact,initdir,nrtp
  real(8) :: vnorm,fnorm
  real(8) :: xmg(3)
  real(8) :: dge(6,3,9),ddge(6,3,3,9,9)
  real(8) :: young_mod, Poisson,rc1,rc2,rk
  integer :: material_type
end module r_common

