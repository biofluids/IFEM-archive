!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     caculation of fp
!     
subroutine r_scalfp(fp,ocpp,i)
  use r_common, only: bpre,cpre,hp
  implicit none 

  real*8 :: fp,ocpp
  integer :: i

  fp=-ocpp*(bpre-cpre)*hp(i)
  return
end subroutine r_scalfp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     caculation of fu
!     
subroutine r_scalfu(fu,i,ni)
  use r_common
  implicit none

  real*8 :: fu
  integer :: i,ni

  integer :: m

  fu=0.0d0
  do m=1,6
     fu=fu+tos(m)*dge(m,i,ni)
  enddo
  
  return
end subroutine r_scalfu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!     caculation of kpp
!     
subroutine r_scalkpp(fkpp,ocpp,k,m)
  use r_common
  implicit none

  real*8 :: fkpp,ocpp
  integer :: k,m

  fkpp = ocpp*hp(k)*hp(m)

  return
end subroutine r_scalkpp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     caculation of kup
!     
subroutine r_scalkup(fkup,ocup,i,k,ni)
  use r_common
  implicit none
  
  real*8 :: fkup,ocup(6)
  integer :: i,k,ni

  integer :: m

  fkup=0.0d0
  do m=1,6
     fkup=fkup+ocup(m)*dge(m,i,ni)*hp(k)
  enddo
  
  return
end subroutine r_scalkup


