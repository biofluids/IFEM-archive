!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     calculation of fp
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
!     calculation of fu
!     
subroutine r_scalfu(fu,isd,ni)
  use r_common
  implicit none

  real*8 :: fu
  integer :: isd,ni

  !integer :: m
  integer :: ksd

  fu=0.0d0
  !do m=1,6
  !   fu=fu+PK2str(m)*dge(m,isd,ni)
  !enddo

  do ksd = 1,3
    fu = fu + bd(ksd,ni)*PK1str_tens(ksd,isd)
  enddo

  return
end subroutine r_scalfu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     calculation of fu in current configuration
!     
subroutine r_scalfu_curr(fu,i,ni,ine,iq)
  use r_common
  implicit none

  real*8 :: fu
  integer :: i,ni,ine,iq

  !integer :: m

  fu=0.0d0
  !do m=1,6
     !fu=fu+PK2str(m)*dge(m,i,ni)
  !enddo
  if (i == 1) then
	 fu=fu+cstr(1,ine,iq)*bd_curr(1,ni)
	 fu=fu+cstr(6,ine,iq)*bd_curr(2,ni)
	 fu=fu+cstr(5,ine,iq)*bd_curr(3,ni)
  elseif (i == 2) then
	 fu=fu+cstr(6,ine,iq)*bd_curr(1,ni)
	 fu=fu+cstr(2,ine,iq)*bd_curr(2,ni)
	 fu=fu+cstr(4,ine,iq)*bd_curr(3,ni)
  elseif (i == 3) then
	 fu=fu+cstr(5,ine,iq)*bd_curr(1,ni)
	 fu=fu+cstr(4,ine,iq)*bd_curr(2,ni)
	 fu=fu+cstr(3,ine,iq)*bd_curr(3,ni)
  endif
  
  return
end subroutine r_scalfu_curr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!     calculation of kpp
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
!     calculation of kup
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


