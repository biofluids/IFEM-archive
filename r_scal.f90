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
  !do m=1,6
  !   fu=fu+PK2str(m)*dge(m,i,ni)
  !enddo

  if (i == 1) then
	 fu=fu+PK1str(1)*bd(1,ni)
	 fu=fu+PK1str(6)*bd(2,ni)
	 fu=fu+PK1str(5)*bd(3,ni)
  elseif (i == 2) then
	 fu=fu+PK1str(6)*bd(1,ni)
	 fu=fu+PK1str(2)*bd(2,ni)
	 fu=fu+PK1str(4)*bd(3,ni)
  elseif (i == 3) then
	 fu=fu+PK1str(5)*bd(1,ni)
	 fu=fu+PK1str(4)*bd(2,ni)
	 fu=fu+PK1str(3)*bd(3,ni)
  endif



  
  return
end subroutine r_scalfu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     caculation of fu in current configuration
!     
subroutine r_scalfu_curr(fu,i,ni,ine,lx,ly,lz)
  use r_common
  implicit none

  real*8 :: fu
  integer :: i,ni,ine,lx,ly,lz

  integer :: m

  fu=0.0d0
  !do m=1,6
     !fu=fu+PK2str(m)*dge(m,i,ni)

  if (i == 1) then
	 fu=fu+cstr(1,ine,lx,ly,lz)*bd_curr(1,ni)
	 fu=fu+cstr(6,ine,lx,ly,lz)*bd_curr(2,ni)
	 fu=fu+cstr(5,ine,lx,ly,lz)*bd_curr(3,ni)
  elseif (i == 2) then
	 fu=fu+cstr(6,ine,lx,ly,lz)*bd_curr(1,ni)
	 fu=fu+cstr(2,ine,lx,ly,lz)*bd_curr(2,ni)
	 fu=fu+cstr(4,ine,lx,ly,lz)*bd_curr(3,ni)
  elseif (i == 3) then
	 fu=fu+cstr(5,ine,lx,ly,lz)*bd_curr(1,ni)
	 fu=fu+cstr(4,ine,lx,ly,lz)*bd_curr(2,ni)
	 fu=fu+cstr(3,ine,lx,ly,lz)*bd_curr(3,ni)
  endif
  !enddo
  
  return
end subroutine r_scalfu_curr


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


