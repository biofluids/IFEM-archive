!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     pressure interpolation
!
subroutine r_spress(rs,ne)
  use r_common
  implicit none 

  real*8 :: rs(3)
  integer,intent(in) :: ne

  real*8 :: r,s,t
  integer :: i,ntt

  r=rs(1)
  s=rs(2)
  t=rs(3)
 
  hp(1)=1.0d0
  hp(2)=r
  hp(3)=s
  hp(4)=t
  cpre=0.0d0
  do i=1,nump
     ntt=(ne-1)*nump+i
     cpre=cpre+prec(ntt)*hp(i)
  enddo

  return
end subroutine r_spress

