!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     calculation cauchy stress
!     
subroutine r_scauchy(det,todet,xto,lx,ly,lz,ne)
  use r_common
  implicit none 

  real*8,intent(in) :: det,todet
  real*8 xto(3,3)
  integer,intent(in) :: lx,ly,lz,ne
  
  real*8 :: ssb(3,3),ss(3,3),tt(3,3)
  integer :: i,j,m,n

  ss(1,1)=tos(1)
  ss(2,2)=tos(2)
  ss(1,3)=tos(5)
  ss(3,1)=tos(5)
  ss(2,3)=tos(4)
  ss(3,2)=tos(4)
  ss(1,2)=tos(6)
  ss(2,1)=tos(6)
  ss(3,3)=tos(3)

  do i=1,6
     cstr(i,ne,lx,ly,lz)=0.0d0
  enddo

  do i=1,3
     do j=1,3
		tt(i,j)=xto(i,j)
     enddo
  enddo

  do m=1,3
     do n=1,3
        cstr(1,ne,lx,ly,lz)=cstr(1,ne,lx,ly,lz) + todet/det*tt(1,m)*ss(m,n)*tt(1,n)
        cstr(2,ne,lx,ly,lz)=cstr(2,ne,lx,ly,lz) + todet/det*tt(2,m)*ss(m,n)*tt(2,n)
        cstr(3,ne,lx,ly,lz)=cstr(3,ne,lx,ly,lz) + todet/det*tt(3,m)*ss(m,n)*tt(3,n)
        cstr(4,ne,lx,ly,lz)=cstr(4,ne,lx,ly,lz) + todet/det*tt(3,m)*ss(m,n)*tt(2,n)
	    cstr(5,ne,lx,ly,lz)=cstr(5,ne,lx,ly,lz) + todet/det*tt(1,m)*ss(m,n)*tt(3,n)
        cstr(6,ne,lx,ly,lz)=cstr(6,ne,lx,ly,lz) + todet/det*tt(1,m)*ss(m,n)*tt(2,n)
	 enddo
  enddo

  return
end subroutine r_scauchy




