!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     calculation cauchy stress
!     
subroutine r_scauchy(det,todet,xto,iq,ne)
  use r_common
  implicit none 

  real*8,intent(in) :: det,todet
  real*8 xto(3,3)
  integer,intent(in) :: iq,ne
  
  !real*8 :: ssb(3,3)
  real*8 :: ss(3,3)
  integer :: isd,jsd

  ss(1,1)=PK2str(1)
  ss(2,2)=PK2str(2)
  ss(1,3)=PK2str(5)
  ss(3,1)=PK2str(5)
  ss(2,3)=PK2str(4)
  ss(3,2)=PK2str(4)
  ss(1,2)=PK2str(6)
  ss(2,1)=PK2str(6)
  ss(3,3)=PK2str(3)

  cstr(1:6,ne,iq) = 0.0d0
  
  do isd=1,3
     do jsd=1,3
        cstr(1,ne,iq)=cstr(1,ne,iq) + todet/det*xto(1,isd)*ss(isd,jsd)*xto(1,jsd)
        cstr(2,ne,iq)=cstr(2,ne,iq) + todet/det*xto(2,isd)*ss(isd,jsd)*xto(2,jsd)
        cstr(3,ne,iq)=cstr(3,ne,iq) + todet/det*xto(3,isd)*ss(isd,jsd)*xto(3,jsd)
        cstr(4,ne,iq)=cstr(4,ne,iq) + todet/det*xto(3,isd)*ss(isd,jsd)*xto(2,jsd)
	    cstr(5,ne,iq)=cstr(5,ne,iq) + todet/det*xto(1,isd)*ss(isd,jsd)*xto(3,jsd)
        cstr(6,ne,iq)=cstr(6,ne,iq) + todet/det*xto(1,isd)*ss(isd,jsd)*xto(2,jsd)
	 enddo
  enddo

  return
end subroutine r_scauchy




