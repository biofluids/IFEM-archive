!     
!     calculate deformation gradient
!     
subroutine r_stoxc(xto,xot,xj,xji,toxj,toxji,toc,ne)
  use solid_variables, only: nsd_solid
  implicit none
      
  real(8) :: xto(nsd_solid,nsd_solid),xot(nsd_solid,nsd_solid)
  real(8) :: xj(nsd_solid,nsd_solid),xji(nsd_solid,nsd_solid)
  real(8) :: toxj(nsd_solid,nsd_solid),toxji(nsd_solid,nsd_solid)
  real(8) :: toc(nsd_solid,nsd_solid)
  integer :: xtoI(nsd_solid,nsd_solid)
  integer :: xotI(nsd_solid,nsd_solid)
  integer :: tocI(nsd_solid,nsd_solid)
  integer :: i,j,k,ne

  do i=1,nsd_solid
     do j=1,nsd_solid
        xto(i,j)=0.0d0
        xot(i,j)=0.0d0

        do k=1,nsd_solid
           xto(i,j)=xto(i,j)+xj(i,k)*toxji(k,j)
           xot(i,j)=xot(i,j)+toxj(i,k)*xji(k,j)
	   	enddo
	
	!	xto(i,j)=xto(i,j)*1.0d5
	!	xtoI(i,j)=xto(i,j)
	!	xto(i,j)=xtoI(i,j)*1.0d-5
	!	if (xto(i,j)<=1.0e-16) xto(i,j)=0.0
		
	!	xot(i,j)=xot(i,j)*1.0d5
	!	xotI(i,j)=xot(i,j)
	!	xot(i,j)=xotI(i,j)*1.0d-5
	!	if (xot(i,j)<=1.0e-16) xot(i,j)=0.0
	enddo
  enddo
  do i=1,nsd_solid
     do j=1,nsd_solid
        toc(i,j)=0.0d0
        do k=1,nsd_solid
           toc(i,j)=toc(i,j)+xto(k,i)*xto(k,j)

        enddo

	!	toc(i,j)=toc(i,j)*1.0d5
	!	tocI(i,j)=toc(i,j)
	!	toc(i,j)=tocI(i,j)*1.0d-5
	!	if (toc(i,j)<=1.0e-16) toc(i,j)=0.0
     enddo
  enddo

  return
end subroutine r_stoxc
