!     
!     calculate deformation gradient
!     
subroutine r_stoxc(xto,xot,xj,xji,toxj,toxji,toc)
  use solid_variables, only: nsd_solid
  implicit none
      
  real(8) :: xto(nsd_solid,nsd_solid),xot(nsd_solid,nsd_solid)
  real(8) :: xj(nsd_solid,nsd_solid),xji(nsd_solid,nsd_solid)
  real(8) :: toxj(nsd_solid,nsd_solid),toxji(nsd_solid,nsd_solid)
  real(8) :: toc(nsd_solid,nsd_solid)
      
  integer :: i,j,k

  do i=1,nsd_solid
     do j=1,nsd_solid
        xto(i,j)=0.0d0
        xot(i,j)=0.0d0

        do k=1,nsd_solid
           xto(i,j)=xto(i,j)+xj(i,k)*toxji(k,j)
           xot(i,j)=xot(i,j)+toxj(i,k)*xji(k,j)
 ! if (i==2.and.j==1) write(*,*) 'xj,toxji=',xj(i,k)*toxji(k,j),xto(i,j)
       	enddo
			if (xto(i,j)<=1.0e-16) xto(i,j)=0.0
	enddo
  enddo
 ! write(*,*) 'xto(2,1)=',xto(2,1),xto(1,2)
  do i=1,nsd_solid
     do j=1,nsd_solid
        toc(i,j)=0.0d0
        do k=1,nsd_solid
           toc(i,j)=toc(i,j)+xto(k,i)*xto(k,j)

!			write(*,*) 'k,i,j=',k,i,j,xto(k,i),xto(k,j)

        enddo
     enddo
  enddo

  return
end subroutine r_stoxc
