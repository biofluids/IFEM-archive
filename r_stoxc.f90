!     
!     calculate deformation gradient
!     
subroutine r_stoxc(xto,xot,xj,xji,toxj,toxji,toc)
  implicit none
      
  real(8) :: xto(3,3),xot(3,3),xj(3,3),xji(3,3),toxj(3,3),toxji(3,3),toc(3,3)
      
  integer :: i,j,k

  do i=1,3
     do j=1,3
        xto(i,j)=0.0d0
        xot(i,j)=0.0d0
        do k=1,3
           xto(i,j)=xto(i,j)+xj(i,k)*toxji(k,j)
           xot(i,j)=xot(i,j)+toxj(i,k)*xji(k,j)
        enddo
     enddo
  enddo

  do i=1,3
     do j=1,3
        toc(i,j)=0.0d0
        do k=1,3
           toc(i,j)=toc(i,j)+xto(k,i)*xto(k,j)
        enddo
     enddo
  enddo
  !write(*,*) 'toc=',toc(1,1),toc(2,2),toc(3,3),toc(1,2)
 
  return
end subroutine r_stoxc
