!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!     green-lagrangian strain and derivative
!     

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     dge(i,j,k)      -- i -- strain, j   -- dir,   k -- node
!    ddge(i,j,m,k,n)  -- i -- strain, j,m -- dir, k,n -- node
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine r_sstrain(toc,xto,lx,ly,lz,ne)
  use r_common
  implicit none

  real*8 :: toc(3,3),xto(3,3)
  integer,intent(in) :: lx,ly,lz,ne

  integer :: i,j,k,m,n

  do i=1,3
     ge(i,ne,lx,ly,lz) = 0.5d0*(toc(i,i)-1.0d0)
  enddo

  ge(4,ne,lx,ly,lz) = toc(2,3)
  ge(5,ne,lx,ly,lz) = toc(1,3)
  ge(6,ne,lx,ly,lz) = toc(1,2)
 
 !...for 3-D only
  do i=1,6
     do k=1,nis
	    do j=1,3
           if (i == 4) then
              dge(i,j,k)=xto(2,j)*bd(3,k)+xto(3,j)*bd(2,k)
           elseif (i == 5) then
              dge(i,j,k)=xto(1,j)*bd(3,k)+xto(3,j)*bd(1,k)
		   elseif (i == 6) then
              dge(i,j,k)=xto(1,j)*bd(2,k)+xto(2,j)*bd(1,k)
           elseif (i <= 3) then
              dge(i,j,k)=xto(i,j)*bd(i,k)
           endif

		   do m=1,3
			  do n=1,nis
			     if (m == j) then
                    if (i == 4) then
	                   ddge(i,j,m,k,n)=bd(2,n)*bd(3,k)+bd(3,n)*bd(2,k)
                    elseif (i == 5) then
				       ddge(i,j,m,k,n)=bd(1,n)*bd(3,k)+bd(3,n)*bd(1,k)
                    elseif (i == 6) then
				       ddge(i,j,m,k,n)=bd(1,n)*bd(2,k)+bd(2,n)*bd(1,k)
				    else
				       ddge(i,j,m,k,n)=bd(i,n)*bd(i,k)
                    endif
			     else
				    ddge(i,j,m,k,n) = 0.0d0
			     endif
			  enddo
		   enddo
		enddo
	 enddo
  enddo
  
  return
end subroutine r_sstrain



