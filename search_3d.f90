!
!  search_3d.f90
!  
!
!  Created by Xingshi on 1/30/08.
!  Copyright 2008 __MyCompanyName__. All rights reserved.
!

       subroutine search_3d(finf, x, xna, nn, nsd, ne, nen, ien, inf, ninf,maxconn)        
           implicit none
	   integer :: ninf, nn, nsd, ne, nen, maxconn ! Fluid variable from main
	   integer ien(nen,ne) ! transfer matrix connectted local information to global
	   integer inf(maxconn) ! subroutine variable to store the index of element which contains the solid node
	   real(8) x(nsd), xna(nsd, nn) ! x - coordinates of solid node   xna - coordinates of fluid nodes in global view
           real(8) tetrahedran(3,4)
	   integer answer
	   integer answer_hex
	   real(8) sub_tra(3,4)
	   integer finf
	   integer i, j
	   
		answer=0   	   
	   do j=1,ninf
	      
	   
	   if (nen==4) then ! for terahedran
!	   write(*,*) 'I am now in the 3d search sub'
           do i=1,nen
	   tetrahedran(1:nsd,i)=xna(1:nsd, ien(i,inf(j)))
	   end do
	   call check_tra(tetrahedran,x,answer)
	   if (answer==1) goto 8888
	  	   
	   
	   
	   
	   
	      
	  
		   	   
	   else if (nen==8) then ! for hexahedral
	   	                    ! the mesh number has to be specific declare
	   ! the hexahedral is devided to six tetrahendral and do six times
	   ! if the point is in the hexahedral one of these six tertrahendarl should contain this point.
	   sub_tra(1:nsd,1)=xna(1:nsd,ien(1,inf(j)))
	   sub_tra(1:nsd,2)=xna(1:nsd,ien(3,inf(j)))
	   sub_tra(1:nsd,3)=xna(1:nsd,ien(4,inf(j)))
	   sub_tra(1:nsd,4)=xna(1:nsd,ien(7,inf(j)))
	   call check_tra(sub_tra,x,answer_hex)
	   if (answer_hex==1) goto 8888
	   sub_tra(1:nsd,1)=xna(1:nsd,ien(1,inf(j)))
	   sub_tra(1:nsd,2)=xna(1:nsd,ien(4,inf(j)))
	   sub_tra(1:nsd,3)=xna(1:nsd,ien(7,inf(j)))
	   sub_tra(1:nsd,4)=xna(1:nsd,ien(8,inf(j)))
	   call check_tra(sub_tra,x,answer_hex)
	   if (answer_hex==1) goto 8888
	   sub_tra(1:nsd,1)=xna(1:nsd,ien(1,inf(j)))
	   sub_tra(1:nsd,2)=xna(1:nsd,ien(2,inf(j)))
	   sub_tra(1:nsd,3)=xna(1:nsd,ien(3,inf(j)))
	   sub_tra(1:nsd,4)=xna(1:nsd,ien(6,inf(j)))
	   call check_tra(sub_tra,x,answer_hex)
	   if (answer_hex==1) goto 8888
	   sub_tra(1:nsd,1)=xna(1:nsd,ien(1,inf(j)))
	   sub_tra(1:nsd,2)=xna(1:nsd,ien(3,inf(j)))
	   sub_tra(1:nsd,3)=xna(1:nsd,ien(6,inf(j)))
	   sub_tra(1:nsd,4)=xna(1:nsd,ien(7,inf(j)))
	   call check_tra(sub_tra,x,answer_hex)
	   if (answer_hex==1) goto 8888
	   sub_tra(1:nsd,1)=xna(1:nsd,ien(1,inf(j)))
	   sub_tra(1:nsd,2)=xna(1:nsd,ien(5,inf(j)))
	   sub_tra(1:nsd,3)=xna(1:nsd,ien(6,inf(j)))
	   sub_tra(1:nsd,4)=xna(1:nsd,ien(8,inf(j)))
	   call check_tra(sub_tra,x,answer_hex)
	   if (answer_hex==1) goto 8888
	   sub_tra(1:nsd,1)=xna(1:nsd,ien(1,inf(j)))
	   sub_tra(1:nsd,2)=xna(1:nsd,ien(6,inf(j)))
	   sub_tra(1:nsd,3)=xna(1:nsd,ien(7,inf(j)))
	   sub_tra(1:nsd,4)=xna(1:nsd,ien(8,inf(j)))
	   call check_tra(sub_tra,x,answer_hex)
	   if (answer_hex==1) goto 8888
	   
				  
			
	   
	   else
           write(*,*) 'Sorry right now we can do only for terahedran or hexahedral'
	   end if
	   
	   
	   end do
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
8888 continue
        finf=inf(j)	   
	   return
       end
	   
subroutine check_tra(tet,x,answer)  ! subroutine to check if a point is in the tetrahendal
real(8) tet(3,4)
real(8) x(3)
real(8) matrix(4,4)
real(8) det_tra(5)
integer i
integer answer
integer nsd

nsd=3
answer=0


          do i=1,4
	      matrix(i,4)=1
		  end do
       matrix(1,1:nsd)=tet(1:3,1)
	   matrix(2,1:nsd)=tet(1:3,2)
	   matrix(3,1:nsd)=tet(1:3,3)
	   matrix(4,1:nsd)=tet(1:3,4)
	   call determinant(matrix,4,4,det_tra(5))
	   
	      do i=1,4
	      matrix(i,4)=1
		  end do
	   matrix(1,1:nsd)=x(1:3)
	   matrix(2,1:nsd)=tet(1:3,2)
	   matrix(3,1:nsd)=tet(1:3,3)
	   matrix(4,1:nsd)=tet(1:3,4)
	   call determinant(matrix,4,4,det_tra(1))! subroutine to calculate the determinant of a matrix
	   
	      do i=1,4
	      matrix(i,4)=1
		  end do
	   matrix(2,1:nsd)=x(1:3)
	   matrix(1,1:nsd)=tet(1:3,1)
	   matrix(3,1:nsd)=tet(1:3,3)
	   matrix(4,1:nsd)=tet(1:3,4)
	   call determinant(matrix,4,4,det_tra(2))
	     
		  do i=1,4
	      matrix(i,4)=1
		  end do
	   matrix(3,1:nsd)=x(1:3)
	   matrix(1,1:nsd)=tet(1:3,1)
	   matrix(2,1:nsd)=tet(1:3,2)
	   matrix(4,1:nsd)=tet(1:3,4)
	   call determinant(matrix,4,4,det_tra(3))
	   
	      do i=1,4
	      matrix(i,4)=1
		  end do
	   matrix(4,1:nsd)=x(1:3)
	   matrix(1,1:nsd)=tet(1:3,1)
	   matrix(2,1:nsd)=tet(1:3,2)
	   matrix(3,1:nsd)=tet(1:3,3)
	   call determinant(matrix,4,4,det_tra(4))
	   
	       if (  det_tra(5)*det_tra(1)>0 .and. &
	             det_tra(5)*det_tra(2)>0 .and. &
		     det_tra(5)*det_tra(3)>0 .and. &
		     det_tra(5)*det_tra(4)>0) then
		         answer=1
		   end if
		   return
		   end
