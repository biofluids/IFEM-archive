!
!  ludcmp.f90
!  
!
!  Created by Xingshi on 1/30/08.
!  Copyright 2008 __MyCompanyName__. All rights reserved.
!  a the matrix put in , n=np give the dimension, d returns the determinant
	   subroutine determinant(a,n,np,d)! LU decomposition to calculate dertiminant
	   integer n, np
	   real(8) d, a(np,np)
	   integer indx(n)
	   integer j
	   call ludcmp(a,n,np,indx,d)
	   do j=1,n
	      d=d*a(j,j)
	   end do
	   return
	   end
       
	   
	   subroutine ludcmp(a,n,np,inx,d)
	   integer n,np,indx(n), NMAX
	   real(8) d, a(np,np), TINY
	   parameter (NMAX=500, TINY=1.0e-20)
	   integer i, imax,j,k
	   real(8) aamax, dum, summ, vv(NMAX)
	   d=1
	   do i=1,n
	      aamax=0.0d0
		  do j=1,n
		     if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j))
		  end do
		  if (aamax .eq. 0) pause 'sigular matrix in ludcmp'
		  vv(i)=1.d0/aamax
	   end do
	   do j=1,n
	      do i=1,j-1
		     summ=a(i,j)
			 do k=1, i-1
			 summ=summ-a(i,k)*a(k,j)
			 end do
			 a(i,j)=summ
              end do
              aamax=0
              do i=j,n
                 summ=a(i,j)
                 do k=1,j-1
                    summ=summ-a(i,k)*a(k,j)
                 end do
                 a(i,j)=summ
			 dum=vv(i)*abs(summ)
			 if (dum .ge. aamax) then
			    imax=i
				aamax=dum
			 end if
	       end do
		   
		   if (j .ne. imax) then
		      do k=1,n
			     dum=a(imax,k)
				 a(imax,k)=a(j,k)
				 a(j,k)=dum
			  end do
			  d=-d
			  vv(imax)=vv(j)
		    end if
			indx(j)=imax
			if(a(j,j) .eq. 0) a(j,j)=TINY
			if (j .ne. n) then
			   dum=1./a(j,j)
			   do i= j+1, n
			      a(i,j)=a(i,j)*dum
			   end do
			end if
		end do
          return
       end
