c
       subroutine gjinv(a,ainv,n,np,flag)
c
c--------------------------------------------
c
c      Gauss-Jordan elimination
c
c--------------------------------------------
c
       implicit double precision (a-h,o-z)
c
       dimension a(np,np),ainv(np,np),b(10,20)
c
       flag = 1.
       eps  = 1.0e-20
c
       if (n .gt. 10)  then
	  pause 'dimension range error !' 
	  go to 999
       endif
c
       do i = 1,n
	  do j= 1, n
	     b(i,j) = a(i,j)
	  enddo
       enddo
c
       np1 = n +1
       npn = n + n
c
       do i = 1,n
	  do j = np1, npn
	     b(i,j) = 0.0d0
	  enddo
       enddo
c
       do i = 1,n
	  j = n + i
	  b(i,j) = 1.0d0
       enddo
c
       do 610 k = 1,n
          kp1 = k +1
	  if (k .eq. n ) go to 500
	  ipt = k
c
	  do i = kp1, n
	     if(dabs(b(i,k)) .gt. dabs(b(ipt,k))) ipt = i
	  enddo
c
	  if (dabs(b(ipt,k)) .le. eps ) then
	      flag = - flag
	      pause 'inversion fails'
	      go to 999
	  endif
c
	  if (ipt .eq. k) go to 500
c
          do j = k, npn
	     sa       = b(k,j)
	     b(k,j)   = b(ipt,j)
	     b(ipt,j) = sa
	  enddo
c
  500     do j = kp1, npn
	     b(k,j) = b(k,j)/b(k,k)
	  enddo
c
	  if ( k .eq. 1) go to 600
	  km1 = k -1
c
	  do i = 1, km1
	     do j = kp1, npn
		b(i,j) = b(i,j) - b(i,k)*b(k,j)
	     enddo
	  enddo
c
	  if (k .eq. n) go to 700
c
  600    do i = kp1, n
	    do j = kp1, npn
	       b(i,j) = b(i,j) - b(i,k)*b(k,j)
	    enddo
	 enddo
  610  continue
c
c  
  700   do i = 1,n
	  do j= 1,n
	     k = j+n
	     ainv(i,j) = b(i,k)
	  enddo
       enddo
c
  999  return
       end
c       
