       subroutine gjinv(a,ainv,n,np,det,flag)
c
c-----------------------------------------------
c
c      Gauss-Jordan elimination
c
c      a(n,n) .........the original matrix;
c      ainv(n,n) ......the computed inverse matrix;
c      n        .......the order of the matrix;
c      np      ........the dimension of the matrix;
c      det     ........the determinant of the matrix;
c      flag          = 1: the normal value
c                    =-1: singular value
c
c      June, 1994
c
c      Written and tested by Shaofan Li
c
c-----------------------------------------------
c
       implicit real* 8 (a-h,o-z)
c
       integer n,np
       real*8 a(np,np),ainv(np,np),b(10,20)
       real*8 det, flag
c
       flag = 1.
       eps  = 1.0e-46
c
       if (n .gt. np)  then
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
	     b(i,j) = 0.0
	  enddo
       enddo
c
       do i = 1,n
	  j = n + i
	  b(i,j) = 1.0
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
c.....Find the determinant det
c
      det = 1.0
      do 850 k = 1, n
	 p = a(k,k)
	 do j = k, n
	    a(k,j) = a(k,j)/p
	 enddo
c
	 do 800 i = 1, n
            if(i-k .lt. 0) then
	       aik = a(i,k)
	    elseif(i-k .eq. 0) then
	       go to 800
	    elseif(i-k .gt. 0) then
	       aik = a(i,k)
	    endif
c
	    do j = k,n
	       a(i,j) = a(i,j) - aik*a(k,j)
	    enddo
 800    continue    
	det = det * p
 850  continue
c
  999  return
       end
       
