c
       subroutine GJInv2(a,ainv,n,np,flag)
c
c--------------------------------------------
c
c      Inversion by 
c      Gauss-Jordan elimination
c
c--------------------------------------------
c
       implicit none
       include 'parameter.h'
c
       integer i,ipt,j,k,km1,kp1,n,np,npn,np1
       real*8 a(np,np),ainv(np,np)
       real*8 b(maxEssbc,2*maxEssbc)
       real*8 flag,eps,sa
c
       flag = 1.0
       eps  = 1.0d-20
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
       np1 = n + 1
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
	     if(dabs(b(i,k)) .gt. abs(b(ipt,k))) ipt = i
	  enddo
c
	  if (abs(b(ipt,k)) .le. eps ) then
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
