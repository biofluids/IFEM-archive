subroutine gjinv(a,ainv,n,np,det,flag)
!
!-----------------------------------------------
!
!      Gauss-Jordan elimination
!
!      a(n,n) .........the original matrix;
!      ainv(n,n) ......the computed inverse matrix;
!      n        .......the order of the matrix;
!      np      ........the dimension of the matrix;
!      det     ........the determinant of the matrix;
!      flag          = 1: the normal value
!                    =-1: singular value
!
!      June, 1994
!
!      Written and tested by Shaofan Li
!
!-----------------------------------------------
!
  implicit none

  integer n,np
  real*8 a(np,np),ainv(np,np),b(10,20)
  real*8 det, flag


  real*8 :: eps,sa,p,aik

  integer :: i,j,k,np1,npn,ipt,kp1,km1

  flag = 1.
  eps  = 1.0e-46


  if (n .gt. np)  then
     pause 'dimension range error !' 
	 go to 999
  endif

  do i = 1,n
     do j= 1, n
	    b(i,j) = a(i,j)
	 enddo
  enddo

  np1 = n + 1
  npn = n + n

  do i = 1,n
     do j = np1, npn
        b(i,j) = 0.0
	 enddo
  enddo

  do i = 1,n
     j = n + i
	 b(i,j) = 1.0
  enddo

  do k = 1,n
     kp1 = k +1
     if (k .eq. n ) go to 500
     ipt = k

     do i = kp1, n
        if(dabs(b(i,k)) .gt. dabs(b(ipt,k))) ipt = i
     enddo

     if (dabs(b(ipt,k)) .le. eps ) then
        flag = - flag
        pause 'inversion fails'
        go to 999
     endif

     if (ipt .eq. k) go to 500

     do j = k, npn
        sa       = b(k,j)
	    b(k,j)   = b(ipt,j)
	    b(ipt,j) = sa
	 enddo

  500     do j = kp1, npn
	     b(k,j) = b(k,j)/b(k,k)
	 enddo

	 if ( k .eq. 1) go to 600
	 km1 = k -1

     do i = 1, km1
        do j = kp1, npn
           b(i,j) = b(i,j) - b(i,k)*b(k,j)
        enddo
	 enddo

	 if (k .eq. n) go to 700

  600    do i = kp1, n
	    do j = kp1, npn
	       b(i,j) = b(i,j) - b(i,k)*b(k,j)
	    enddo
	 enddo
  enddo

  
  700   do i = 1,n
     do j= 1,n
        k = j+n
        ainv(i,j) = b(i,k)
     enddo
  enddo
!
!.....Find the determinant det
!
  det = 1.0
  do k = 1, n
	 p = a(k,k)
	 do j = k, n
	    a(k,j) = a(k,j)/p
	 enddo

	 do i = 1, n
        if(i-k .lt. 0) then
	       aik = a(i,k)
	    elseif(i-k .eq. 0) then
	       go to 800
	    elseif(i-k .gt. 0) then
	       aik = a(i,k)
	    endif

	    do j = k,n
	       a(i,j) = a(i,j) - aik*a(k,j)
	    enddo
 800 enddo   
	 det = det * p
  enddo

  999  return
end
       
