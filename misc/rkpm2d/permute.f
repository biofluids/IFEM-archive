      subroutine permute(array,m,n)
c
c---------------------------------------------
c
c....Arrange the real array: array according
c....to the desending order
c
c    n: the fixed dimension of array(n)
c
c    m: the first m entries that needed to be
c       exchanged
c
c---------------------------------------------
c
      implicit none
      integer m,n,i,j
      real*8  array(n)
      real*8  aj,aj1
c      
      do i = 1, m - 1
	 do j = 1, m - 1
	    aj   = array(j)
	    aj1  = array(j+1)
	    if (aj .lt. aj1) then
	       array(j)   = aj1
	       array(j+1) = aj
	    endif
	 enddo
      enddo
c      
      return
      end     
c      
c
