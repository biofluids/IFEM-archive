      subroutine permuteTF(array,Lmap,m,n)
c
c---------------------------------------------
c
c....Arrange the real array: array according
c....to the aesending order
c
c    n: the fixed dimension of array(n)
c
c    m: the first m entries that needed to be
c       exchanged
c
c---------------------------------------------
c
      implicit none
      integer m,n,i,j,mj,mj1,Lmap(n)
      real*8  array(n)
      real*8  aj,aj1
c      
      do i = 1, m - 1
	 do j = 1, m - 1
	    aj   = array(j)
	    aj1  = array(j+1)
c
	    mj   = Lmap(j)
	    mj1  = Lmap(j+1)
c
	    if (aj .gt. aj1) then
	       array(j)   = aj1
	       array(j+1) = aj
c
	       Lmap(j)   = mj1
	       Lmap(j+1) = mj
c
	    endif
	 enddo
      enddo
c      
      return
      end     
c      
c
