ccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine equal(f,d,m)
	 implicit none
       real* 8 f(m),d(m)
	 integer i,m

       do i=1,m
       d(i) = f(i) 
       enddo
       return
       end
