	subroutine cut(d)                              
	implicit none
      include "global.h"

	real* 8 d(nn_loc),eps
	integer i

	eps = 0.01 

      do i=1,nn_loc
	if(d(i).lt.eps) d(i) = 0.0
	if(d(i).gt.1.0-eps) d(i) = 1.0
	enddo

      return
      end
