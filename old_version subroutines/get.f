	subroutine getu(d,do,u)

	implicit none
	include "global.h"
	integer i,j

	real*8 d(ndf,nn_loc),do(ndf,nn_loc),u(nsd,nn_loc)

	do i=1,nn_loc
	do j=1,nsd
	u(j,i) = alpha*d(j,i)+(1-alpha)*do(j,i)
	enddo
	enddo

	return
	end

	subroutine getfi(d,do,f)
	implicit none
	include "global.h"
	integer i,j

	real*8 d(nn_loc),do(nn_loc),f(nn_loc)

	do i=1,nn_loc
	f(i) = alpha*d(i)+(1-alpha)*do(i)
	enddo

	return
	end
