	subroutine getu(d,do,u)

	implicit none
	include "global.h"
	integer i,j

	real*8 d(ndf,nn),do(ndf,nn),u(nsd,nn)

	do i=1,nn
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

	real*8 d(nn),do(nn),f(nn)

	do i=1,nn
	f(i) = alpha*d(i)+(1-alpha)*do(i)
	enddo

	return
	end
