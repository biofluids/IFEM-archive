	subroutine updatex(xn,dn)
	implicit none
	include "global.h"
	real*8 xn(nsd,nnc),dn(nsd,nnc)
	integer i,j

	do i=1,nnc
		do j=1,nsd
		xn(j,i)=xn(j,i)+dn(j,i)
		enddo
	enddo
	return
	end
