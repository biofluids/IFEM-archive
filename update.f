c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	updated.fcm                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update(dinc, d, dg, mn)

	implicit none
	include "global.h"

	integer i,j,mn
	real* 8 d(mn,nn),dg(mn,nn)
	real* 8 dinc(mn,nn)
	real* 8 hn(nn),hm(nn)
	
	
	call fclear (dinc,mn*nn)
	call equal(dg,dinc,mn*nn)
c	call gather(dinc, dg, mn, hn, hm)

	do i=1,nn
	   do j=1,mn 
	      d(j,i) = d(j,i) + dinc(j,i)
	   enddo
	enddo

	return
	end
