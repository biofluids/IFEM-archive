c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	updated.fcm                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update(dinc, d, dg, mn, hn, hm)

      implicit none
	include "global.h"

	real* 8 d(mn,nn), dg(mn,nn)
	real* 8 dinc(mn,nn)
	real* 8 hn(nn),hm(nn)
	integer i,j,mn
	
	call fclear (dinc,mn*nn)
	call gather(dinc, dg, mn, hn, hm)

	do i=1,nn
	do j=1,mn 
	d(j,i) = d(j,i) + dinc(j,i)
	enddo
	enddo

	return
	end
