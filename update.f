c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	updated.fcm                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update(dinc, d, dg, mn, hn, hm)

      implicit none
	include "global.h"

	real* 8 d(mn,nn_loc), dg(mn,nnc)
	real* 8 dinc(mn,nn_loc)
	real* 8 hn(nnc),hm(nn_loc)
	integer i,j,mn
	
	call fclear (dinc,mn*nn_loc)
	call gather(dinc, dg, mn, hn, hm)

	do i=1,nn_loc
	do j=1,mn 
	d(j,i) = d(j,i) + dinc(j,i)
	enddo
	enddo

	return
	end
