c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	updated.fcm                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updatem(dinc, d, dg, mn, hn, hm)

	implicit none
	include "global.h"

	real* 8 d(mn,nn_on), dg(mn,nnc)
	real* 8 dinc(mn,nn_on)
	real* 8 hn(nnc),hm(nn_on)
	integer i,j,mn
	
	call fclear (dinc,mn*nn_on)
	call grab_all(dinc, dg, mn, hn, hm)
c	call gather(dinc,dg,mn,hn,hm)
	do i=1,nn_on
	do j=1,mn 
	d(j,i) = d(j,i) + dinc(j,i)
	enddo
	enddo

	return
	end
