c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	updated.fcm                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine updatenode(stn,cnn,d,mdf)

	include "global.h"

	real* 8 dg(mdf,nn_on), d(mdf,nnc)
	integer stn(nnc),cnn(2,nnc)
	
	dg = 0.0
	call grab_all(dg,d,mdf)    

	do i=1,nnc
	   if(stn(i).eq.3) then
		do k=1,mdf 
		   d(k,i) = 0.0
		enddo
		   do j=1,2
	         node= cnn(j,i)
	         do k=1,mdf
		      d(k,i) = d(k,i) + dg(k,node)/2
		   enddo
		   enddo
	    endif
	enddo

	return
	end
