c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formd(d)
	
	implicit none
	include "global.h"
	
	real* 8  d(ndf,nnc)
	integer idf,inn

	d(:,:) = 0.0
	do inn = 1,nnc
	   do idf = 1,ndf
	      d(idf,inn) = ic(idf)
	   end do
	end do
      
	return
	end
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine formdm(d)
	
	implicit none
	include "global.h"
	
	real* 8  d(nsd,nnc)
	integer isd,inn
	
	d(:,:) = 0.0
	do inn = 1,nnc
	   do isd = 1,nsd
	      d(isd,inn) = ic(isd)
c	      d(isd,inn) = -999
	   end do
	end do
	
	return
	end
