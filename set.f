c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine setid(d,id,mn)

	include "global.h"

	integer id(mn,nn)
	real* 8  d(mn,nn)

	do inc=1,nn
	do idf=1,mn 
	if (id(idf,inc).eq.0) d(idf,inc) = 0.0                
	end do
	end do

	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine setd(d,f,id,mn)

	include "global.h"

	integer id(mn,nn)
	real* 8  d(mn,nn), f(mn,nn)

	do inc=1,nn
	do idf=1,mn
	if (id(idf,inc).eq.0) d(idf,inc) = f(idf,inc)         
	end do
	end do

	return
	end
