c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine diskout(d)

      include "global.h"
	real* 8 d(ndf,nn)

      open(unit=1,file="outfile",status="unknown",access="DIRECT",
     $             recl=2*ndf*nn)
      write(unit=1,rec=1) d
      close(1)

	return
	end
