c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readx(xn)
	include "global.h"
	real* 8 xn(nsd,nn)

      open(unit=1,file="mxyz",status="OLD",access="DIRECT",
     $               recl=2*nsd*nn)
      read(unit=1,rec=1) xn
	close(1)

      return
      end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readien(ien)
	include "global.h"
	integer ien(nen,ne)

      open(unit=1,file="mien",status="OLD",access="DIRECT",
     $               recl=1*nen*ne)
      read(unit=1,rec=1) ien
	close(1)

	return
	end
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine readrng(rngface)
	include "global.h"
	integer rngface(neface,ne)

	open(unit=1,file="mrng",status="OLD",access="DIRECT",
     $               recl=1*neface*ne)
      read(unit=1,rec=1) rngface 
	close(1)


	mynrng = 0
	do ieface=1,neface
		do ie=1,ne
                if (rngface(ieface,ie).lt.0) rngface(ieface,ie)=0
        	mynrng = max(mynrng, rngface(ieface,ie))
		end do
	end do

	if (debug) write(0,1000) nrng
 1000	format("mrng: ",i2," boundaries")

	return
	end
