c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	M. Behr [AHPCRC]                                                     c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_hg(x, hg)

	implicit none

	include "global.h"

	real* 8 x(nsd,nen)
	real* 8 hg

      if((etype.eq.tet).or.(etype.eq.tets))then
		hg = max(sqrt((x(1,2)-x(1,1))**2
     &			+(x(2,2)-x(2,1))**2
     &			+(x(3,2)-x(3,1))**2
     &		), sqrt((x(1,3)-x(1,2))**2
     &			+(x(2,3)-x(2,2))**2
     &			+(x(3,3)-x(3,2))**2
     &		), sqrt((x(1,1)-x(1,3))**2
     &			+(x(2,1)-x(2,3))**2
     &			+(x(3,1)-x(3,3))**2
     &		), sqrt((x(1,4)-x(1,1))**2
     &			+(x(2,4)-x(2,1))**2
     &			+(x(3,4)-x(3,1))**2
     &		), sqrt((x(1,4)-x(1,2))**2
     &			+(x(2,4)-x(2,2))**2
     &			+(x(3,4)-x(3,2))**2
     &		), sqrt((x(1,4)-x(1,3))**2
     &			+(x(2,4)-x(2,3))**2
     &			+(x(3,4)-x(3,3))**2))

      else if((etype.eq.hex).or.(etype.eq.hexs))then

			hg = 0.577350269 * max(
     &		sqrt((x(1,7)-x(1,1))**2
     &		+ (x(2,7)-x(2,1))**2
     &		+ (x(3,7)-x(3,1))**2),
     &		sqrt((x(1,5)-x(1,3))**2
     &		+ (x(2,5)-x(2,3))**2
     &		+ (x(3,5)-x(3,3))**2),
     &		sqrt((x(1,8)-x(1,2))**2
     &		+ (x(2,8)-x(2,2))**2
     &		+ (x(3,8)-x(3,2))**2),
     &		sqrt((x(1,6)-x(1,4))**2
     &		+ (x(2,6)-x(2,4))**2
     &		+ (x(3,6)-x(3,4))**2))
	else
		call error("get_hg: unsupported element type", -999, .true.)
	end if

	return
	end
