c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine lump(xx, ien, p)

	implicit none
	include "global.h"

	real* 8 xx(nsd,nn_loc)      ! nodal spatial coordinates
	real* 8  p(nsd,nn_loc)   ! element level residual
	integer ien(nen,ne)

	real* 8 eft0,det
	real* 8 sh(0:3,8)
	real* 8 xr(3,3),cf(3,3),sx(3,3)

	real* 8 x(3,8), vol_el

	integer isd,inl,node,idf,ie,iq,iflag

	iflag = 0
        do ie=1,nec
	   do isd=1,nsd
	      do inl=1,nen
		 x(isd,inl) = xx(isd,ien(inl,ie))
	      enddo
	   enddo

	   do iq=1,nquad
	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if

	      eft0 = abs(det) * wq(iq)
c	      stop
	      do isd=1,nsd
		 do inl=1,nen
		    node=ien(inl,ie)
		    p(isd,node) = p(isd,node) + sh(0,inl)*eft0
		 enddo
	      enddo	      
	   enddo

	enddo

	return
	end

