c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine vol(xloc, ien, floc) 

      implicit none
	include "global.h"

      integer ien(nen,ne)
	real* 8 xloc(nsd,nn),floc(nn)
      real* 8 x(nsdpad,nenpad),f(nenpad),fi

      real* 8 eft0,det
      real* 8 sh(0:nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad),cf(nsdpad,nsdpad),sx(nsdpad,nsdpad)

	real* 8 e_gas,e_liq,p_gas,p_liq
	integer i,inl,ie,iq,isd

      integer ir

      p_gas =  0.0
	p_liq =  0.0

      do ie=1,ne

        do inl=1,nen
        do isd=1,nsd
        x(isd,inl) = xloc(isd,ien(inl,ie))
        enddo
	  f(inl) = floc(ien(inl,ie))
        enddo

        e_gas = 0.0
	  e_liq = 0.0
        do iq=1,nquad
	  if (nen.eq.4) then
		include "sh3d4n.h"
	  else if (nen.eq.8) then
		include "sh3d8n.h"
	  end if
        eft0 = abs(det) * wq(iq)  

	  fi = 0.0
	  do inl=1,nen
	  fi = fi+ sh(0,inl)*f(inl)
	  enddo

	  e_gas = e_gas + (1.0-fi)*eft0
	  e_liq = e_liq + fi*eft0

        enddo

	  p_gas = p_gas + e_gas
	  p_liq = p_liq + e_liq
        
        enddo

	gas=p_gas
	liq=p_liq
          
      return
      end
