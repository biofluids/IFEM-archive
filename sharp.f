	subroutine sharp(d,m,iflag)                              
	implicit none
      include "global.h"

	real* 8 d(m),t,a,c
	integer m, i, iflag

	a = delta(5)
	c = delta(6)
	if(iflag.ne.0) then
	a = 2.0
	c = 0.5
	endif

      do i=1,m  
	if(d(i).le.c) then
	t = c*((d(i)/c)**a)
	d(i) = t
	else
	t = (1.0-c)*(((1.0-d(i))/(1.0-c))**a)
	d(i) = 1.0-t
	endif
	enddo

      return
      end
