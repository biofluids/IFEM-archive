c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	S. Aliabadi                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine shape
	  
	  implicit none
	  include "global.h"

	  integer iq

	  if (nen.eq.4) then
		call quad3d4n(iquad, nquad, xq, wq, nsdpad, nquadpad)
	  else if (nen.eq.8) then
		call quad3d8n(iquad, nquad, xq, wq, nsdpad, nquadpad)
	  end if

      do iq=1,nquad

		if(nen.eq.4) then

c  shape function values
		  sq(0,1,iq) = xq(1,iq)
		  sq(0,2,iq) = xq(2,iq)
		  sq(0,3,iq) = xq(3,iq)
		  sq(0,4,iq) = 1 - xq(1,iq) - xq(2,iq) - xq(3,iq)
		  
		else
		  
c  shape function values
		  sq(0,1,iq) = (1 - xq(1,iq))
	1		   * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,2,iq) = (1 + xq(1,iq))
	1		   * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,3,iq) = (1 + xq(1,iq))
	1		   * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,4,iq) = (1 - xq(1,iq))
	1		   * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,5,iq) = (1 - xq(1,iq))
	1		   * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,6,iq) = (1 + xq(1,iq))
	1		   * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,7,iq) = (1 + xq(1,iq))
	1		   * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,8,iq) = (1 - xq(1,iq))
	1		   * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  
c  local first derivatives
		  sq(1,1,iq) = - (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,2,iq) = + (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,3,iq) = + (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,4,iq) = - (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,5,iq) = - (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,6,iq) = + (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,7,iq) = + (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,8,iq) = - (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  
		  sq(2,1,iq) = - (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,2,iq) = - (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,3,iq) = + (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,4,iq) = + (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,5,iq) = - (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,6,iq) = - (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,7,iq) = + (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,8,iq) = + (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
		  
		  sq(3,1,iq) = - (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,2,iq) = - (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,3,iq) = - (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,4,iq) = - (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,5,iq) = + (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,6,iq) = + (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,7,iq) = + (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,8,iq) = + (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
		  
		endif
		
	  enddo
	  
	  return
	  end
	  
