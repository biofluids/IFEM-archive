	subroutine nondimension  
	implicit none
	include "global.h"
	
	integer i,j

	den_gas = den_gas/ref_den
	den_liq = den_liq/ref_den
	vis_gas = vis_gas/ref_den/ref_vel/ref_lgt
	vis_liq = vis_liq/ref_den/ref_vel/ref_lgt

	do i=1,nsd
	   gravity(i) = gravity(i)*ref_lgt/ref_vel/ref_vel
	enddo
	return
	end
