	subroutine nondimension  
	use fluid_variables
	implicit none
	
	integer isd

	den_gas = den_gas/ref_den
	den_liq = den_liq/ref_den
	vis_gas = vis_gas/ref_den/ref_vel/ref_lgt
	vis_liq = vis_liq/ref_den/ref_vel/ref_lgt

	do isd=1,nsd
	   gravity(isd) = gravity(isd)*ref_lgt/ref_vel/ref_vel
	enddo
	return
	end
