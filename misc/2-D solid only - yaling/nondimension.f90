subroutine nondimension  
  use fluid_variables
  implicit none

  den_gas = den_gas/ref_den
  den_liq = den_liq/ref_den
  vis_gas = vis_gas/ref_den/ref_vel/ref_lgt
  vis_liq = vis_liq/ref_den/ref_vel/ref_lgt

  gravity(1:nsd) = gravity(1:nsd)*ref_lgt/ref_vel/ref_vel

  return
end subroutine nondimension
