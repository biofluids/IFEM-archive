subroutine solid_calcaccel(dt)
  use solid_variables
  implicit none

  real*8 :: dt

 !...acceleration
  solid_accel(1:nsd_solid,1:nn_solid) = (solid_vel(1:nsd_solid,1:nn_solid) - solid_prevel(1:nsd_solid,1:nn_solid))/dt     

 !...save velocity from previous timestep
  solid_prevel(1:nsd_solid,1:nn_solid) = solid_vel(1:nsd_solid,1:nn_solid)

  return
end 
      
