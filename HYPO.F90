!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   *.fi files are used to shorten hypo.f (keeping the overview)
!   the include command reads these files and replaces the include line
!   with the content of these files

subroutine hypo
  use global_simulation_parameter
  use global_constants
  use run_variables
  use delta_nonuniform
  use solid_variables
  use fluid_variables
  use r_common, only: ninit
  use meshgen_fluid
  use meshgen_solid
  use solid_fem_BC
  use form
!  use formm
!  use gmres_fluid_mesh
!  use gmres_fluid_ale
!  use deformation_gradient
  use ensight_output
  implicit none

!==============================	  
! Definition of variables

  integer :: klok

  !integer :: naxx1,naxx2

  !integer time
  !external time


!============================
! Define local variables

  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information

  include "hypo_restart_file_check.fi"

  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"

  if (restart == 0) then
     include 'hypo_write_output.fi'
  else
     include "hypo_restart_read.fi"
  endif

  !naxx1=time()

!=================================================================
!						 time loop	
!=================================================================
  time_loop: do its = nts_start,nts !.....count from 1 or restart-timestep to number of timesteps

     write (6,*) ' '
     write (6,*) 'TIME STEP = ', its
     write (6,*) ' '
     write (7,*) ' '
     write (7,*) 'TIME STEP = ', its
     write (7,*) ' '

!=================================================================
! Write restart information in binary file

     include "hypo_restart_write.fi"

     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output

     write (6,'("  physical time = ",f7.3," s")') tt
     write (7,'("  physical time = ",f7.3," s")') tt

!=================================================================
! Construction of the dirac deltafunctions at actual solid and fluid node positions

!     call delta_initialize(nn_solid,solid_coor_curr,x,ien,dvolume)

!=================================================================
! Solid solver

!     call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
!                       solid_pave,solid_stress,solid_strain,solid_force_FSI)

!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)

!     call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
!                         delta_exchange_solid_to_fluid)

!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)

     include "hypo_fluid_solver.fi"
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)

!     call delta_exchange(solid_vel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
!						  delta_exchange_fluid_to_solid)

!=================================================================
! Update solid domain

!     call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
!                       solid_vel,solid_prevel,solid_accel)

!=================================================================
! Write output file every ntsbout steps

     include "hypo_write_output.fi"

  enddo time_loop


 !...stops time counting and write output to screen
  !naxx2=time()
  !write(*,*) naxx1,naxx2,naxx2-naxx1

end subroutine hypo