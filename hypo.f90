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
  use r_common, only: ninit,nprestress
  use meshgen_fluid
  use meshgen_solid
  use solid_fem_BC

  use ensight_output
  implicit none

!==============================	  
! Definition of variables

  integer :: klok,j

!============================
! Define local variables

  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information

  include "hypo_prepare_solid.fi"

  if (restart == 0) then
     include 'hypo_write_output.fi'
  else
  endif

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

     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output

     write (6,'("  physical time = ",f7.3," s")') tt
     write (7,'("  physical time = ",f7.3," s")') tt

!=================================================================
! Solid solver

     call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
                      solid_pave,solid_stress,solid_strain,solid_force_FSI)


!=================================================================
! Update solid domain
    call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
                      solid_vel,solid_prevel,solid_accel)

!=================================================================
! Write output file every ntsbout steps

     include "hypo_write_output.fi"

  enddo time_loop


end subroutine hypo