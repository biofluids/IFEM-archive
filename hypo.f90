!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   *.fi files are used to shorten hypo.f (keeping the overview)
!   the include command reads these files and replaces the include line
!   with the content of these files

subroutine hypo
  use global_simulation_parameter
  use global_constants
  use run_variables
  use solid_variables
  use fluid_variables
  use meshgen_fluid
  use form
  use ensight_output
  implicit none

!==============================	  
! Definition of variables
 real run_time,timer1,timer2
!============================
! Define local variables

  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"

!===============================================================
! Prepare for calculation, read in inputs or restart information

  include "hypo_prepare_fluid.fi"
  include 'hypo_write_output.fi'
  nts_start=nts_start+1
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
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
    timer1=time()
     include "hypo_fluid_solver.fi"
     timer2=time()
     write(*,*) 'total running time=', timer2-timer1
!=================================================================
! Write output file every ntsbout steps

     include "hypo_write_output.fi"

  enddo time_loop

end subroutine hypo
