!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   *.fi files are used to shorten hypo.f (keeping the overview)
!   the include command reads these files and replaces the include line
!   with the content of these files

subroutine hypo
  use adaptive_meshing
  use global_simulation_parameter
  use global_constants
  use run_variables
  use delta_nonuniform
  use solid_variables
  use fluid_variables
  use r_common, only: ninit
  use meshgen_fluid
  use meshgen_gmsh
  use meshgen_solid
  use solid_fem_BC
  !use form
  use form_gmsh
  use gmres_fluid_mesh
  use gmres_fluid_ale
  use deformation_gradient
  use output_ensight
  use output_gmsh
  use restart_lib
  implicit none

!ccccccccccccccccccccccccccccc	  
! Definition of variables


  integer :: naxx1,naxx2

  integer time
  external time



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define local variables

  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare for calculation, read in inputs or restart information

  !include "hypo_restart_file_check.fi"

  if (restart == 0) then
     nts_start = 1
     tt = t_start
     klok = 0
     restart_klok = 0
  else
    call restart_file_check
  endif

  open(unit = 7,file = "echoinput.dat",status="unknown")
  
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write initial

  if (restart == 0) then
     include 'hypo_write_output.fi'
  else
     call restart_read_data(nts_start,nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface,             &
                            nn_solid,ne_solid,nen_solid,nn_solid_1,ne_solid_1,n_solid,nsurface,   &
                            nsd,ndf,nsd_solid,                                                    &
                            tt,klok,                                                              &
                            ien,ien_surf,rng,x,xref,hg,d,bid,                                     &
                            solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)
  endif

  naxx1=time()


  !call write_background_mesh(nn,x,ne,nen,ien,ne_surf,nen_surf,ien_surf,d(1:3,:))
  !write(*,*) "backgroundmesh written!!!"
  !stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  time_loop: do its = nts_start,nts !.....count from 1 (or restart-timestep) to number of timesteps

     write (6,*) ' '
     write (6,*) 'TIME STEP = ', its
     write (6,*) ' '
     write (7,*) ' '
     write (7,*) 'TIME STEP = ', its
     write (7,*) ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write restart information in binary file

     if (mod(its,restart_freq) == 0) then

        call restart_write_data(its,nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface,                   &
                                nn_solid,ne_solid,nen_solid,nn_solid_1,ne_solid_1,n_solid,nsurface,   &
                                nsd,ndf,nsd_solid,                                                    &
                                tt,klok,                                                              &
                                ien,ien_surf,rng,x,xref,hg,d,bid,                                     &
                                solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)

        restart_klok = restart_klok + 1

     endif

          
     call write_2d_scalar(nn,ne_surf,nen_surf,nsd,ien_surf,x,sqrt(d(1,:)**2 + d(2,:)**2 + d(3,:)**2),18,"Velocity Magnitude",31,"ifem_output_scalar_velocity.pos")
     call write_2d_scalar(nn,ne_surf,nen_surf,nsd,ien_surf,x,d(4,:),19,"Pressure (relative)",32,"ifem_output_scalar_pressure.pos")
	  call write_2d_vector(nn,ne_surf,nen_surf,nsd,ien_surf,x,d(1:3,1:nn),8,"Velocity",31,"ifem_output_vector_velocity.pos")
     write(*,*) "Gmsh output written"


     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output

     write (6,'("  physical time = ",f7.3," s")') tt
     write (7,'("  physical time = ",f7.3," s")') tt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALE mesh update

     include "hypo_fluid_solver_mesh.fi"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construction of the dirac deltafunctions at actual solid and fluid node positions

     !call delta_initialize(nn_solid,solid_coor_curr,nn,ne,x,ien,dvolume,cnn,ncnn,shrknode)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solid solver

     !call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
     !                  solid_pave,solid_stress,solid_strain,solid_force_FSI)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)

     !call delta_exchange(nsd,solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,  &
     !                    cnn,ncnn,shrknode,delta_exchange_solid_to_fluid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)

     include "hypo_fluid_solver.fi"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)

     !call delta_exchange(nsd,solid_vel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,   &
     !                    cnn,ncnn,shrknode,delta_exchange_fluid_to_solid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update solid domain

     !call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
     !                  solid_vel,solid_prevel,solid_accel)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write output file every ntsbout steps

     include "hypo_write_output.fi"

  enddo time_loop


 !...stops time counting and write output to screen
  naxx2=time()
  write(*,*) naxx1,naxx2,naxx2-naxx1

end subroutine hypo