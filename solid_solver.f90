subroutine solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,  &
       solid_vel,solid_accel,solid_pave,solid_stress,solid_strain,solid_force_FSI)
  use r_common, only: xmg,density_solid,predrf
  use solid_variables, only: nn_solid,ne_solid,nen_solid,nsd_solid,nsurface
  use solid_fem_BC
  implicit none

  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con   !...connectivity for solid FEM mesh

  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_force_FSI   !...fluid structure interaction force
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_vel         !...velocity
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_accel       !...acceleration


  real(8),dimension(nn_solid)   :: solid_pave  !...averaged solid pressure (from mixed formulation -> ???)

  real(8),dimension(6,nn_solid) :: solid_stress  !...solid stress (Voigt notation)
  real(8),dimension(6,nn_solid) :: solid_strain  !...solid strain (Voigt notation)

  integer :: ipt

  write(*,*) '*** Solving Solids ***'

  select case (nsd_solid)
     case (3)    !...3D structure

     write(*,*) ' --> solving for 3-d structure'

    !...calculate internal + inertial forces + gravity/bouyancy forces   (initial configuration)  
     call r_stang(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel, &
                  solid_pave,solid_stress,solid_strain)

    !...calculate timefunction
     call r_timefun
    !...calculate nodal forces by timefunction  (initial configuration)
     call r_load
    !...apply concentrated force                (initial configuration)
     call r_nodalf

     do ipt = 1,nn_solid
        solid_force_FSI(1, ipt) = predrf(ipt           ) !+ drf(ipt           )
        solid_force_FSI(2, ipt) = predrf(ipt+  nn_solid) !+ drf(ipt+  nn_solid)
        solid_force_FSI(3, ipt) = predrf(ipt+2*nn_solid) !+ drf(ipt+2*nn_solid)
     enddo


    !...apply Dirichlet BC in current configuration (penalty stiffness)
     call solid_fem_BC_apply_essential(solid_force_FSI,solid_coor_init,solid_coor_curr)

     write(*,*) ' netto interaction force in solid (x-dir) =',sum(solid_force_FSI(1,:))

  case (0) ! 0D point

     write(*,*) ' --> solving for a point'

     solid_force_FSI(1,1)=-xmg(1)*density_solid
     solid_force_FSI(2,1)=-xmg(2)*density_solid
     solid_force_FSI(3,1)=-xmg(3)*density_solid
  end select
  
  
  return
end subroutine solid_solver