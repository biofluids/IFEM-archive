subroutine solid_solver
  use r_common, only: xmg,density_solid,predrf,du
  use solid_variables, only: solid_force_FSI,nn_solid,nsd_solid
  use fluid_variables, only:ndf,nn
  use solid_fem_BC
  use run_variables, only:dt
  implicit none

  integer :: ipt
  real*8 :: dn(ndf,nn)

  write(*,*) '*** Solving Solids ***'

  select case (nsd_solid)
	 case (3) ! 3D structure

     write(*,*) ' --> solving for 3-d structure'

    !...calculate internal + inertial forces + gravity/bouyancy forces   (initial configuration)  
     call r_stang

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
     call solid_fem_BC_apply_essential


     write(*,*) ' netto interaction force in solid (x-dir) =',sum(solid_force_FSI(1,:))

  case (0) ! 0D point

     write(*,*) ' --> solving for a point'

     solid_force_FSI(1,1)=-xmg(1)*density_solid
     solid_force_FSI(2,1)=-xmg(2)*density_solid
     solid_force_FSI(3,1)=-xmg(3)*density_solid
  end select

  return
end subroutine solid_solver