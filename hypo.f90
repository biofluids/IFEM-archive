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
  use interface_variables
  use centermesh_variables
  use r_common, only: ninit
  use meshgen_fluid
  use meshgen_solid
  use meshgen_interface
  use form
  use ensight_output
  implicit none

!==============================	  
! Definition of variables
  integer :: klok,j,icount

  integer infdomain(nn_solid)
  real(8) mass_center(2),temp
  real(8) I_fluid_center_temp(ne)
  integer ne_inter_temp
  real(8) I_fluid_temp(nn)
  real(8) res_I
  real(8) I_Laplace(nn)
  real(8) I_Laplace_center(ne)
!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
  include "hypo_declaration_interface.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information

  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
  include "hypo_prepare_interface.fi"
!=============================
  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
! define the influence domain matrix
 ! integer infdomain(nn_solid)
  if (restart == 0) then
     include 'hypo_write_output.fi'
  else
     include "hypo_restart_read.fi"
  endif

!====get center coordinate=============
  call get_center_coor(x,x_center,ien,ien_center)
!=======================================
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


! choise of the interpolation method
if (ndelta==1) then
!=================================================================
! Construction of the dirac deltafunctions at actual solid and fluid node positions
     call delta_initialize(nn_solid,solid_coor_curr,x,ien,dvolume)

!=================================================================
! Solid solver
 write(*,*) 'starting solid solver'
    call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
                     solid_pave,solid_stress,solid_strain,solid_force_FSI)

!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
 write(*,*) 'calculating delta'
     call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
                         delta_exchange_solid_to_fluid)

!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
     include "hypo_fluid_solver.fi"

!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
    call delta_exchange(solid_vel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
					  delta_exchange_fluid_to_solid)
else if (ndelta==2) then
!=================================================================
! Construction of the FEM influence domain
!     call search_inf(solid_coor_curr,x,nn,nn_solid,nsd,ne,nen,ien,infdomain)
!      call get_center_coor(x,x_center,ien)
!=================================================================
    nn_inter_ini=nn_inter
    x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
    call search_inf_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter)
    call get_inter_ele(infdomain_inter,inter_ele,ne_inter,inter_ele_nn)

!==================================================================
! solve indicator for fluid nodes
    include 'hypo_indicator_solver.fi'
    call corr_I_fluid_center(I_fluid_center,ne_inter,x_inter,x_center,inter_ele,hg,corr_Ip)
    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
     Ic_inter=-999.0
    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)

    I_Laplace(:)=I_fluid(:)
    call points_regen(I_fluid,inter_ele_nn,x_center,I_fluid_center,inter_ele,ne_inter,Ic_inter,corr_Ip,x,x_inter,ien,nn_inter_regen,x_inter_regen,hg,infdomain_inter)

    nn_inter=nn_inter_regen
    x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter_regen)
111 format(f14.10,f14.10,f14.10)

    call search_inf_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter)
    call get_inter_ele(infdomain_inter,inter_ele,ne_inter,inter_ele_nn)
!    include 'hypo_indicator_solver.fi'
!    call reset_Icenter(I_fluid,I_fluid_center,ien,Ic_inter)
!   do icount=1,10
!    I_fluid_temp(:)=I_fluid(:)
    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
    Ic_inter=0.5
    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)
    call corr_I_fluid_center(I_fluid_center,ne_inter,x_inter,x_center,inter_ele,hg,corr_Ip)
    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
     Ic_inter=0.5
    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)




!    call get_centerpoint_Ia(x,I_fluid,x_center,I_fluid_center,corr_Ip,x_inter,hg,infdomain_inter,ne_inter,inter_ele)
!     call getnorm(I_fluid(:)-I_fluid_temp(:),I_fluid(:)-I_fluid_temp(:),nn,res_I)
!     res_I=sqrt(res_I)
!     write(*,*)'residual for I_fluid=',res_I
!   end do   
!====================================================================
! get the boundary nodes for solving normal and curvature!!!!!!!!!!!
!    call get_id_inter(id_inter,I_fluid)
!    id_inter(:)=0
!    call form_inter_bc(id_inter,rng,ien,1)
!    include 'hypo_curv_solver.fi'
!====================================================================
! start norm solver
     call get_normal_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,norm_inter,hg)
     call get_curv_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,curv_inter,hg,norm_inter)
     arc_inter(1:nn_inter)=hg(infdomain_inter(1:nn_inter))
!     arc_inter(1:nn_inter)=arc_exact(1:nn_inter)
     call get_arc_Bspline(x_inter,arc_inter,infdomain_inter,hg)
     arc_inter(1:nn_inter)=arc_inter(1:nn_inter)*3
     call get_arc_Bspline(x_inter,arc_inter,infdomain_inter,hg)
!     arc_inter(1:nn_inter)=arc_inter(1:nn_inter)*3
!     call get_arc_Bspline(x_inter,arc_inter,infdomain_inter,hg)

open(69,file='curv_inter.dat',status='unknown')
do icount=1,nn_inter
   write(69,*)icount,curv_inter(icount)
end do
close(69)
!write(*,*)'norm_inter=',norm_inter(1:2,1:nn_inter)
do icount=1,nn_inter
   temp=sqrt(norm_inter(1,icount)**2+norm_inter(2,icount)**2)
  
   norm_inter(1:2,icount)=norm_inter(1:2,icount)/temp
!   write(*,*)icount,norm_inter(1,icount),norm_inter(2,icount)
end do

   call get_surten_Bspline(x,x_inter,norm_inter,curv_inter,arc_inter,sur_fluid,hg,infdomain_inter,I_fluid,Ic_inter)
write(*,*)'ne_inter=',ne_inter
write(*,*)'# of points per ele=',nn_inter/ne_inter
write(*,*)'ne_inner=',ne_inner
write(*,*)'ne_outer=',ne_outer
!write(*,*)'Ic=',Ic_inter
! Solid solver
 write(*,*) 'starting solid solver'
!    call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
!                     solid_pave,solid_stress,solid_strain,solid_force_FSI)

!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
 write(*,*) 'calculating delta', solid_fem_con(1,3)
!     call data_exchange_FEM(solid_force_FSI,nn_solid,f_fluids,nn,dvolume,nsd,  &
!                         2,ne,nen,ne_solid,nen_solid,&
!                        solid_coor_curr,solid_fem_con,x,ien,infdomain)

!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
     include "hypo_fluid_solver.fi"
!================================================================
!interpolate fluid velocity to interfacial points
     call get_intervel_Bspline(x,x_inter,infdomain_inter,d(1:nsd,:),vel_inter,hg)

!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
! swith button should be added , right now use 1 or 2 first
!    call data_exchange_FEM(solid_vel,nn_solid,d(1:nsd,:),nn,dvolume,nsd, &
!			1,ne,nen,ne_solid,nen_solid,&
!                       solid_coor_curr,solid_fem_con,x,ien,infdomain)
end if
!=================================================================
!uPDAte solid domain
!    call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
!                     solid_vel,solid_prevel,solid_accel)

    open(unit=8406, file='masscenter.txt', status='unknown')

    mass_center(1)=sum(solid_coor_curr(1,:))/nn_solid
    mass_center(2)=sum(solid_coor_curr(2,:))/nn_solid
    write(8406,*)  mass_center(:)
!===========================================================
!=================================================================
! Write output file every ntsbout steps
!d(ndf,:)=curv_fluid(:)
!d(1:nsd,:)=norm_fluid(1:nsd,:)
!d(1:nsd,:)=sur_fluid(1:nsd,:)
!norm_inter(1:nsd,1:nn_inter)=vel_inter(1:nsd,1:nn_inter)
!write(*,*)'p_norm=',p_norm(:,:)
!I_fluid(1:nn)=curv_fluid(1:nn)
!I_fluid(1:nn)=I_Laplace(1:nn)
!I_fluid_center(:)=I_Laplace_center(:)
     include "hypo_write_output.fi"
  x_inter(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)+vel_inter(1:nsd,1:nn_inter)*dt
  enddo time_loop



end subroutine hypo
