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
  use denmesh_variables
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
!  real(8) I_fluid_center_temp(ne)
  integer ne_inter_temp
!  real(8) I_fluid_temp(nn)
  real(8) res_I
!  real(8) I_Laplace(nn)
!  real(8) I_Laplace_center(ne)
!============================
  integer nn_mesh
  real(8) mesh_count
  real(8) pai
  real(8) x_mesh(nsd,5000)
  integer icount_mesh
  integer regen_point_flag(maxmatrix),regen_ele(maxmatrix),nn_regen_ele
  real(8) r_l,r_s

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
  curv_bound=999.0
  nn_inter_ini=nn_inter
!  x(1:nsd,:)=0.7*x(1:nsd,:) 
!  hg(:)=hg(:)*0.7
 x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
! define the influence domain matrix
 ! integer infdomain(nn_solid)
  if (restart == 0) then
     include 'hypo_write_output.fi'
  else
     include "hypo_restart_read.fi"
  endif
  pai=3.1415926
  nn_mesh=400
  mesh_count=400.0
!  do icount_mesh=1,nn_mesh
!     x(1,icount_mesh)=-2.0*pai+(icount_mesh-1)*4.0*pai/mesh_count
!     x(2,icount_mesh)=sin(x(1,icount_mesh)+pai*1.5+0.5)
!     x(1,icount_mesh+nn_mesh)=-2.0*pai+(icount_mesh-1)*4.0*pai/mesh_count
!     x(2,icount_mesh+nn_mesh)=-sqrt(4*pai**2-x(1,icount_mesh+nn_mesh)**2)
!  end do
!  open(117,file='sinmesh.dat',status='unknown')
!  do icount_mesh=1,nn_inter
!     write(117,111)x(1,icount_mesh)/2.0/pai/4.0,x(2,icount_mesh)/2.0/pai/3.0,0.0
!      write(117,111)x_inter(1,icount_mesh),x_inter(2,icount_mesh),0.0
!  end do
!  close(117)
!stop
 






!====get center coordinate=============
  call get_submesh_info(x,x_center,ien,ien_den,x_den,bcnode_den)
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
!    nn_inter_ini=nn_inter
!    x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
!    x_inter(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)+0.03867
! search infdoman for interface points
!==================================================================================
!find interfacial points in base fluid element
    call search_inf_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter)
!find centermesh nodes in dense mesh element
    call search_inf_den(x_center,x_den,nn_den,ne,nsd,ne_den,nen_den,ien_den,infdomain_den)
!find interfacial points in dense mesh element
    call search_inf_inter(x_inter,x_den,nn_den,nn_inter,nsd,ne_den,nen_den,ien_den,infdomain_inter_den)
! get interface elements
    call get_inter_ele(infdomain_inter,inter_ele,ne_inter,inter_ele_nn,nn,ne)
    call get_inter_ele(infdomain_inter_den,inter_ele_den,ne_inter_den,inter_ele_nn_den,nn_den,ne_den)
write(*,*)'ne_inter=',ne_inter
write(*,*)'ne_inter_den=',ne_inter_den
!    write(*,*)'infdomain_inter',infdomain_inter(1:nn_inter)

!==================================================================
!if(mod(its,2)==0) then
! solve indicator for fluid nodes
    include 'hypo_indicator_solver.fi'
    nn_inter_ini=nn_inter
    x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)

    call get_indicator(x_inter,x,hg,infdomain_inter,I_fluid_center,x_center,corr_Ip,&
                       Ic_inter)
!     call get_normal(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,norm_inter,hg,Sp_sum_0d,Sp_sum_1d)
    call get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,infdomain_inter,norm_inter,curv_inter,hg)

    curv_bound=min(curv_bound*1.5,maxval(abs(curv_inter(1:nn_inter))),1.0/hg(1)*0.9)
write(*,*)'curv_bound=',curv_bound

!    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
!    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)

!     call get_normal_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,norm_inter,hg)
 !    call get_curv_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,curv_inter,hg,norm_inter)
!do icount=1,nn_inter
!write(*,*)curv_inter(icount),curv_inter_temp(icount)
!end do
!stop


!    I_Laplace(:)=I_fluid(:)
!    call points_regen(I_fluid,inter_ele_nn,x_center,I_fluid_center,inter_ele,ne_inter,Ic_inter,corr_Ip,x,x_inter,ien,nn_inter_regen,x_inter_regen,hg,infdomain_inter)
    call points_regen(x,x_inter,x_center,x_inter_regen,&
			inter_ele,ne_inter,Ic_inter,nn_inter_regen, &
			I_fluid_center,corr_Ip,hg,infdomain_inter,ien, &
			regen_point_flag,regen_ele,nn_regen_ele)
    open(448,file='xyz_ini.dat',status='unknown')
  write(448,'(a12)')'#Version 1.0'
  write(448,'(a21)')'#EnSight Point Format'

        do j=1,nn_inter
	  write(448,111)x_inter(1,j),x_inter(2,j),0.0
        end do
    close(448)
    nn_inter=nn_inter_regen
    x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter_regen)
  open(336,file='xyz_regen.dat',status='unknown')
  write(336,'(a12)')'#Version 1.0'
  write(336,'(a21)')'#EnSight Point Format'
  do j=1,nn_inter_regen
     write(336,111)x_inter_regen(1,j),x_inter_regen(2,j),0.0
  end do
  close(336)
111 format(f14.10,f14.10,f14.10)
!stop

    call search_inf_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter)
    call get_inter_ele(infdomain_inter,inter_ele,ne_inter,inter_ele_nn,nn,ne)

!    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
!    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)

!     call get_normal_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,norm_inter,hg)
!     call get_curv_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,curv_inter,hg,norm_inter)


!    call get_indicator(x_inter,x,hg,infdomain_inter,I_fluid_center,x_center,corr_Ip,&
!                       Ic_inter)
!    call get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,infdomain_inter,norm_inter,curv_inter,hg)

    call regen_in_ele(x,x_inter,x_center,Ic_inter,I_fluid_center,ien,&
			infdomain_inter,hg,norm_inter,curv_inter, &
			regen_point_flag,regen_ele,nn_regen_ele)

    curv_bound=maxval(abs(curv_inter(1:nn_inter)))

    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)


     arc_inter(1:nn_inter)=hg(infdomain_inter(1:nn_inter))
     call get_arc_Bspline(x_inter,arc_inter,infdomain_inter,hg)
     arc_inter(1:nn_inter)=arc_inter(1:nn_inter)*3
     call get_arc_Bspline(x_inter,arc_inter,infdomain_inter,hg)

open(69,file='curv_inter.dat',status='unknown')
do icount=1,nn_inter
   write(69,*)icount,curv_inter(icount)
end do
close(69)
!write(*,*)'norm_inter=',norm_inter(1:2,1:nn_inter)
!do icount=1,nn_inter
!   temp=sqrt(norm_inter(1,icount)**2+norm_inter(2,icount)**2)
  
!   norm_inter(1:2,icount)=norm_inter(1:2,icount)/temp
!   write(*,*)icount,norm_inter(1,icount),norm_inter(2,icount)
!end do

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
!     nn_inter=nn_inter_regen
!     x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)
!     call search_inf_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter)
!     call get_intervel_fem(x,x_inter,nn_inter,infdomain_inter,d(1:nsd,:),vel_inter,ien)
     call get_intervel_Bspline(x,x_inter,infdomain_inter,d(1:nsd,:),vel_inter,hg)
!     x_inter_regen(1:nsd,1:nn_inter_regen)=x_inter_regen(1:nsd,1:nn_inter_regen)+vel_inter(1:nsd,1:nn_inter_regen)*dt
!open(8888,file='updateregen.txt',status='unknown')
! do j=1,nn_inter_regen
!     write(8888,111)x_inter_regen(1,j),x_inter_regen(2,j),0.0
!  end do
!  close(8888)
!open(9998,file='origin.txt',status='unknown')
!do j=1,nn_inter
!   write(9998,111)x_inter(1,j),x_inter(2,j),0.0
!end do
!close(9998)

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

     include "hypo_write_output.fi"
  x_inter_ini(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter) 
  x_inter(1:nsd,1:nn_inter)=x_inter_ini(1:nsd,1:nn_inter)+vel_inter(1:nsd,1:nn_inter)*dt

  temp=maxval(vel_inter(1,1:nn_inter)**2+vel_inter(2,1:nn_inter)**2)**0.5
  r_l=1.0/curv_bound+temp*dt
  r_s=(1.0/curv_bound)**2/(1.0/curv_bound+temp*dt)
  curv_bound=(1.0/r_l+1.0/r_s)/2.0

write(*,*)'curv_inter_max_update=',curv_bound



  call get_intervel_Bspline(x,x_inter,infdomain_inter,d(1:nsd,:),vel_inter,hg)
  x_inter(1:nsd,1:nn_inter)=0.75*x_inter_ini(1:nsd,1:nn_inter)+0.25*x_inter(1:nsd,1:nn_inter)+0.25*vel_inter(1:nsd,1:nn_inter)*dt
  call get_intervel_Bspline(x,x_inter,infdomain_inter,d(1:nsd,:),vel_inter,hg)
  x_inter(1:nsd,1:nn_inter)=0.33333*x_inter_ini(1:nsd,1:nn_inter)+0.66667*x_inter(1:nsd,1:nn_inter)+0.66667*vel_inter(1:nsd,1:nn_inter)*dt



  enddo time_loop



end subroutine hypo
