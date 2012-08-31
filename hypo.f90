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
  use solid_bc_var
  use r_common, only: ninit, vis_solid, density_solid
  use meshgen_fluid
  use meshgen_solid
  use meshgen_interface
  use form
  use ensight_output
  use mpi_variables ! call mpi variable module
  use allocate_variables
  implicit none
  include 'mpif.h'
!==============================	  
! Definition of variables
  integer :: klok,j

  integer infdomain(nn_solid)
  real(8) mass_center(nsd_solid)
  real(8) sfxyz(nsd,node_sfcon)
  integer inode_sf
! Variables for different fluid density using by implicit form  
  integer mdata(nn_solid)
  integer n_mdata
  real(8) res_l0
  real(8) del_l0

  integer ie, inen
! For output pressure on the solid nodes
  real(8) pre_inter(nn_solid)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration 
  real(8) time
  real(8) time_com
  real(8) pin_s
! solid 1st type bc temp varaibles
 real(8) omega_tmp
 real(8) theta_tmp
!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
  include "hypo_declaration_interface.fi"
!================================================
!temparily declare the variables.
real(8) norm_node(nsd,nn)  ! normal for all the fluid nodes
integer spbcele(ne_spbc),spbcnode(nn_spbc) !slip bc element and node index
integer index_bcnode(2,ne_spbc) ! for each slip be ele, store the node on the bc. for 2D,only 2 nodes stored
integer pbnode(2,nn_pb) ! index of  periodical nodes
real(8) res_pb(nsd,nn_pb),res_pb_temp(nsd,nn_pb) ! residual for periodical nodes
real(8) res_pb_w(nsd,nn_pb),res_pb_w_temp(nsd,nn_pb) ! precoditioner
integer con_ele(ne_spbc),nn_con_ele ! index of elements that form contact line. here equals to the slipbcele
real(8) norm_con_ele(nsd,ne_spbc) ! normal for contact line ele
real(8) vol_nn(nn)    ! volume for each fluid node
integer itmp
!===============================================
!============================
! Define fluid varibales on each processor
  include "hypo_declaration_part.fi"
! Define solid varibales on each processor
  include "hypo_declaration_part_solid.fi"

!===============================================================
! Prepare for calculation, read in inputs or restart information
  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
  include "hypo_prepare_interface.fi"
!===================================
! Prepare for MPI Fludi solver
  include "hypo_prepare_part.fi"
  include "hypo_prepare_com_node.fi"

! Prepare for MPISolid solver
  include "hypo_prepare_part_solid.fi"
  include "hypo_prepare_com_node_solid.fi"
!=============================
! define the influence domain matrix
 ! integer infdomain(nn_solid)
      call mpi_barrier(mpi_comm_world,ierror)
!  write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger

!=============================
!  vis_solid= - vis_liq 
  vis_solid=  1.0
  I_solid(:)=0.0
  solid_pave(:)=0.0d0
  solid_vel(:,:) = 0.0d0
  solid_accel(:,:) = 0.0d0
  damp_solid = 10000.0
  solid_bcvel(:,:) = 0.0
  solid_bcvel_old(:,:) = 0.0
  outedge=2
  fden(:)=0.0
  fvis(:)=0.0
!=============================

if (edge_inflow .ne. 0) then
call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
end if

!===================================
! save the orignal position of solid nodes at fluid boundary
if (node_sfcon .ne. 0 ) then
omega_tmp=1.0
  do inode_sf=1,node_sfcon
     sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
!     solid_vel(1,sfcon(inode_sf)) = -(solid_coor_init(2,sfcon(inode_sf)) - 1.0)/0.1*0.1*omega_tmp
!     solid_vel(2,sfcon(inode_sf)) = (solid_coor_init(1,sfcon(inode_sf)) - 4.0)/0.1*0.1*omega_tmp

  end do
end if
!======================================
! preprocessing for interface and slip bc

call get_submesh_info(x,x_center,ien,vol_nn,hg)
!nn_inter_ini=nn_inter
!x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
hsp=rkpm_scale*maxval(hg(:))
max_hg=maxval(hg(:))
call normal_node(norm_node,x,ien,spbcele,spbcnode,index_bcnode,pbnode)
nn_con_ele=ne_spbc
con_ele(:)=spbcele(:)
call norm_contact_ele(nn_con_ele,con_ele,norm_con_ele,ien,norm_node)

nn_inter_ex=nn_sbc
if(allocated(x_inter_ex))then
  deallocate(x_inter_ex)
end if
allocate(x_inter_ex(nsd,nn_inter_ex))

do itmp=1,200
 nn_inter=nn_inter+1
 x_inter(1,nn_inter)=0.0975*cos(2.0*3.14159/199.0*(itmp-1))-0.105
 x_inter(2,nn_inter)=0.0975*sin(2.0*3.14159/199.0*(itmp-1))+0.000123
end do
!========================================
  if (restart == 0) then
	if (myid == 0) then
    	 include 'hypo_write_output.fi'
	end if
  else
     include "hypo_restart_read.fi"
  endif

!=================================================================
!						 time loop	
!=================================================================
  time_loop: do its = nts_start,nts !.....count from 1 or restart-timestep to number of timesteps
      call mpi_barrier(mpi_comm_world,ierror)


	if (myid ==0) then

    	 write (6,*) ' '
    	 write (6,*) 'TIME STEP = ', its
    	 write (6,*) ' '
    	 write (7,*) ' '
    	 write (7,*) 'TIME STEP = ', its
    	 write (7,*) ' '
	
!=================================================================
! Write restart information in binary file

    	 include "hypo_restart_write.fi"
	end if


     tt = tt + dt    !....update real time
     klok = klok + 1 !....update counter for output

	if (myid ==0) then
    	 write (6,'("  physical time = ",f14.10," s")') tt
    	 write (7,'("  physical time = ",f14.10," s")') tt
	end if
! choise of the interpolation method
if (ndelta==1) then
!--------------------------------
! Update the solid coor first
        solid_coor_pre2(:,:) = solid_coor_pre1(:,:)
        solid_coor_pre1(:,:) = solid_coor_curr(:,:)
!        solid_prevel(1:nsd_solid,1:nn_solid) = solid_vel(1:nsd_solid,1:nn_solid)
	solid_stress(:,:) = 0.0d0
	do ie=1,nn_solid
!		solid_pave(ie) = 1000
!		solid_bcvel(1,ie) = 1.0 
		solid_stress(1:nsd_solid,ie) = solid_pave(ie)
	end do

if (node_sfcon .ne. 0 ) then

!  do inode_sf=1,node_sfcon
!          solid_accel(1,sfcon(inode_sf))= - 2*3.14*cos(2*3.14*tt)
!	  solid_accel(2,sfcon(inode_sf))= 0
!  end do


 ! do inode_sf=1,node_sfcon
!	solid_accel(1,sfcon(inode_sf))= - (solid_coor_curr(1,sfcon(inode_sf)) - 4.0)*(omega_tmp*omega_tmp) ! omega =10
!        solid_accel(2,sfcon(inode_sf))= - (solid_coor_curr(2,sfcon(inode_sf)) - 1.0)*(omega_tmp*omega_tmp) 
!  end do
end if
!-------------------------------
! correct the curr solid coor by solving the solid mon equations
call mpi_barrier(mpi_comm_world,ierror)

if(its==1) then
  solid_vel(1,:)=4.0
end if


if (its .gt. 3) then
if (myid == 0) write(*,*) '=== Fluid solver have converged for', its,'time steps ==='

time = mpi_wtime()
id_solidbc(:,:)=1
call form_solidid12(id_solidbc,nsd_solid,nn_solid,ien_sbc,ne_sbc,nen_solid,ne_solid,solid_fem_con)
call solve_solid_disp_pa(solid_coor_init,solid_coor_curr,id_solidbc,solid_fem_con,node_sbc, &
                        solid_coor_pre1,solid_vel,solid_accel,ien_sbc,solid_stress,solid_bcvel,mtype,&
	ne_intlocal_solid,ien_intlocal_solid,nn_local_solid,node_local_solid,send_address_solid,ad_length_solid,&
	global_com_solid,nn_global_com_solid,local_com_solid,nn_local_com_solid)



call mpi_barrier(mpi_comm_world,ierror)

time = mpi_wtime()-time
if (myid == 0)  write(*,*) 'Time for solid solver', time


end if


!if (.false.) then ! debug solid disp solver only

time=mpi_wtime()
!--------------------------------
! Find the fluid nodes overlapping with solid domain
call search_inf_re(solid_coor_curr,x,nn,nn_solid,nsd,ne_solid,nen_solid,solid_fem_con,&
                flag_fnode,node_local,nn_local)

! Construction of the dirac deltafunctions at actual solid and fluid node positions
call rkpm_nodevolume(x,nsd,nn,ien,ne,nen,ien_intlocal,ne_intlocal,dvolume,sp_radius)
call rkpm_init(solid_coor_curr,nn_solid,x,nsd,nn,dvolume,sp_radius)

time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for initiate RKPM interpolation function---', time
!==================================================================
! Solve Laplace equation get indicatior field
!if (myid == 0) then
	time=mpi_wtime()
	lp_source(:,:)=0.0
        call source_laplace(x,nn,nsd,solid_coor_curr,nn_solid,solid_fem_con,ne_solid,nen_solid,&
	ien_sbc,ne_sbc,node_sbc,nn_sbc,lp_source)
	
call mpi_barrier(mpi_comm_world,ierror)
        call solve_laplace_pa(lp_source,nsd,nn,nn_solid,ien,ne,nen,x,node_sbc,nn_sbc,I_solid,flag_fnode,&
		ne_intlocal,ien_intlocal,nn_local,node_local,send_address,ad_length,&
		global_com,nn_global_com,local_com,nn_local_com)


 !       call solve_laplace(lp_source,nsd,nn,nn_solid,ien,ne,nen,x,node_sbc,nn_sbc,I_solid,flag_fnode)


	fden(:)=den_liq+density_solid*I_solid(:)
	fvis(:)=vis_liq+vis_solid*I_solid(:)
	time=mpi_wtime()-time
!end if
if (myid == 0) write(*,*) '---Time for update indicator field---', time
!=================================================================
! Solid solver

!if (its .gt. 10) then
if (myid == 0) then
	 write(*,*) 'starting my own litte solid solver'
call res_solid(solid_coor_init,solid_coor_curr,solid_fem_con, solid_force_FSI,&
              solid_coor_pre1,solid_coor_pre2,id_solidbc,solid_accel,solid_pave,solid_bcvel,solid_bcvel_old,solid_vel)

!end if

!solid_force_FSI(:,:) =  solid_force_FSI(:,:) + solid_bcforce(:,:)
!=================================================================
! Set the FSI force for the solid nodes at fluid boundary to be zero
	if (node_sfcon .ne. 0) then
	 do inode_sf=1,node_sfcon
	   solid_force_FSI(1:nsd,sfcon(inode_sf))=0.0
	 end do
	end if 


!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
	 write(*,*) 'calculating delta'
	     call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
	                         delta_exchange_solid_to_fluid,solid_pave,d(ndf,:))

do ie = 1, nn
f_fluids(:,ie) = f_fluids(:,ie)! * fden(ie)
end do

endif
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
if(its==1) then
  do itmp=1,nn
     if(flag_fnode(itmp)==1) then
       d(1,itmp)=4.0
     end if
  end do
end if

      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_bcast(f_fluids(1,1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
!      call mpi_bcast(fden(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
!      call mpi_bcast(fvis(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
!      call mpi_bcast(I_fluid(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!================       INTERFACE PART  =========================!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!================================================================!
!1111111111111111111111111111111111111111111111111111111111111111
!PART 1: assign proximate indicator, find mesh info, get corr term
!1111111111111111111111111111111111111111111111111111111111111111
!================================================================!
!x_inter_ex(1:nsd,1:nn_inter_ex)=solid_coor_curr(1:nsd,node_sbc(1:nn_sbc))

call pre_process_interface(its,x,x_center,x_inter,&
     ien,ne_intlocal,ien_intlocal,hg,nn_local,node_local, &
     I_fluid,I_fluid_center,corr_Ip,I_solid)

!================================================================!
!2222222222222222222222222222222222222222222222222222222222222222!
!PART 2: points regeneration , topology change
!2222222222222222222222222222222222222222222222222222222222222222!
!================================================================!

call pts_regen_topo(its,x,x_center,x_inter, &
     ien,ne_intlocal,ien_intlocal,hg,nn_local,node_local, &
     I_fluid,I_fluid_center,corr_Ip,I_solid, &
     nn_con_ele,con_ele,norm_con_ele)
!================================================================!
!3333333333333333333333333333333333333333333333333333333333333333!
!PART 3: get surface tension force
!3333333333333333333333333333333333333333333333333333333333333333!
!================================================================!
call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                        I_fluid,I_solid)
call find_sur_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local,I_solid,flag_fnode)
call get_sur_cf(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,sur_fluid,hg,I_solid)

!goto 333
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!================       interface part  =========================!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
time=mpi_wtime()
      include "hypo_fluid_solver.fi"

time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for fluid solver---', time
!nn_inter=nn_inter_ini
call update_x_inter(x,x_inter,vel_inter,d(1:nsd,:),vol_nn,dt)

!	if (myid == 0) then
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)


!write(*,*) '***solid velocity interpolation is commended out***'

solid_bcvel_old(:,:) =  solid_bcvel(:,:)

    call delta_exchange(solid_bcvel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
					  delta_exchange_fluid_to_solid,solid_pave,d(ndf,:))



!	endif
333 continue
else if (ndelta==2) then
!=================================================================
! Not an option now
if (myid == 0) write(*,*) 'Not an option now as ndelta == 2 for semi-implicit FSI'

stop

end if


!=================================================================
!uPDAte solid domain
! Already updated at the beginning of the time step
	
!    call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
!                     solid_vel,solid_prevel,solid_accel)


!	include "solid_update_new.fi"


!----------------------------------------
!---------------------------------------
!end if ! debug solid disp solver only
!----------------------------------------
!---------------------------------------

!=================================================================
! Write output file every ntsbout steps
!    Out put the fluid pressure at solid nodals
!if (node_sfcon .ne. 0) then
!  do inode_sf=1,node_sfcon
!     solid_coor_curr(1:nsd,sfcon(inode_sf))=sfxyz(1:nsd,inode_sf)
!  end do
!end if
	if (myid == 0) then
     include "hypo_write_output.fi"
	endif
  enddo time_loop


end subroutine hypo
