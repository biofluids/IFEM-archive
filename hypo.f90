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
  real(8) mass_center(2)
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
!============================
! Variables for interface part
  real(8) norm_node(nsd,nn)
  integer spbcele(ne_spbc),spbcnode(nn_spbc)
  integer pbnode(2,nn_pb) !index of periodical nodes:res(1==res0(2)
  real(8) res_pb(nsd,nn_pb),res_pb_temp(nsd,nn_pb)
  real(8) res_pb_w(nsd,nn_pb),res_pb_w_temp(nsd,nn_pb)
  integer con_ele(ne_spbc),nn_con_ele
  real(8) norm_con_ele(nsd,ne_spbc)
  integer flag_contact
  real(8) vol_nn(nn) !volume for each fluid node
  real(8) curv_nn(nn)
  integer center_mapping(ne,ele_refine**nsd+1)
!  integer denote_cp_con((ele_refine**nsd)*ne_spbc)  ! for center points in contact element
!0 when I<0.5 1 for I>0.5
!  integer cp_near_inter1((ele_refine**nsd)*ne_spbc),cp_near_inter2((ele_refine**nsd)*ne_spbc) !center points near interface when 1
!  integer flag_near_inter
!============================
character(len=14) filename
character(len=14) resfilename
character(len=7) fileroot
!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
  include "hypo_declaration_interface.fi"
!============================
! Define varibales on each processor
  include "hypo_declaration_part.fi"
!===============================================================
! Prepare for calculation, read in inputs or restart information
  include "hypo_restart_file_check.fi"
  include "hypo_prepare_solid.fi"
  include "hypo_prepare_fluid.fi"
  include "hypo_prepare_interface.fi"
!===================================
! Prepare for MPI
!call readpartele(partele)
  include "hypo_prepare_part.fi"
  include "hypo_prepare_com_node.fi"
!=============================
! define the influence domain matrix
 ! integer infdomain(nn_solid)
      call mpi_barrier(mpi_comm_world,ierror)
!  write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger

!=============================
  vis_solid=  vis_liq
!  vis_solid=  1.0
  I_solid(:)=0.0
  solid_pave(:)=0.0d0
  damp_solid = 100.0

  vis_solid=0.01

  solid_vel(:,:) = 0.0
  solid_bcvel(:,:) = 0.0
  solid_bcvel_old(:,:) = 0.0
!=============================

if (edge_inflow .ne. 0) then
call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
end if
!=================================

  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)


  hsp=rkpm_scale*maxval(hg(:))
!   hsp=rkpm_scale*0.015
!  denote_cp_con(:)=0
!  cp_near_inter1(:)=0
!  cp_near_inter2(:)=0
  I_solid_inter=0.5
  max_hg=maxval(hg(:))
!  max_hg=0.015
  hg_sp=max_hg/1.0
  curv_nn(:)=0.0
  if(myid==0)write(*,*)'max_hg=',max_hg
  vol_nn(:)=0.0
  do j=1,ne
        vol_nn(ien(1:nen,j))=vol_nn(ien(1:nen,j))+hg(j)**nsd/real(nen)
  end do

  if(f_slip==1) call normal_node(norm_node,x,ien,spbcele,spbcnode)

  if(f_pb==1) then
    open(776,file='pbnode.in',status='old')
      do j=1,nn_pb
         read(776,'(I8,I8)')pbnode(1:2,j)
      end do
     close(776)
  end if
  nn_con_ele=ne_spbc
  con_ele(:)=spbcele(:)
  if(f_slip==1.and.f_con==1) call norm_contact_ele(nn_con_ele,con_ele,norm_con_ele,ien,norm_node)
  nn_center=ne-ne_spbc+(ele_refine**nsd)*ne_spbc
  call get_submesh_info(x,x_center,ien,ne_spbc,spbcele,hg,center_mapping)


!===================================
! save the orignal position of solid nodes at fluid boundary
if (node_sfcon .ne. 0 ) then
  do inode_sf=1,node_sfcon
     sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
  end do
end if

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
        solid_prevel(1:nsd_solid,1:nn_solid) = solid_vel(1:nsd_solid,1:nn_solid)
	solid_stress(:,:) = 0.0d0
	do ie=1,nn_solid
		solid_stress(1:nsd_solid,ie) = solid_pave(ie)
	end do

!-------------------------------
! correct the curr solid coor by solving the solid mon equations

call form_solidid12(id_solidbc,nsd_solid,nn_solid,ien_sbc,ne_sbc,nen_solid,ne_solid,solid_fem_con)
call solve_solid_disp(solid_coor_init,solid_coor_curr,id_solidbc,solid_fem_con,node_sbc, &
                        solid_coor_pre1,solid_vel,solid_accel,ien_sbc,solid_stress,solid_bcvel,mtype)

!--------------------------------
! Find the fluid nodes overlapping with solid domain
call search_inf_re(solid_coor_curr,x,nn,nn_solid,nsd,ne_solid,nen_solid,solid_fem_con,&
                flag_fnode,node_local,nn_local)

! Construction of the dirac deltafunctions at actual solid and fluid node positions
time=mpi_wtime()
call rkpm_nodevolume(x,nsd,nn,ien,ne,nen,ien_intlocal,ne_intlocal,dvolume,sp_radius)
call rkpm_init(solid_coor_curr,nn_solid,x,nsd,nn,dvolume,sp_radius)

time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for initiate RKPM interpolation function---', time
!==================================================================
! Solve Laplace equation get indicatior field
if (myid == 0) then
	time=mpi_wtime()
	lp_source(:,:)=0.0
        call source_laplace(x,nn,nsd,solid_coor_curr,nn_solid,solid_fem_con,ne_solid,nen_solid,&
			ien_sbc,ne_sbc,node_sbc,nn_sbc,lp_source)
	call solve_laplace(lp_source,nsd,nn,nn_solid,ien,ne,nen,x,node_sbc,nn_sbc,I_solid,flag_fnode)
! Update density and viscosity field for fluid 
	fden(:)=den_liq+density_solid*I_solid(:)
write(*,*)'density solid=',den_liq+density_solid
	fvis(:)=vis_liq+vis_solid*I_solid(:)
	time=mpi_wtime()-time
	if (myid == 0) write(*,*) '---Time for update indicator field---', time
!=================================================================
! Solid solver

!if (its .gt. 10) then

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
f_fluids(:,ie) = f_fluids(:,ie)*I_solid(ie)! * fden(ie)
end do
!f_fluids=0.0

endif
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_bcast(f_fluids(1,1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(fden(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(fvis(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(solid_pave(1),nn_solid,mpi_double_precision,0,mpi_comm_world,ierror)
      call mpi_bcast(I_solid(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
!===============================================================!
!
!                       INTERFACE PART
!
!===============================================================!

nn_inter_ini=nn_inter
x_inter_ini(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)
if(its==1) then
nn_inter_ini=0
do j=1,nn_inter
   if((x_inter(1,j).lt.3.66).or.(x_inter(1,j).gt.4.34)) then
     nn_inter_ini=nn_inter_ini+1
     x_inter_ini(1:nsd,nn_inter_ini)=x_inter(1:nsd,j)
   end if
end do
nn_inter=nn_inter_ini
x_inter(1:nsd,1:nn_inter)=x_inter_ini(1:nsd,1:nn_inter)
end if
   

call center_in_solid(I_solid,center_mapping,ien,its)
call pre_process_interface(x,x_inter,x_center,ne_intlocal,ien_intlocal,nn_local,node_local,infdomain_inter,&
	I_fluid_center,I_fluid,ien,corr_Ip,its)

maxdcurv=20.0

call pts_regen(its,x,x_inter,x_center,I_fluid,I_solid,I_fluid_center,corr_Ip,ien,center_mapping,&
                        ne_intlocal,ien_intlocal,node_local,nn_local, &
                        global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length)

call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,I_fluid)

sur_fluid(:,:)=0.0
!goto 67
call get_sur_cf(x,x_inter,x_center,I_fluid,I_solid,corr_Ip,I_fluid_center,sur_fluid,ien,flag_curv_domain)

    call solve_curvature(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,curv_nn,ien,&
                      ne_intlocal,ien_intlocal,node_local,nn_local, &
                global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
                sur_fluid,flag_curv_domain)

do j=1,nn
if(I_solid(j).lt.0.2) then
   sur_fluid(1:nsd,j)=sur_fluid(1:nsd,j)*curv_nn(j)
else
   sur_fluid(1:nsd,j)=0.0
end if
end do
if(f_slip==1) then
do j=1,nn_spbc
   sur_fluid(:,spbcnode(j))=0.0
end do
end if


time=mpi_wtime()
67 continue
!if(its==1)d(1,:)=1.5*I_fluid(:)
      include "hypo_fluid_solver.fi"

time=mpi_wtime()-time
if (myid == 0) write(*,*) '---Time for fluid solver---', time

nn_inter_ini=nn_inter

call update_center_indicator(x,x_inter,x_center,d(1:nsd,:),vol_nn,dt,I_fluid_center,corr_Ip,I_solid)
call update_x_inter(x,x_inter,x_inter_ini,vel_inter,d(1:nsd,:),dold(1:nsd,:),vol_nn,dt,I_solid)


!	if (myid == 0) then
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)


!write(*,*) '***solid velocity interpolation is commended out***'

!solid_bcvel_old(:,:) =  solid_bcvel(:,:)

    call delta_exchange(solid_bcvel_old,nn_solid,dold(1:nsd,:),nn,ndelta,dvolume,nsd, &
                                          delta_exchange_fluid_to_solid,solid_pave,dold(ndf,:))

    call delta_exchange(solid_bcvel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
					  delta_exchange_fluid_to_solid,solid_pave,d(ndf,:))


!	endif

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
if (node_sfcon .ne. 0) then
  do inode_sf=1,node_sfcon
     solid_coor_curr(1:nsd,sfcon(inode_sf))=sfxyz(1:nsd,inode_sf)
  end do
end if
222 continue
	if (myid == 0) then
     include "hypo_write_output.fi"
	endif
  enddo time_loop


end subroutine hypo
