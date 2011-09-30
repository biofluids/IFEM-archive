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
  use r_common, only: ninit, vis_solid
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
  integer :: klok,j,inl

  integer infdomain(nn_solid)
  real(8) mass_center(2)
! Variables for fixed solid points on fluid boundary
!  integer :: node_sfcon, node_sfcon1  
!  integer sfcon_1(201) 
!  integer sfcon(402)
!  real(8) sfxyz(nsd,402)
!  integer inode_sf
! Variables for different fluid density using by implicit form  
!  integer mdata(nn_solid)
!  integer n_mdata

  real(8) res_l0
  real(8) del_l0

  integer ie, inen
  real(8) var1, var2, temp
  integer nn_inter_temp
  real(8) x_inter_temp(nsd,maxmatrix)
  real(8) R_K(nsd,maxmatrix) !used for points advection
!  real(8) max_hg
!  integer tmp_index(ne) ----> move to hypo_declare_part
! For output pressure on the solid nodes
  real(8) pre_inter(nn_solid)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
!  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration ---> save space use p instead
  real(8) time
  integer nt_regen,max_regen,nn_bound
  real(8) vol_nn(nn) !volume for each fluid node
!============================
! Define local variables
  include "hypo_declaration_solid.fi"
  include "hypo_declaration_fluid.fi"
  include "hypo_declaration_interface.fi"
!============================
! Define varibales on each processor
  include "hypo_declaration_part.fi"
  include "hypo_declaration_part_den.fi"
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
!  include "hypo_prepare_den.fi"
!=============================
! define the influence domain matrix
 ! integer infdomain(nn_solid)
      call mpi_barrier(mpi_comm_world,ierror)
  write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger

!=============================
! save the orignal position of solid nodes at fluid boundary
!  do inode_sf=1,node_sfcon
!     sfxyz(1:nsd,inode_sf)=solid_coor_init(1:nsd,sfcon(inode_sf))
!  end do
  vis_solid=vis_liq
if (edge_inflow .ne. 0) then
call edgeele(edge_inflow,rng,neface,ne,bc4el,ne_inflow)
end if

!  if (restart == 0) then
!	if (myid == 0) then
 !   	 include 'hypo_write_output.fi'
!	end if
!  else
!     include "hypo_restart_read.fi"
!  endif
!================================================================
!find the coor for center points,read den mesh info
  call get_submesh_info(x,x_center,ien,bcnode)
!================================================================
!  call search_inf_pa_den(x_center,x_den,nn_den,ne,nsd,ne_den,nen_den,&
!			ien_den,infdomain_den,ne_local_den,ien_local_den)
!  do j=1,ne
!     infdomain_den(j)=j
!  end do
!================================================================
  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
  hsp=rkpm_scale*maxval(hg(:))
  max_hg=maxval(hg(:))
if(myid==0)write(*,*)'max_hg=',max_hg
  vol_nn(:)=0.0
  do j=1,ne
	vol_nn(ien(1:nen,j))=vol_nn(ien(1:nen,j))+hg(j)**nsd/real(nen)
  end do


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
    	 write (6,'("  physical time = ",f7.3," s")') tt
    	 write (7,'("  physical time = ",f7.3," s")') tt
	end if

! choise of the interpolation method
if (ndelta==1) then

! call the communication and solid solver only on processor 0
	if (myid == 0) then

!=================================================================
! Construction of the dirac deltafunctions at actual solid and fluid node positions
     call delta_initialize(nn_solid,solid_coor_curr,x,ien,dvolume)

!=================================================================
! Solid solver
 write(*,*) 'starting solid solver'
    call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
                     solid_pave,solid_stress,solid_strain,solid_force_FSI,mtype)

!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
 write(*,*) 'calculating delta'
     call delta_exchange(solid_force_FSI,nn_solid,f_fluids,nn,ndelta,dvolume,nsd,  &
                         delta_exchange_solid_to_fluid)
	endif
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_bcast(f_fluids(1,1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
      include "hypo_fluid_solver.fi"


	if (myid == 0) then
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
    call delta_exchange(solid_vel,nn_solid,d(1:nsd,:),nn,ndelta,dvolume,nsd, &
					  delta_exchange_fluid_to_solid)
	endif

else if (ndelta==2) then
!=================================================================
! Construction of the FEM influence domain
!     call search_inf_pa(solid_coor_curr,x,nn,nn_solid,nsd,ne,nen,ien,infdomain,&
!			ne_intlocal,ien_intlocal)

!================================================================
! Merge the auxillary finf array
!    call mergefinf(infdomain,nn_solid,mdata,n_mdata)

!=================================================================
! Solid solver
 !   call solid_solver(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,  &
 !                    solid_pave,solid_stress,solid_strain,solid_force_FSI,mtype)

!=================================================================
! Set the FSI force for the solid nodes at fluid boundary to be zero
 !do inode_sf=1,node_sfcon
 !  solid_force_FSI(1:nsd,sfcon(inode_sf))=0.0
 !end do


!=================================================================
! Distribution of the solid forces to the fluid domain
!   f^fsi(t)  ->  f(t)
!if (myid ==0) then
! write(*,*) 'calculating delta'
! write(*,*) 'solid fsi force', solid_force_FSI(1,:)
!end if
!     call data_exchange_FEM(solid_force_FSI,nn_solid,f_fluids,nn,dvolume,nsd,  &
!                         2,ne,nen,ne_solid,nen_solid,&
!                        solid_coor_curr,solid_fem_con,x,ien,infdomain,d(nsd+1,:),pre_inter)
!=================================================================
! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)

!=================================================================
!=========used for interface part=================================
      call mpi_barrier(mpi_comm_world,ierror)
time=mpi_wtime()
  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)
! find the center domain and dense mesh domain.both are the narrow band near the interface
  call find_domain_pa(x_center,x_inter,ne_intlocal,ien_intlocal,&
			hg)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
				ne_intlocal,ien_intlocal)

if(myid==0)write(*,*)'begin get_inter_ele'
  call get_inter_ele(infdomain_inter,ien)!, &
if(myid==0)write(*,*)'indicator denmesh'
!  call indicator_denmesh(I_fluid,x,ien,bcnode,its,nn,ne,nen,nbc)
if(its==1) then
  I_fluid_center(:)=0.0
  do j=1,ne
!     temp=sqrt(x_center(1,j)**2+x_center(2,j)**2+(x_center(3,j)-0.3)**2)
!     if(temp.lt.0.15) I_fluid_center(j)=1.0
!     temp=sqrt(x_center(1,j)**2+x_center(2,j)**2+(x_center(3,j)-0.7)**2)
!     if(temp.lt.0.15) I_fluid_center(j)=1.0
    if( (x_center(1,j).lt.1.5) .and.(x_center(3,j).lt.0.5)) then
      I_fluid_center(j)=0.0
    end if
  end do    
else

  call set_center_after(I_fluid_center,I_fluid,ien)

end if
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

time=mpi_wtime()-time
if(myid==0)write(*,*)'time before regen=',time

time=mpi_wtime()
if(mod(its,20)==0) then
maxdcurv=10.0

if(mod(its,60)==0) then
  if(nsd==3) then
  call points_regen_3D(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                        I_fluid_center,corr_Ip,hg,ien,2)
  else
    call points_regen(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                            I_fluid_center,corr_Ip,hg,ien,2)
  end if
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)
  call find_domain_pa(x_center,x_inter,ne_intlocal,ien_intlocal,&
                        hg)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                        I_fluid)

!  call get_inter_ele(infdomain_inter,ien)
  call set_center_after(I_fluid_center,I_fluid,ien)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

end if



  if(nsd==3) then
  call points_regen_3D(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                        I_fluid_center,corr_Ip,hg,ien,1)
  else
    call points_regen(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                            I_fluid_center,corr_Ip,hg,ien,1)
  end if
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)

  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)


  if(nsd==3) then
  call points_regen_3D(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                        I_fluid_center,corr_Ip,hg,ien,1)
  else
    call points_regen(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                            I_fluid_center,corr_Ip,hg,ien,1)
  end if
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)

  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                        I_fluid)

  call set_center_after(I_fluid_center,I_fluid,ien)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
end if

time=mpi_wtime()-time
if(myid==0)write(*,*)'time for regen=',time

time=mpi_wtime()
  call get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,&
                        norm_inter,curv_inter,hg,dcurv)


  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
			I_fluid)
  if(nsd==2) then
  call get_arc_2D(arc_inter,norm_inter,x_inter)
  else
  call get_arc_3D(arc_inter,norm_inter,x_inter)
  end if
  call get_sur_nu(x,x_inter,hg,vol_nn,arc_inter,curv_inter,norm_inter,sur_fluid,I_fluid)
!sur_fluid(:,:)=0.0
time=mpi_wtime()-time
if(myid==0)write(*,*)'time before fluid solver=',time


time=mpi_wtime()
     f_fluids(:,:)=0.0d0
234 continue
     include "hypo_fluid_solver.fi"
time=mpi_wtime()-time

if (myid == 0) write(*,*) 'Time for fluid solver', time

time=mpi_wtime()
  call get_inter_vel(x,x_inter,d(1:nsd,1:nn),vel_inter,hg,vol_nn)
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
! swith button should be added , right now use 1 or 2 first
!    call data_exchange_FEM(solid_vel,nn_solid,d(1:nsd,:),nn,dvolume,nsd, &
!			1,ne,nen,ne_solid,nen_solid,&
!                       solid_coor_curr,solid_fem_con,x,ien,infdomain,d(nsd+1,:),pre_inter)
end if


!=================================================================
!uPDAte solid domain
!    call solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,  &
!                     solid_vel,solid_prevel,solid_accel)

!-----------------------------------------------------------------
! Set the solid nodes at the fluid boundary at their original position
!  do inode_sf=1,node_sfcon
!     solid_coor_curr(1:nsd,sfcon(inode_sf))=sfxyz(1:nsd,inode_sf)
!  end do
!    open(unit=8406, file='masscenter.txt', status='unknown')

!    mass_center(1)=sum(solid_coor_curr(1,:))/nn_solid
!    mass_center(2)=sum(solid_coor_curr(2,:))/nn_solid
!    write(8406,*)  mass_center(:)
!=================================================================
! Volume correction  
!   if (mod(its,10) .eq. 9) then
!   call volcorr(solid_coor_curr,nsd_solid,nen_solid,solid_fem_con,nn_solid,ne_solid, &
!                solid_coor_init)
!   write(*,*) 'volume correction applied'
!   call energy_fluid(x,d(1:nsd,:),ien)
!   end if
!=================================================================
! Write output file every ntsbout steps
!    Out put the fluid pressure at solid nodals
     solid_pave(:)=0.0d0
222 continue
   if (myid == 0) then
     include "hypo_write_output.fi"
	endif
call mpi_barrier(mpi_comm_world,ierror)
!==============update interface position==========================
  x_inter_ini(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)
  R_K(1:nsd,1:nn_inter)=vel_inter(1:nsd,1:nn_inter)*dt
  x_inter(1:nsd,1:nn_inter)=x_inter_ini(1:nsd,1:nn_inter)+1.0/6.0*R_K(1:nsd,1:nn_inter)

  x_inter_temp(1:nsd,1:nn_inter)=x_inter_ini(1:nsd,1:nn_inter)+0.5*R_K(1:nsd,1:nn_inter)
  call get_inter_vel(x,x_inter_temp,d(1:nsd,:),vel_inter,hg,vol_nn)
  R_K(1:nsd,1:nn_inter)=vel_inter(1:nsd,1:nn_inter)*dt
  x_inter(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)+1.0/3.0*R_K(1:nsd,1:nn_inter)

  x_inter_temp(1:nsd,1:nn_inter)=x_inter_ini(1:nsd,1:nn_inter)+0.5*R_K(1:nsd,1:nn_inter)
  call get_inter_vel(x,x_inter_temp,d(1:nsd,:),vel_inter,hg,vol_nn)
  R_K(1:nsd,1:nn_inter)=vel_inter(1:nsd,1:nn_inter)*dt
  x_inter(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)+1.0/3.0*R_K(1:nsd,1:nn_inter)

  x_inter_temp(1:nsd,1:nn_inter)=x_inter_ini(1:nsd,1:nn_inter)+R_K(1:nsd,1:nn_inter)
  call get_inter_vel(x,x_inter_temp,d(1:nsd,:),vel_inter,hg,vol_nn)
  R_K(1:nsd,1:nn_inter)=vel_inter(1:nsd,1:nn_inter)*dt
  x_inter(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)+1.0/6.0*R_K(1:nsd,1:nn_inter)
!==============do volume correction===============================
time=mpi_wtime()-time
if(myid==0)write(*,*)'time after fluid solver=',time  
  enddo time_loop


end subroutine hypo
