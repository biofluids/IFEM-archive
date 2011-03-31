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
  use r_common, only: ninit, vis_solid
  use meshgen_fluid
  use meshgen_solid
  use meshgen_interface
  use form
  use ensight_output
  use mpi_variables ! call mpi variable module
  implicit none
  include 'mpif.h'
!==============================	  
! Definition of variables
  integer :: klok,j

  integer infdomain(nn_solid)
  real(8) mass_center(2)
 
  integer part_ele_loc
  integer ne_inter_loc
  integer inter_ele_loc(ne),ne_inter_temp
  real(8) temp_ne_inter
  real(8) temp_ncpus
  integer index_regen_loc(ncpus)
  integer index_regen_loc_temp(ncpus)
  integer lower, upper,node
  real(8) x_regen_temp(nsd,maxmatrix)
  real(8) temp
  integer flag_regen
  real(8) vol_corr,arc_total
  real(8) curv_max_1, curv_max_2
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

  integer ie, inen,inl
!  integer tmp_index(ne) ----> move to hypo_declare_part
! For output pressure on the solid nodes
  real(8) pre_inter(nn_solid)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
!  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration ---> save space use p instead
  real(8) time
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
  include "hypo_prepare_den.fi"
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

  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)

  if (restart == 0) then
	if (myid == 0) then
    	 include 'hypo_write_output.fi'
	end if
  else
     include "hypo_restart_read.fi"
  endif
!====get center coordinate=============
  call get_submesh_info(x,x_center,ien,ien_den,x_den,bcnode_den)
!=======================================
!find centermesh nodes in dense mesh element
!    call search_inf_den(x_center,x_den,nn_den,ne,nsd,ne_den,nen_den,ien_den,infdomain_den)
    call search_inf_pa_den(x_center,x_den,nn_den,ne,nsd,ne_den,nen_den,ien_den,infdomain_den,ne_local_den,ien_local_den)
  ne_den_domain=ne_den
  do icount=1,ne_den
     den_domain(icount)=icount
  end do
  vol_corr=0.0
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
!find interfacial points in base fluid element
!    call search_inf_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter)
time=mpi_wtime()
    call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,ne_intlocal,ien_intlocal)
!find interfacial points in dense mesh element
!    call search_inf_inter(x_inter,x_den,nn_den,nn_inter,nsd,ne_den,nen_den,ien_den,infdomain_inter_den)
    call search_inf_pa_inter(x_inter,x_den,nn_den,nn_inter,nsd,ne_den,nen_den,ien_den,infdomain_inter_den,ne_local_den,ien_local_den)
! get interface elements
    call get_inter_ele(infdomain_inter,inter_ele,ne_inter,inter_ele_nn,nn,ne)
    call get_inter_ele(infdomain_inter_den,inter_ele_den,ne_inter_den,inter_ele_nn_den,nn_den,ne_den)
if(myid==0) then
write(*,*)'ne_inter=',ne_inter
write(*,*)'ne_inter_den=',ne_inter_den
end if
time=mpi_wtime()-time
if(myid==0) write(*,*)'time for search',time
time=mpi_wtime()
    include 'hypo_indicator_solver.fi'
    nn_inter_ini=nn_inter
    x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
    Ic_inter=0.5
!I_fluid_center(3589)=0.0
time=mpi_wtime()-time
if(myid==0)write(*,*)'time for indicator',time
time=mpi_wtime()
      call set_element_index(ne_inter_temp,ne_den_domain,den_domain,infdomain_den,ne_inter,ne_inner,ne_outer,inter_ele,inner_ele,outer_ele,I_fluid_center,I_fluid_den,ien_den,ien,inter_ele_den,ne_inter_den,ne_outer_den,ne_inner_den,outer_ele_den,inner_ele_den)
    Ic_inter=0.5
!    call get_Ia_ls(x_inter,x_center,infdomain_inter,ne_inter_temp,inter_ele,hg,I_fluid_center,Ic_inter)
    call get_Ia_distance(ne_inter,inter_ele,x_inter,x_center,hg,I_fluid_center)
if(myid==0) write(*,*)'===========begin volume correction=============='
    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)
    call get_normal_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,norm_inter,hg)
    arc_inter(1:nn_inter)=hg(infdomain_inter(1:nn_inter))
    call get_arc_Bspline(arc_total,x_inter,arc_inter,infdomain_inter,hg)
    do icount=1,nn_inter
       temp=sqrt(norm_inter(1,icount)**2+norm_inter(2,icount)**2)
       x_inter(1:nsd,icount)=x_inter(1:nsd,icount)-vol_corr/arc_total*norm_inter(1:nsd,icount)/temp
    end do
    call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,ne_intlocal,ien_intlocal)
300 continue
    flag_regen=0
200 continue
    flag_regen=flag_regen+1


    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)
!==============
     call get_normal_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,norm_inter,hg)

     call get_curv_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,curv_inter,hg,norm_inter)
     if(flag_regen==1) then
	curv_max_1=maxval(abs(curv_inter(1:nn_inter)))
     else if(flag_regen==2) then
        curv_max_2=maxval(abs(curv_inter(1:nn_inter)))
     end if


time=mpi_wtime()-time
if(myid==0)write(*,*)'time for normal curvature',time
time=mpi_wtime()
index_regen_loc(:)=0
index_regen_loc_temp(:)=0
if(ne_inter.gt.ncpus) then
   temp_ne_inter=dble(ne_inter)
   temp_ncpus=dble(ncpus)
   part_ele_loc=floor(temp_ne_inter/temp_ncpus)
   if(myid.ne.(ncpus-1)) then
     ne_inter_loc=part_ele_loc
   else
     ne_inter_loc=ne_inter-(ncpus-1)*part_ele_loc
   end if
!   index_regen_loc(:)=0
!   index_regen_loc_temp(:)=0
else
  if(myid+1.le.ne_inter) then
    part_ele_loc=1
    ne_inter_loc=1
  else
    part_ele_loc=0
    ne_inter_loc=0
  end if
end if
   call mpi_barrier(mpi_comm_world,ierror)

    call points_regen(flag_regen,part_ele_loc,I_fluid,inter_ele_nn,x_center,I_fluid_center,inter_ele,ne_inter_loc,Ic_inter,corr_Ip,x,x_inter,ien,nn_inter_regen_loc,x_inter_regen_loc,hg,infdomain_inter)

   index_regen_loc_temp(myid+1)=nn_inter_regen_loc
if(nn_inter_regen_loc.gt.0) then
! write(*,*)'myid=',myid,'nn_regen=',nn_inter_regen_loc
end if
   call mpi_barrier(mpi_comm_world,ierror)
   call mpi_reduce(index_regen_loc_temp(1),index_regen_loc(1),ncpus,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
   call mpi_bcast(index_regen_loc(1),ncpus,mpi_double_precision,0,mpi_comm_world,ierror)
   x_inter_regen(:,:)=0.0
   x_regen_temp(:,:)=0.0
   call mpi_barrier(mpi_comm_world,ierror)
   nn_inter_regen=0
   do icount=1,ncpus
      nn_inter_regen=nn_inter_regen+index_regen_loc(icount)
   end do
   lower=0
   upper=0
   if(myid==0) then
     lower=0
     upper=index_regen_loc(1)
   else
     do icount=1,myid
	lower=lower+index_regen_loc(icount)
     end do
     upper=lower+index_regen_loc(myid+1)
     
   end if
   do icount=1,index_regen_loc(myid+1) 
      x_regen_temp(1:nsd,lower+icount)=x_inter_regen_loc(1:nsd,icount) 
   end do

   call mpi_barrier(mpi_comm_world,ierror)
   call mpi_reduce(x_regen_temp(1,1),x_inter_regen(1,1),nsd*maxmatrix,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
   call mpi_bcast(x_inter_regen(1,1),nsd*maxmatrix,mpi_double_precision,0,mpi_comm_world,ierror)
   call mpi_barrier(mpi_comm_world,ierror)

if(myid==0) then
  open(337,file='xyz_ini.dat',status='unknown')
  write(337,'(a12)')'#Version 1.0'
  write(337,'(a21)')'#EnSight Point Format'
  do j=1,nn_inter
     write(337,111)x_inter(1,j),x_inter(2,j),0.0
  end do  
  close(337)
end if




    nn_inter=nn_inter_regen
    x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter_regen)
if(myid==0) then
  write(*,*)'nn_inter_regen=',nn_inter_regen
  open(336,file='xyz_regen.dat',status='unknown')
  write(336,'(a12)')'#Version 1.0'
  write(336,'(a21)')'#EnSight Point Format'
  do j=1,nn_inter_regen  
     write(336,111)x_inter_regen(1,j),x_inter_regen(2,j),0.0
  end do  
  close(336)
end if
111 format(f14.10,f14.10,f14.10)
time=mpi_wtime()-time
if(myid==0)write(*,*)'time for regeneration',time
    call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,ne_intlocal,ien_intlocal)  
if(flag_regen==1) then
goto 200
end if
if(curv_max_2.gt.2.0*curv_max_1) then
goto 300
end if

    call get_inter_ele(infdomain_inter,inter_ele,ne_inter,inter_ele_nn,nn,ne)

!    call get_Ia_ls(x_inter,x_center,infdomain_inter,ne_inter,inter_ele,hg,I_fluid_center,Ic_inter)

    call get_interpoint_Ia(x_inter,x,I_inter,I_fluid,hg,infdomain_inter,I_fluid_center,x_center)
    call correct_Ip(x_inter,x,I_inter,I_fluid,hg,Ic_inter,infdomain_inter,corr_Ip,2)
    call get_normal_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,norm_inter,hg)  
    call get_curv_Bspline(x_inter,x_center,I_fluid_center,infdomain_inter,corr_Ip,curv_inter,hg,norm_inter)
    arc_inter(1:nn_inter)=hg(infdomain_inter(1:nn_inter))
!     arc_inter(1:nn_inter)=arc_exact(1:nn_inter)
    call get_arc_Bspline(arc_total,x_inter,arc_inter,infdomain_inter,hg)
    do icount=1,nn_inter
       temp=sqrt(norm_inter(1,icount)**2+norm_inter(2,icount)**2)
       norm_inter(1:2,icount)=norm_inter(1:2,icount)/temp
    end do
   call get_surten_Bspline(x,x_inter,norm_inter,curv_inter,arc_inter,sur_fluid,hg,infdomain_inter,I_fluid,Ic_inter)



! FEM Navier-Stokes Solver (GMRES) - calculates v(t+dt),p(t+dt)
      call mpi_barrier(mpi_comm_world,ierror)
time=mpi_wtime()
     f_fluids(:,:)=0.0d0
     include "hypo_fluid_solver.fi"
time=mpi_wtime()-time
if (myid == 0) write(*,*) 'Time for fluid solver', time

!if(myid==0) then
open(123,file='fluid_vel.dat',status='unknown')
do j=1,nn
   read(123,121)d(1,j),d(2,j)
end do
close(123)
!end if
121 format(f14.10 f14.10)
d(1,1:nn)=-1.0*d(1,1:nn)
d(2,1:nn)=-1.0*d(2,1:nn)
!=================================================================
! Interpolation fluid velocity -> immersed material points
!     v^f(t+dt)  ->  v^s(t+dt)
! swith button should be added , right now use 1 or 2 first
!    call data_exchange_FEM(solid_vel,nn_solid,d(1:nsd,:),nn,dvolume,nsd, &
!			1,ne,nen,ne_solid,nen_solid,&
!                       solid_coor_curr,solid_fem_con,x,ien,infdomain,d(nsd+1,:),pre_inter)
end if
     call get_intervel_Bspline(x,x_inter,infdomain_inter,d(1:nsd,:),vel_inter,hg)

 !  call get_intervel_fem(x,x_inter,nn_inter,infdomain_inter,d(1:nsd,:),vel_inter,ien)
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
	if (myid == 0) then
     include "hypo_write_output.fi"
	endif
!===============================
!  nn_inter=nn_inter_ini
!  x_inter(2,1:nn_inter)=x_inter_ini(2,1:nn_inter)+0.0001234     
!  x_inter(1,1:nn_inter)=x_inter_ini(1,1:nn_inter)
  x_inter_ini(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)
  x_inter(1:nsd,1:nn_inter)=x_inter_ini(1:nsd,1:nn_inter)+vel_inter(1:nsd,1:nn_inter)*dt
  call get_intervel_Bspline(x,x_inter,infdomain_inter,d(1:nsd,:),vel_inter,hg)

!     call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,ne_intlocal,ien_intlocal)   
!     call get_intervel_fem(x,x_inter,nn_inter,infdomain_inter,d(1:nsd,:),vel_inter,ien)

  x_inter(1:nsd,1:nn_inter)=0.75*x_inter_ini(1:nsd,1:nn_inter)+0.25*x_inter(1:nsd,1:nn_inter)+0.25*vel_inter(1:nsd,1:nn_inter)*dt
  call get_intervel_Bspline(x,x_inter,infdomain_inter,d(1:nsd,:),vel_inter,hg)

!     call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,ne_intlocal,ien_intlocal)   
!     call get_intervel_fem(x,x_inter,nn_inter,infdomain_inter,d(1:nsd,:),vel_inter,ien)

  x_inter(1:nsd,1:nn_inter)=0.33333*x_inter_ini(1:nsd,1:nn_inter)+0.66667*x_inter(1:nsd,1:nn_inter)+0.66667*vel_inter(1:nsd,1:nn_inter)*dt
!==========used for volume correction==============================
   vol_corr=0.0
   do icount=1,nn_inter
!         vol_corr=vol_corr+arc_inter(icount)*dt*(vel_inter(1,icount)*norm_inter(1,icount)+vel_inter(2,icount)*norm_inter(2,icount))
       vol_corr=vol_corr+arc_inter(icount)*((x_inter(1,icount)-x_inter_ini(1,icount))*norm_inter(1,icount)+&
					    (x_inter(2,icount)-x_inter_ini(2,icount))*norm_inter(2,icount))
   end do
   if(myid==0) then
      write(*,*)'vol_corr=',vol_corr
   end if

  enddo time_loop


end subroutine hypo
