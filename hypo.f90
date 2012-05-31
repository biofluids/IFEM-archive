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

  real(8) res_l0
  real(8) del_l0

  integer ie, inen
  real(8) var1, var2, temp
  integer nn_inter_temp
  real(8) norm_node(nsd,nn)
  integer spbcele(ne_spbc),spbcnode(nn_spbc),index_bcnode(2,ne_spbc)
  integer pbnode(2,nn_pb) !index of periodical nodes:res(1==res0(2)
  real(8) res_pb(nsd,nn_pb),res_pb_temp(nsd,nn_pb)
  real(8) res_pb_w(nsd,nn_pb),res_pb_w_temp(nsd,nn_pb)
  real(8) pre_inter(nn_solid)
  integer con_ele(ne_spbc),nn_con_ele
  real(8) norm_con_ele(nsd,ne_spbc)
  integer flag_contact
  integer ne_regen,ele_regen(100)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
!  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration ---> save space use p instead
  real(8) time
  integer nt_regen,max_regen,nn_bound
  real(8) vol_nn(nn) !volume for each fluid node
  real(8) I_fluid_temp(nn)
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

!================================================================
!find the coor for center points,read den mesh info
  call get_submesh_info(x,x_center,ien,bcnode)
!================================================================
  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
  hsp=rkpm_scale*0.025!maxval(hg(:))
  max_hg=0.025!maxval(hg(:))
if(myid==0)write(*,*)'max_hg=',max_hg
  vol_nn(:)=0.0
  do j=1,ne
	vol_nn(ien(1:nen,j))=vol_nn(ien(1:nen,j))+hg(j)**nsd/real(nen)
  end do
 call normal_node(norm_node,x,ien,spbcele,spbcnode,index_bcnode)
open(776,file='pbnode.in',status='old')
do j=1,nn_pb
  read(776,'(I8,I8)')pbnode(1:2,j)
end do
close(776)
nn_con_ele=ne_spbc
con_ele(:)=spbcele(:)
call norm_contact_ele(nn_con_ele,con_ele,norm_con_ele,ien,norm_node)
d(:,:)=0.0
!do j=1,nn
!  if(x(2,j).gt.0) d(1,j)=x(2,j)/2.0
!end do
  if (restart == 0) then
        if (myid == 0) then
         include 'hypo_write_output.fi'
        end if
  else
     include "hypo_restart_read.fi"
  endif
call mpi_barrier(mpi_comm_world,ierror)
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


else if (ndelta==2) then
call mpi_barrier(mpi_comm_world,ierror)
!goto 234
time=mpi_wtime()

call mpi_barrier(mpi_comm_world,ierror)
  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)
!  call points_removal(x_inter,nn_inter)
! find the center domain and dense mesh domain.both are the narrow band near the interface
  call find_fluid_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local)
  call find_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
				ne_intlocal,ien_intlocal)

  call get_inter_ele(infdomain_inter,ien)!, &
!if(myid==0)write(*,*)'indicator denmesh'
!  call indicator_denmesh(I_fluid,x,ien,bcnode,its,nn,ne,nen,nbc)
if(its==1) then
  I_fluid_center(:)=0.0
  I_fluid(:)=0.0
  do j=1,ne
     temp=sqrt((x_center(1,j)+3.0)**2+(x_center(2,j)-3.0)**2)
     if(temp.le.0.8) I_fluid_center(j)=1.0
  end do 
  do j=1,nn
     temp=sqrt((x(1,j)+3.0)**2+(x(2,j)-3.0)**2)
     if(temp.le.0.8)I_fluid(j)=1.0
  end do

else

  call set_center_after(I_fluid_center,I_fluid,ien)
end if
call mass_conserve(x,x_inter,x_center,hg,I_fluid_center,I_fluid,ien,corr_Ip,its)
!call points_removal(x_inter,nn_inter)
call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
				ne_intlocal,ien_intlocal)
call get_inter_ele(infdomain_inter,ien)

call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

time=mpi_wtime()-time
if(myid==0)write(*,*)'time before regen=',time

time=mpi_wtime()
if(mod(its,10)==0) then
maxdcurv=20.0


if(myid==0) then
write(*,*)'***********************************************'
write(*,*)'regenerate points for topology change'
write(*,*)'***********************************************'
end if
if(mod(its,20)==0) then
  if(nsd==3) then
  call points_regen_3D(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                        I_fluid_center,corr_Ip,hg,ien,2)
  else
    call points_regen_topo(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                            I_fluid_center,corr_Ip,hg,ien,2,nn_con_ele,con_ele,norm_con_ele)
  end if
else
    if(nsd==3) then
  call points_regen_3D(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                        I_fluid_center,corr_Ip,hg,ien,3)
  else
    call points_regen_topo(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                            I_fluid_center,corr_Ip,hg,ien,3,nn_con_ele,con_ele,norm_con_ele)
  end if
end if

  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)
   call find_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                        I_fluid)

  call set_center_after(I_fluid_center,I_fluid,ien)

  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
if(myid==0) then
write(*,*)'***********************************'
write(*,*)'regular points regeneration'
write(*,*)'************************************'
end if
  if(nsd==3) then
  call points_regen_3D(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                        I_fluid_center,corr_Ip,hg,ien,1)
  else
    call points_regen(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                            I_fluid_center,corr_Ip,hg,ien,1)
  end if
  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
!goto 222
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                        I_fluid)

  call set_center_after(I_fluid_center,I_fluid,ien)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
if(myid==0) then
write(*,*)'!!!!!!!!!!!!!end of regeneration!!!!!!!!!!!!!!!'
end if

end if
!goto 222
time=mpi_wtime()-time
if(myid==0)write(*,*)'time for regen=',time
!=======================================================================!
!******************construct contact line*******************************!
!if(its.lt.500) then

time=mpi_wtime()
 call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                               ne_intlocal,ien_intlocal)
call contactline(x_inter,x,x_center,hg,I_fluid_center,corr_Ip,ien,d(1:nsd,:),vol_nn,&
		infdomain_inter,nn_con_ele,con_ele,norm_con_ele,flag_contact,index_bcnode,I_fluid,ne_regen,ele_regen)

time=mpi_wtime()-time
if(myid==0)write(*,*)'time for contact line=',time
time=mpi_wtime()
if(flag_contact.ne.0) then
  call regulate_points(x_inter,x,nn_inter,ien,hg,ne_intlocal,ien_intlocal)
call find_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local)
 call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                               ne_intlocal,ien_intlocal)

call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                           I_fluid)
!goto 222
call set_center_after(I_fluid_center,I_fluid,ien)
call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
time=mpi_wtime()-time
if(myid==0)write(*,*)'time after contact line',time
!end if
!end if
!***********************************************************************!
if(mod(its,5)==0) then
time=mpi_wtime()
  if(nsd==3) then
  call points_regen_3D(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                        I_fluid_center,corr_Ip,hg,ien,1)
  else
    call points_regen_con(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
                            I_fluid_center,corr_Ip,hg,ien,1,ne_regen,ele_regen)
  end if
  do j=1,nn_inter_regen
     nn_inter=nn_inter+1
     x_inter(:,nn_inter)=x_inter_regen(:,j)
  end do

  call regulate_points(x_inter,x,nn_inter,ien,hg,ne_intlocal,ien_intlocal)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
time=mpi_wtime()-time
if(myid==0)write(*,*)'time for point regen after contact=',time

end if



end if
!=======================================================================!

time=mpi_wtime()
!if(mod(its,1)==0.0) then
!  call get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,&
!                        norm_inter,curv_inter,hg,dcurv)
!end if

  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
			I_fluid)
!goto 222
!  if(nsd==2) then
!  call get_arc_2D(arc_inter,norm_inter,x_inter)
!  else
!  call get_arc_3D(arc_inter,norm_inter,x_inter)
!  end if
!  call get_sur_nu(x,x_inter,hg,vol_nn,arc_inter,curv_inter,norm_inter,sur_fluid,I_fluid)
   call get_sur_cf(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,sur_fluid,hg)
!sur_fluid(:,:)=0.0
time=mpi_wtime()-time
if(myid==0)write(*,*)'time before fluid solver=',time

time=mpi_wtime()
     f_fluids(:,:)=0.0d0
234 continue
     include "hypo_fluid_solver.fi"



!goto 222
time=mpi_wtime()-time

if (myid == 0) write(*,*) 'Time for fluid solver', time

time=mpi_wtime()
! call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
 !                              ne_intlocal,ien_intlocal)
!call contactline(x_inter,x,x_center,hg,I_fluid_center,corr_Ip,ien,d(1:nsd,:),vol_nn,&
!                infdomain_inter,nn_con_ele,con_ele,norm_con_ele)
!call contactline(x_inter,x,x_center,hg,I_fluid_center,corr_Ip,ien,d(1:nsd,:),vol_nn)
nn_inter_ini=nn_inter
call update_x_inter(x,x_inter,x_inter_ini,vel_inter,d(1:nsd,:),vol_nn,dt)
time=mpi_wtime()-time
if(myid==0)write(*,*)'time after fluid solver=',time
!call get_inter_vel(x,x_inter,d(1:nsd,1:nn),vel_inter,hg,vol_nn)
!=================================================================
end if


!=================================================================
222 continue
   if (myid == 0) then
     include "hypo_write_output.fi"
	endif
call mpi_barrier(mpi_comm_world,ierror)
  enddo time_loop


end subroutine hypo
