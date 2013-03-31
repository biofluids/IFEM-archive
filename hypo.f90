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
  integer :: klok,j,inl,k

  integer infdomain(nn_solid)
  real(8) mass_center(2)

  real(8) res_l0
  real(8) del_l0

  integer ie, inen
  real(8) var1, var2, temp
  integer nn_inter_temp
  real(8) norm_node(nsd,nn)
  integer spbcele(ne_spbc),spbcnode(nn_spbc)
  integer pbnode(2,nn_pb) !index of periodical nodes:res(1==res0(2)
  real(8) res_pb(nsd,nn_pb),res_pb_temp(nsd,nn_pb)
  real(8) res_pb_w(nsd,nn_pb),res_pb_w_temp(nsd,nn_pb)
  real(8) pre_inter(nn_solid)
  integer con_ele(ne_spbc),nn_con_ele
  real(8) norm_con_ele(nsd,ne_spbc)
  integer flag_contact
!  integer ne_regen,ele_regen(100)
!============================
! Variables for boudary equations
  integer bc4el(ne_inflow) ! 10 is the number of nodes on edge 4
!  real(8) res_bc(nsd,nn) ! residual comming from nature B.C. integration ---> save space use p instead
  real(8) time
!  integer nt_regen,max_regen,nn_bound
  real(8) vol_nn(nn) !volume for each fluid node
!  real(8) I_fluid_temp(nn)
  real(8) curv_nn(nn)

  integer finf,index_con,contact_remove(ne_spbc),nn_contact_remove
  integer denote_cp_con(4*ne_spbc)  ! for center points in contact element
!0 when I<0.5 1 for I>0.5
  integer cp_near_inter(4*ne_spbc) !center points near interface when 1
  integer flag_near_inter
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
!  write(*,*) 'myid', myid, 'nn_local', nn_local, 'ne_local', ne_local !id for debuger

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
!  nn_center=ne-ne_spbc+8*ne_spbc
!  call get_submesh_info(x,x_center,ien,ne_spbc,spbcele)
!================================================================
  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter_ini)=x_inter(1:nsd,1:nn_inter)
  hsp=rkpm_scale*maxval(hg(:))
  denote_cp_con(:)=0
  cp_near_inter(:)=0
  max_hg=maxval(hg(:))
  hg_sp=max_hg/1.0
  curv_nn(:)=0.0
if(myid==0)write(*,*)'max_hg=',max_hg
  vol_nn(:)=0.0
  do j=1,ne
	vol_nn(ien(1:nen,j))=vol_nn(ien(1:nen,j))+hg(j)**nsd/real(nen)
  end do
 call normal_node(norm_node,x,ien,spbcele,spbcnode)!,index_bcnode)
open(776,file='pbnode.in',status='old')
do j=1,nn_pb
  read(776,'(I8,I8)')pbnode(1:2,j)
end do
close(776)
nn_con_ele=ne_spbc
con_ele(:)=spbcele(:)
call norm_contact_ele(nn_con_ele,con_ele,norm_con_ele,ien,norm_node)
  nn_center=ne-ne_spbc+(2**nsd)*ne_spbc
  call get_submesh_info(x,x_center,ien,ne_spbc,spbcele,hg)

solid_pave(:)=0.0



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
  nn_inter_ini=nn_inter
  x_inter_ini(1:nsd,1:nn_inter)=x_inter(1:nsd,1:nn_inter)
  call find_fluid_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local)
  call find_domain_pa(x,x_center,x_inter)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
				ne_intlocal,ien_intlocal)
  call get_inter_ele(infdomain_inter,ien)
if(its==1) then
  I_fluid_center(:)=0.0
  I_fluid(:)=0.0
  do j=1,nn_center
     temp=sqrt((x_center(1,j)+0.0)**2+(x_center(2,j)-0.54)**2)
     if(temp.le.0.5) I_fluid_center(j)=1.0
  end do
  do j=1,nn
     temp=sqrt((x(1,j)+0.0)**2+(x(2,j)-0.54)**2)
     if(temp.le.0.5)I_fluid(j)=1.0
  end do



!  do j=1,nn_center
!     temp=sqrt((x_center(1,j)+1.0)**2+(x_center(2,j)+0.1)**2+(x_center(3,j)-1.5)**2)
!     if(temp.lt.0.4-max_hg) I_fluid_center(j)=1.0
!     if(abs(temp-0.4).le.max_hg) I_fluid_center(j)=-0.5*(temp-0.4)/max_hg+0.5
!  end do 
!  do j=1,nn
!     temp=sqrt((x(1,j)+1.0)**2+(x(2,j)+0.1)**2+(x(3,j)-1.5)**2)
!     if(temp.le.0.4-max_hg)I_fluid(j)=1.0
!     if(abs(temp-0.4).le.max_hg)I_fluid(j)=-0.5*(temp-0.4)/max_hg+0.5
!  end do

end if
if(nsd==2) then
call mass_conserve(x,x_inter,x_center,hg,I_fluid_center,I_fluid,ien,corr_Ip,its)
call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
				ne_intlocal,ien_intlocal)
call get_inter_ele(infdomain_inter,ien)

end if

call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
call set_center_after(I_fluid_center,ien,x_center,x_inter,corr_Ip)

call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

do j=1,4*ne_spbc
   if(I_fluid_center(nn_center-4*ne_spbc+j).gt.0.5)then
      denote_cp_con(j)=1
   else
      denote_cp_con(j)=0
   end if
end do


maxdcurv=20.0
if(mod(its,1)==0) then

call wall_connect(x_inter,I_fluid_center,nn_con_ele,con_ele,infdomain_inter)
call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

if(myid==0) then
write(*,*)'***********************************************'
write(*,*)'regenerate points for topology change'
write(*,*)'***********************************************'
end if

if(mod(its,2)==0) then
  call points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,hg,ien,2,&
                        regen_ele,ne_regen_ele)
else
  call points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,hg,ien,3,&
                        regen_ele,ne_regen_ele)
end if

  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)
  call find_domain_pa(x,x_center,x_inter)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

  call set_center_after(I_fluid_center,ien,x_center,x_inter,corr_Ip)


  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
if(myid==0) then
write(*,*)'***********************************'
write(*,*)'regular points regeneration'
write(*,*)'************************************'
end if
  call points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,hg,ien,1,&
                        regen_ele,ne_regen_ele)
  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
if(myid==0) then
write(*,*)'!!!!!!!!!!!!!end of regeneration!!!!!!!!!!!!!!!'
end if

end if

!++++++++++++++++contact line for 3D+++++++++++++++++++++!
  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,&
                        I_fluid)
!goto 222
  call construct_contactline(x,x_inter,x_center,I_fluid,I_fluid_center,hg,corr_Ip,ien,&
               nn_con_ele,con_ele,norm_con_ele,rng,ne_intlocal,ien_intlocal,vol_nn,d,flag_contact)

  call mpi_barrier(mpi_comm_world,ierror)
if(flag_contact==0) goto 888
  cp_near_inter(:)=0
  do j=1,4*ne_spbc
     do k=1,ne_den_domain
        if(nn_center-4*ne_spbc+j==den_domain(k))cp_near_inter(j)=1
     end do
  end do

  call find_domain_pa(x,x_center,x_inter)


  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                               ne_intlocal,ien_intlocal)
  call find_inter(infdomain_inter,ien,nn_inter)

  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
  call set_center_after(I_fluid_center,ien,x_center,x_inter,corr_Ip)

  do j=1,4*ne_spbc
     if(cp_near_inter(j)==1) then
       flag_near_inter=0
       do k=1,ne_den_domain
          if(nn_center-4*ne_spbc+j==den_domain(k))flag_near_inter=1
          goto 567
       end do
567 continue
       if(flag_near_inter==0) then
!          if(denote_cp_con(j)==1)I_fluid_center(nn_center-4*ne_spbc+j)=0.0
!          if(denote_cp_con(j)==0)I_fluid_center(nn_center-4*ne_spbc+j)=1.0
       end if
      end if
  end do




  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
  call points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,hg,ien,1,&
                        regen_ele,ne_regen_ele)

  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

!=======================================================================!


  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,&
			I_fluid)
888 continue
sur_fluid(:,:)=0.0
   call get_sur_cf(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,sur_fluid,hg,ien,flag_curv_domain)

  call solve_curvature(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,curv_nn,hg,ien,&
                      ne_intlocal,ien_intlocal,node_local,nn_local, &
                global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
                sur_fluid,flag_curv_domain)

do j=1,nn
   sur_fluid(1:nsd,j)=sur_fluid(1:nsd,j)*curv_nn(j)
end do
!I_fluid=I_fluid_temp

do j=1,nn_spbc
   sur_fluid(:,spbcnode(j))=0.0
end do

     f_fluids(:,:)=0.0d0
234 continue
     include "hypo_fluid_solver.fi"


!goto 222


nn_inter_ini=nn_inter

call update_center_indicator(x,x_inter,x_center,d(1:nsd,:),vol_nn,dt,I_fluid_center,corr_Ip)
call update_x_inter(x,x_inter,x_inter_ini,vel_inter,d(1:nsd,:),dold(1:nsd,:),vol_nn,dt)
!call get_inter_vel(x,x_inter,d(1:nsd,1:nn),vel_inter,hg,vol_nn)
!=================================================================
end if


!=================================================================
222 continue
d(nsd+1,:)=curv_nn(:)
   if (myid == 0) then
     include "hypo_write_output.fi"
	endif
d(nsd+1,:)=0.0
call mpi_barrier(mpi_comm_world,ierror)
  enddo time_loop


end subroutine hypo
