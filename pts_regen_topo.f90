
subroutine pts_regen_topo(its,x,x_center,x_inter,&
           ien,ne_intlocal,ien_intlocal,hg,nn_local,node_local, &
           I_fluid,I_fluid_center,corr_Ip,I_solid, &
           nn_con_ele,con_ele,norm_con_ele)

  use fluid_variables
  use interface_variables
  use mpi_variables
  use solid_bc_var
  include 'mpif.h'

  integer its
  real(8) x(nsd,nn),x_center(nsd,ne),x_inter(nsd,maxmatrix)
  integer ien(nen,ne),ne_intlocal,ien_intlocal(ne_intlocal)
  real(8) hg(ne)
  integer nn_local,node_local(nn_local)
  real(8) I_fluid(nn),I_fluid_center(ne),corr_Ip(maxmatrix),I_solid(nn)
  integer nn_con_ele,con_ele(nn_con_ele)
  real(8) norm_con_ele(nsd,nn_con_ele)

  integer nn_inter_regen,infdomain_inter(maxmatrix)
  real(8) x_inter_regen(nsd,maxmatrix)

  integer i,j,icount,jcount
  real(8) temp

  nn_inter_regen=0
  x_inter_regen(:,:)=0.0
if(mod(its,10)==0) then
maxdcurv=20.0
!=================================================!
! first regenerate poitns to treat topology change
  if(mod(its,20)==0) then
    call points_regen_topo(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
            I_fluid_center,corr_Ip,hg,ien,2,nn_con_ele,con_ele,norm_con_ele,I_solid)
  else
    call points_regen_topo(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
            I_fluid_center,corr_Ip,hg,ien,3,nn_con_ele,con_ele,norm_con_ele,I_solid)
  end if

  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)

  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)

!  nn_inter_ini=nn_inter
!  nn_inter=nn_inter_ini+nn_inter_ex
!  x_inter(1:nsd,nn_inter_ini+1:nn_inter)=x_inter_ex(1:nsd,1:nn_inter_ex)

   call find_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local,I_solid,ien)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
! update the approximate indicator and re-calculate the correction term in order to smooth the indicator
  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                        I_fluid,I_solid)
  call set_center_after(I_fluid_center,I_fluid,ien)

!do i=1,nn
!   if(I_solid(i).gt.0.5)I_fluid(i)=0.5
!end do
!
!do i=1,ne
!   temp=0.0
!   do j=1,4
!      temp=temp+I_solid(ien(j,i))
!   end do
!   temp=temp/4.0
!   if(temp.ge.0.4999) I_fluid_center(i)=0.5
!end do


  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

!=====================================================!
! general regenerate the points to smooth the interface
200 continue
  call points_regen(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
       I_fluid_center,corr_Ip,hg,ien,1,I_solid)

  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)
  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)

!nn_inter_ini=nn_inter

!nn_inter=nn_inter_ini+nn_inter_ex
!x_inter(1:nsd,nn_inter_ini+1:nn_inter)=x_inter_ex(1:nsd,1:nn_inter_ex)


  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg,&
                        I_fluid,I_solid)

  call set_center_after(I_fluid_center,I_fluid,ien)

  
!do i=1,nn
!   if(I_solid(i).gt.0.5)I_fluid(i)=0.5
!end do
!
!do i=1,ne
!   temp=0.0
!   do j=1,4
!      temp=temp+I_solid(ien(j,i))
!   end do
!   temp=temp/4.0
!   if(temp.ge.0.4999) I_fluid_center(i)=0.5
!end do


  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)
end if

end subroutine pts_regen_topo
