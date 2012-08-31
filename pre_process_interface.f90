

subroutine pre_process_interface(its,x,x_center,x_inter, &
	      ien,ne_intlocal,ien_intlocal,hg,nn_local,node_local, &
	      I_fluid,I_fluid_center,corr_Ip,I_solid)

  use fluid_variables
  use interface_variables
  use mpi_variables
  use solid_bc_var
  include 'mpif.h'

  integer its   !time step

  real(8) x(nsd,nn),x_center(nsd,ne),x_inter(nsd,maxmatrix) !coor for fluid,center and interface

  integer ien(nen,ne),ne_intlocal,ien_intlocal(ne_intlocal)
  real(8) hg(ne)
  integer nn_local,node_local(nn_local)

  real(8) I_fluid(nn),I_fluid_center(ne),corr_Ip(maxmatrix),I_solid(nn)

  integer i,j,icount,jcount
  real(8) temp
  integer infdomain_inter(maxmatrix)


!related to mesh
!do i=1,nn
!   if(I_solid(i).gt.0.5)I_fluid(i)=0.5
!end do
!nn_inter_ini=nn_inter
!nn_inter=nn_inter_ini+nn_inter_ex
!x_inter(1:nsd,nn_inter_ini+1:nn_inter)=x_inter_ex(1:nsd,1:nn_inter_ex)


  call find_fluid_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local,I_solid)
  call find_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local,I_solid,ien)
!nn_inter=nn_inter_ini
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                                ne_intlocal,ien_intlocal)
  call get_inter_ele(infdomain_inter,ien)


!  call mass_conserve(x,x_inter,x_center,hg,I_fluid_center,I_fluid,ien,corr_Ip,its)
!  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
!                                ne_intlocal,ien_intlocal)
!  call get_inter_ele(infdomain_inter,ien)


!assign initial indicator and update indicate

if(its==1) then
  I_fluid_center(:)=0.0
  I_fluid(:)=0.0
!  do j=1,ne
!     temp=sqrt((x_center(1,j)+3.0)**2+(x_center(2,j)-3.0)**2)
!     if(temp.le.0.8) I_fluid_center(j)=1.0
!  end do
!  do j=1,nn
!     temp=sqrt((x(1,j)+3.0)**2+(x(2,j)-3.0)**2)
!     if(temp.le.0.8)I_fluid(j)=1.0
!  end do
   do j=1,ne
     if(x_center(1,j).gt.0.0) I_fluid_center(j)=1.0
   end do
   do j=1,nn
     if(x(1,j).gt.0.0) I_fluid(j)=1.0
   end do

do i=1,nn
!   if(I_solid(i).gt.0.4999)I_fluid(i)=1.0
    temp=sqrt((x(1,i)+0.105)**2+(x(2,i)-0.0)**2)
    if(temp.lt.0.0975)I_fluid(i)=1.0

end do

do i=1,ne
    temp=sqrt((x_center(1,i)+0.105)**2+(x_center(2,i)-0.0)**2)
    if(temp.lt.0.0975)I_fluid_center(i)=1.0
end do

!do i=1,ne
!   temp=0.0
!   do j=1,4
!      temp=temp+I_solid(ien(j,i))
!   end do
!   temp=temp/4.0
!   if(temp.ge.0.4999) I_fluid_center(i)=1.0
!end do
else

  call set_center_after(I_fluid_center,I_fluid,ien)
end if

!nn_inter=nn_inter+nn_inter_ex
!x_inter(1:nsd,nn_inter_ini+1:nn_inter)=x_inter_ex(1:nsd,1:nn_inter_ex)
  call get_correction_mf(x_inter,x_center,hg,corr_Ip,I_fluid_center)

!do i=1,ne
!   temp=0.0
!   do j=1,4
!      temp=temp+I_solid(ien(j,i))
!   end do
!   temp=temp/4.0
!   if(temp.ge.0.4999) I_fluid_center(i)=0.5
!end do

end subroutine pre_process_interface








