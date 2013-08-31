

subroutine pre_process_interface(x,x_inter,x_center,ne_intlocal,ien_intlocal,nn_local,node_local,infdomain_inter,&
		I_fluid_center,I_fluid,ien,corr_Ip,its)

  use fluid_variables, only:nn,ne,nen,nsd
  use interface_variables,only:nn_inter,nn_center,maxmatrix
  use mpi_variables

  implicit none

  real(8) x(nsd,nn),x_center(nsd,nn_center),x_inter(nsd,maxmatrix)
  integer nn_local,node_local(nn_local),ne_intlocal,ien_intlocal(ne_intlocal)
  integer infdomain_inter(maxmatrix)
  real(8) I_fluid_center(nn_center),I_fluid(nn),corr_Ip(maxmatrix)
  integer ien(nen,ne),its

  call find_fluid_domain_pa(x,x_inter,nn_local,node_local)
  call find_domain_pa(x,x_center,x_inter)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
		ne_intlocal,ien_intlocal)
if(myid==0)write(*,*)'nn_interfaere=',nn_inter
  call get_inter_ele(infdomain_inter,ien)

  if(its==1) call indicator_initialize(I_fluid,I_fluid_center,x,x_center,nn,nn_center,nsd)

  if(nsd==-999) then
    call mass_conserve(x,x_inter,x_center,I_fluid_center,I_fluid,ien,corr_Ip,its)
    call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
              ne_intlocal,ien_intlocal)
    call get_inter_ele(infdomain_inter,ien)
  end if

  call get_correction_mf(x_inter,x_center,corr_Ip,I_fluid_center)
  call set_center_after(I_fluid_center,ien,x_center,x_inter,corr_Ip)
  call get_correction_mf(x_inter,x_center,corr_Ip,I_fluid_center)

end subroutine pre_process_interface
