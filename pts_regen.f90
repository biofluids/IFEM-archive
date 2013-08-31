

subroutine pts_regen(its,x,x_inter,x_center,I_fluid,I_solid,I_fluid_center,corr_Ip,ien,center_mapping,&
                        ne_intlocal,ien_intlocal,node_local,nn_local, &
                        global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length)


  use fluid_variables,only:nn,ne,nen,nsd
  use interface_variables,only:nn_inter,nn_center,maxmatrix,ele_refine,I_solid_inter
  use allocate_variables,only:regen_ele,ne_regen_ele
  use mpi_variables
  implicit none

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix),I_fluid(nn),I_solid(nn)
  integer ien(nen,ne),ne_intlocal,ien_intlocal(ne_intlocal)
  integer node_local(nn_local),nn_local
  integer nn_global_com,global_com(nn_global_com)
  integer nn_local_com,local_com(nn_local_com)
  integer ad_length,send_address(ad_length,2)
  integer nn_inter_regen,its,infdomain_inter(maxmatrix)
  real(8) x_inter_regen(nsd,maxmatrix)
  integer center_mapping(ne,ele_refine**nsd+1)


if(mod(its,10)==0) then
  if(myid==0)write(*,*)'--------points regeneration begin----------'

if(mod(its,20)==0) then
  call points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,ien,2,&
                        regen_ele,ne_regen_ele,I_solid)
else
  call points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,ien,3,&
                        regen_ele,ne_regen_ele,I_solid)
end if

if(myid==0)write(*,*)'nn_regen=',nn_inter_regen

  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,ne_intlocal,ien_intlocal)

  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)

  call find_domain_pa(x,x_center,x_inter)
  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                ne_intlocal,ien_intlocal)
!  call find_inter(infdomain_inter,ien,nn_inter)
  call get_inter_ele(infdomain_inter,ien)
  call solve_indicator_laplace(x,x_inter,x_center,I_fluid,I_solid,I_fluid_center,ien,center_mapping, &
	                        ne_intlocal,ien_intlocal,node_local,nn_local, &
                        global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length)
  call get_correction_mf(x_inter,x_center,corr_Ip,I_fluid_center)
  call set_center_after(I_fluid_center,ien,x_center,x_inter,corr_Ip)
  call get_correction_mf(x_inter,x_center,corr_Ip,I_fluid_center)

   if(myid==0)write(*,*)'------------regular points regeneration--------'

  call points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,ien,1,&
                        regen_ele,ne_regen_ele,I_solid)

if(myid==0)write(*,*)'nn_regen=',nn_inter_regen

  call regulate_points(x_inter_regen,x,nn_inter_regen,ien,ne_intlocal,ien_intlocal)

  nn_inter=nn_inter_regen
  x_inter(1:nsd,1:nn_inter)=x_inter_regen(1:nsd,1:nn_inter)

  call search_inf_pa_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter,&
                ne_intlocal,ien_intlocal)

  call get_correction_mf(x_inter,x_center,corr_Ip,I_fluid_center)
if(myid==0)write(*,*)'--------end of regeneration-------------'
end if

990 continue
end subroutine pts_regen
