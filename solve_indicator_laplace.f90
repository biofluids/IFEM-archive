

  subroutine solve_indicator_laplace(x,x_inter,x_center,I_fluid,I_fluid_center,ien,center_mapping,&
			ne_local,ien_local,node_local,nn_local, &
			global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length)

  use fluid_variables,only:nsd,ne,nn,nen
  use interface_variables
  use allocate_variables,only:ne_inter,inter_ele,den_domain,ne_den_domain
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) I_fluid(nn),I_fluid_center(nn_center)
  integer ien(nen,ne),center_mapping(ne,ele_refine**nsd+1)

  integer ne_local,ien_local(ne_local)
  integer node_local(nn_local),nn_local
  integer nn_global_com,global_com(nn_global_com)
  integer nn_local_com,local_com(nn_local_com)
  integer ad_length,send_address(ad_length,2)
  real(8) II_center(ne),II_den_domain(ne_den_domain)
  integer node_id(nn),ele_id(ne)

  integer i,j,node,inl,isd,jsd,icount,jcount,ie
  integer flag
  real(8) p(nn),w(nn),dg(nn)

  do i=1,ne_den_domain
     ie=den_domain(i)
     II_den_domain(i)=I_fluid_center(ie)
  end do

  node_id(:)=1
  ele_id(:)=0

  do ie=1,ne
     flag=0
     do inl=1,nen
	node=ien(inl,ie)
	if(abs(I_fluid(node)-0.5).gt.0.4) then
	  node_id(node)=0
	  flag=flag+1
	  if(I_fluid(node).gt.0.5)I_fluid(node)=1.0
	  if(I_fluid(node).lt.0.5)I_fluid(node)=0.0
	end if
     end do
     if(flag.ne.nen) ele_id(ie)=1
  end do

  do i=1,ne_inter
     ie=inter_ele(i)
     do inl=1,nen
        node=ien(inl,ie)
        node_id(node)=0
        I_fluid(node)=0.5
     end do
  end do
  p(:)=0.0
  w(:)=0.0
  dg(:)=0.0
  call block_Laplace_inter(x,p,w,ien,ne_local,ien_local,I_fluid,ele_id)
  call mpi_barrier(mpi_comm_world,ierror)
  call communicate_res_ad(p,1,nn,send_address,ad_length)
  call communicate_res_ad(w,1,nn,send_address,ad_length)

  call setid_pa(p,1,nn,node_id,node_local,nn_local)
  do i=1,nn
     if(abs(w(i)).lt.1.0e-3) w(i)=1.0
  end do

  call gmres_Laplace_inter(x,w,p,dg,ien,node_id,ele_id, &
                ne_local,ien_local,node_local,nn_local, &
                global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
                I_fluid)
  call setid_pa(dg,1,nn,node_id,node_local,nn_local)
  I_fluid(:)=I_fluid(:)+dg(:)

  do ie=1,ne
     II_center(ie)=0.0
     do inl=1,nen
        node=ien(inl,ie)
        II_center(ie)=II_center(ie)+1.0/real(nen)*I_fluid(node)
     end do
     if(II_center(ie).gt.0.5) II_center(ie)=1.0
     if(II_center(ie).lt.0.5) II_center(ie)=0.0

     j=center_mapping(ie,1)
    
     do i=2,j+1
        node=center_mapping(ie,i)
        I_fluid_center(node)=II_center(ie)
     end do

  end do

  do i=1,ne_den_domain
     ie=den_domain(i)
     I_fluid_center(ie)=II_den_domain(i)
  end do

  end subroutine solve_indicator_laplace


