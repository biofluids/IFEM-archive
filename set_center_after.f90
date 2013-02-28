

subroutine set_center_after(I_fluid_center,ien,x_center,x_inter,corr_Ip)

  use mpi_variables
  use fluid_variables, only:ne,nen,nn,nsd
  use allocate_variables, only:den_domain, ne_den_domain
  use interface_variables, only:nn_center,maxmatrix
  include 'mpif.h'
  real(8) I_fluid_center(nn_center),I_fluid_temp1(nn_center),I_fluid_temp2(nn_center)
  integer ien(nen,ne)
  real(8) x_center(nsd,nn_center),corr_Ip(maxmatrix)
  real(8) x_inter(nsd,maxmatrix)

  integer ie,je,icount,inl,node,flag

  integer ne_loc,base,top,loc_index
  real(8) II

  I_fluid_temp1(:)=0.0
  I_fluid_temp2(:)=0.0
  if(ne_den_domain.le.ncpus) then
    if(myid+1.le.ne_den_domain) then
      ne_loc=1
    else
      ne_loc=0
    end if
  else
     base=floor(real(ne_den_domain)/real(ncpus))
     top=ne_den_domain-base*ncpus
     if(myid+1.le.top) then
        ne_loc=base+1
     else
        ne_loc=base
     end if
  end if

  do loc_index=1,ne_loc
     ie=den_domain(myid+1+(loc_index-1)*ncpus)
     if(nsd==2) then
	call get_indicator_2D(x_center(1:nsd,ie),x_inter,x_center,I_fluid_center,corr_Ip,II)
     else if(nsd==3) then
	call get_indicator_3D(x_center(1:nsd,ie),x_inter,x_center,I_fluid_center,corr_Ip,II)
     end if
     if(II.gt.1.0) II=1.0
     if(II.lt.0.0) II=0.0
     I_fluid_temp1(ie)=II
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(I_fluid_temp1(1),I_fluid_temp2(1),nn_center,mpi_double_precision, &
		mpi_sum,mpi_comm_world,ierror)

  do ie=1,nn_center
     if(I_fluid_center(ie).gt.0.5)I_fluid_center(ie)=1.0
     if(I_fluid_center(ie).lt.0.5)I_fluid_center(ie)=0.0
  end do

  do icount=1,ne_den_domain
     ie=den_domain(icount)
     I_fluid_center(ie)=I_fluid_temp2(ie)
  end do


end subroutine set_center_after
