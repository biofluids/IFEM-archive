

subroutine set_center_after(I_fluid_center,I_fluid,ien)

  use mpi_variables
  use fluid_variables, only:ne,nen,nn
  use allocate_variables, only:den_domain, ne_den_domain
  include 'mpif.h'
  real(8) I_fluid_center(ne),I_fluid_temp(ne)
  real(8) I_fluid(nn)
  integer ien(nen,ne)

  integer ie,je,icount,inl,node,flag

  integer ne_loc,base,top,loc_index

  if(ne.le.ncpus) then
    if(myid+1.le.ne) then
      ne_loc=1
    else
      ne_loc=0
    end if
  else
     base=floor(real(ne)/real(ncpus))
     top=ne-base*ncpus
     if(myid+1.le.top) then
        ne_loc=base+1
     else
        ne_loc=base
     end if
  end if

  I_fluid_temp(:)=0.0
  I_fluid_center(:)=0.0


!  do ie=1,ne
  do loc_index=1,ne_loc
     ie=myid+1+(loc_index-1)*ncpus
     I_fluid_temp(ie)=0.0
     do inl=1,nen
        I_fluid_temp(ie)=I_fluid_temp(ie)+1.0/real(nen)*I_fluid(ien(inl,ie))
     end do
     if(I_fluid_temp(ie).ge.0.5) then
	I_fluid_temp(ie)=1.0
     else
	I_fluid_temp(ie)=0.0
     end if
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(I_fluid_temp(1),I_fluid_center(1),ne,mpi_double_precision, &
		mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

  do icount=1,ne_den_domain
     ie=den_domain(icount)
     I_fluid_center(ie)=0.0
     do inl=1,nen
	I_fluid_center(ie)=I_fluid_center(ie)+1.0/real(nen)*I_fluid(ien(inl,ie))
     end do
  end do

end subroutine set_center_after
