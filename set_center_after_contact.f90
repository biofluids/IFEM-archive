

subroutine set_center_after_contact(I_fluid_center,ien,contact_remove,nn_contact_remove,spbcele)

  use mpi_variables
  use fluid_variables, only:ne,nen,nn,ne_spbc
  use allocate_variables, only:den_domain, ne_den_domain
  use interface_variables, only:nn_center
  include 'mpif.h'
  real(8) I_fluid_center(nn_center),I_fluid_temp(nn_center)
!  real(8) I_fluid(nn)
  integer ien(nen,ne)
  integer contact_remove(ne_spbc),nn_contact_remove
  integer spbcele(ne_spbc)

  integer ie,je,icount,inl,node,flag

  integer ne_loc,base,top,loc_index

  if(ne_spbc.le.ncpus) then
    if(myid+1.le.ne_spbc) then
      ne_loc=1
    else
      ne_loc=0
    end if
  else
     base=floor(real(ne_spbc)/real(ncpus))
     top=ne_spbc-base*ncpus
     if(myid+1.le.top) then
        ne_loc=base+1
     else
        ne_loc=base
     end if
  end if


  I_fluid_temp(:)=I_fluid_center(:)

!  do loc_index=1,ne_loc
!     ie=myid+1+(loc_index-1)*ncpus
  do ie=1,ne_spbc
     do icount=1,nn_contact_remove
        if(spbcele(ie)==contact_remove(icount)) then
	   do inl=1,8
	      node=ne-ne_spbc+(ie-1)*8+inl
	      if(I_fluid_center(node).lt.0.4) I_fluid_temp(node)=1.0
	      if(I_fluid_center(node).gt.0.6) I_fluid_temp(node)=0.0
	   end do
	end if
     end do
  end do

  do icount=1,ne_den_domain
     ie=den_domain(icount)
     I_fluid_temp(ie)=I_fluid_center(ie)
  end do
  I_fluid_center(:)=I_fluid_temp(:)

end subroutine set_center_after_contact
