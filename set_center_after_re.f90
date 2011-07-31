

subroutine set_center_after_re(I_fluid_den,I_fluid_center,I_fluid,ien)

  use mpi_variables
  use fluid_variables, only:ne,nen,nn
  use denmesh_variables,only:ne_den,nen_den,nn_den
  use allocate_variables, only:den_domain, ne_den_domain

  real(8) I_fluid_center(ne)
  real(8) I_fluid(nn)
  real(8) I_fluid_den(nn_den)
  integer ien(nen,ne)

  integer ie,je,icount,inl,node,flag

!  do ie=1,ne
!     if(I_fluid_center(ie).ge.0.5) then
!	I_fluid_center(ie)=1.0
!	I_fluid_den(ien(1:nen,ie))=1.0
!     else
!	I_fluid_center(ie)=0.0
!	I_fluid_den(ien(1:nen,ie))=0.0
!     end if
!  end do

  do ie=1,ne_den
     I_fluid_center(ie)=0.0
     do inl=1,nen_den
	I_fluid_center(ie)=I_fluid_center(ie)+0.25*I_fluid_den(ien(inl,ie))
     end do
     if(I_fluid_center(ie).gt.0.5) then
	I_fluid_center(ie)=1.0
     else
	I_fluid_center(ie)=0.0
     end if
  end do
     
	
  do icount=1,ne_den_domain
     ie=den_domain(icount)
     I_fluid_center(ie)=0.0
     do inl=1,nen
	I_fluid_center(ie)=I_fluid_center(ie)+0.25*I_fluid(ien(inl,ie))
     end do
  end do

end subroutine set_center_after_re
