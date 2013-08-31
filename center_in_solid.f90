
  subroutine center_in_solid(I_solid,center_mapping,ien,its)

  use interface_variables, only:nn_center,ele_refine,I_solid_inter
  use fluid_variables,only:nsd,nen,ne,nn
  use allocate_variables,only:flag_center_in_solid
  use mpi_variables
  include 'mpif.h'

  real(8) I_solid(nn)
  integer center_mapping(ne,ele_refine**nsd+1)
  integer ien(nen,ne),its
  integer i,j,k,node,ie,flag,inl

  if(allocated(flag_center_in_solid)) deallocate(flag_center_in_solid)
  allocate(flag_center_in_solid(nn_center))

  flag_center_in_solid(:)=0

  do ie=1,ne
     flag=0
     do inl=1,nen
        node=ien(inl,ie)
        if(I_solid(node).gt.I_solid_inter) flag=flag+1
     end do
     if(flag==nen) then
       j=center_mapping(ie,1)
       do i=2,j+1
          node=center_mapping(ie,i)
          flag_center_in_solid(node)=1
       end do
     end if
  end do
  end subroutine center_in_solid
