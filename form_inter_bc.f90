!!!!!!!!!!!

!!!!!!!!!!set the Indicator of the nodes of the boundary to be 0

subroutine form_inter_bc(I_var,bcnode_den,flag,id_den)

!  use fluid_variables, only:nn,nen,ne,neface
!  use centermesh_variables
  use denmesh_variables, only:nbc_den,nn_den
  real(8) I_var(nn_den)
  integer i,inl,node
  real(8) flag
  integer bcnode_den(nbc_den)
  integer id_den(nn_den)

  do i=1,nbc_den
     I_var(bcnode_den(i))=flag
     id_den(bcnode_den(i))=1
  end do
return
end subroutine
