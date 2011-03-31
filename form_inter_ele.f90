!!!!!!!!!!!!!!!!

!!!!!!!!!set the Indicator of the nodes of interfacial elements to be 1

subroutine form_inter_ele(inter_ele_den,ne_inter_den,I_var,ien_den,flag,id_den)

  use interface_variables,only: nn_inter,maxmatrix
!  use fluid_variables, only:nn,nen,ne
!  use centermesh_variables
  use denmesh_variables, only:nn_den,nen_den,ne_den

  integer ne_inter_den
  integer inter_ele_den(ne_den)
  integer ien_den(nen_den,ne_den)
  real(8) I_var(nn_den)
  integer id_den(nn_den)

  integer i,j
  integer inl,node
  real(8) flag

  do i=1,ne_inter_den
     do inl=1,nen_den
	node=ien_den(inl,inter_ele_den(i))
	I_var(node) = flag
	id_den(node)=1
     enddo
!     I_var(inter_ele(i))=flag
  end do
end subroutine form_inter_ele



