!!!!!!!!!!!!!!!!

!!!!!!!!!set the Indicator of the nodes of interfacial elements to be 1

subroutine form_inter_ele(inter_ele,ne_inter,I_var,ien,flag)

  use interface_variables,only: nn_inter,maxmatrix
  use fluid_variables, only:nn,nen,ne

  integer ne_inter
  integer inter_ele(nn)
  integer ien(nen,ne)
  real(8) I_var(nn)

  integer i,j
  integer inl,node
  real(8) flag

  do i=1,ne_inter
     do inl=1,nen
	node=ien(inl,inter_ele(i))
	I_var(node) = flag
     enddo
  end do

end subroutine form_inter_ele



