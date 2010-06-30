!!!!!!!!!!!

!!!!!!!!!!set the Indicator of the nodes of the boundary to be 0

subroutine form_inter_bc(I_var_c,rng,ien,flag)

  use fluid_variables, only:nn,nen,ne,neface

  real(8) I_var_c(nn)
  integer rng(neface,ne)
  integer ien(nen,ne)
  integer i,inl,node
  real(8) flag

  do i=1,ne
     do inl=1,nen
	if (rng(inl,i) .gt. 0) then
	   node=ien(inl,i)
	   I_var_c(node) = flag
	end if
     end do
  end do

return
end subroutine
