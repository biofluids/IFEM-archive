!!!!!!!!!!!!!!!Get interface elements!!!!!!!!!!!!!!


subroutine get_inter_ele(infdomain_sub,inter_ele,ne_inter)

  use interface_variables,only:nn_inter,maxmatrix
  use fluid_variables,only:nn

  integer infdomain_sub(maxmatrix)
  integer inter_ele(nn)
  integer i,j,icount,flag
  integer ne_inter

  ncount = 0
  inter_ele(:) = 0

  do i=1,nn_inter
     do j=1,ncount
	if (infdomain_sub(i) == inter_ele(j)) then
	   go to 100
	end if
     end do

	ncount = ncount + 1
	inter_ele(ncount) = infdomain_sub(i)

100 continue
  end do
  ne_inter = ncount

return
end subroutine get_inter_ele









