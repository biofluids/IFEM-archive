!!!!!!!!!!!!!!!Get interface elements!!!!!!!!!!!!!!


subroutine get_inter_ele(infdomain_sub,inter_ele,ne_inter,inter_ele_nn)

  use interface_variables,only:nn_inter,maxmatrix
  use fluid_variables,only:nn,ne

  integer infdomain_sub(maxmatrix)
  integer inter_ele(ne)
  integer i,j,icount,flag
  integer ne_inter
  integer inter_ele_nn(ne)

  ncount = 0
  inter_ele(:) = 0
  inter_ele_nn(:) = 1
  do i=1,nn_inter
     do j=1,ncount
	if (infdomain_sub(i) == inter_ele(j)) then
	   inter_ele_nn(j)=inter_ele_nn(j)+1
	   go to 100
	end if
     end do

	ncount = ncount + 1
	inter_ele(ncount) = infdomain_sub(i)

100 continue
  end do
  ne_inter = ncount
!  do i=1,ne_inter
!     do j=1,nn_inter
!	if(infdomain_sub(j)==inter_ele(i)) then
!	   inter_ele_nn(i)=inter_ele_nn(i)+1
!	end if
!     end do
!  end do
!write(*,*)'inter_ele_nn=',inter_ele_nn(1:ne_inter)

!write(*,*)'inter_ele=',inter_ele(1:ne_inter)
return
end subroutine get_inter_ele









