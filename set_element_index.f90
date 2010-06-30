!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!differetiate inner outer elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine set_element_index(I_var,ne_inter,ne_inner,ne_outer,inter_ele,inner_ele,outer_ele,ien)
  use fluid_variables,only:nen,nn,ne

  integer ie,je,icount,inl,node
  integer ne_inter,ne_inner,ne_outer
  integer inter_ele(nn),inner_ele(nn),outer_ele(nn)
  integer ien(nen,ne)
  real(8) I_var(nn)
  real(8) flag(nen)
  real(8) norm1
  real(8) eps

  eps = 1.0e-4
  ne_inner = 0
  ne_outer = 0
  do ie=1,ne
     do je=1,ne_inter
	if (ie==inter_ele(je))then
	   go to 100
	end if
     end do
     
     do inl=1,nen
	node=ien(inl,ie)
	flag(inl)=I_var(node)-1.0
     end do

     call getnorm(flag,flag,nen,norm1)

     norm1=sqrt(norm1)

     if (norm1 .le. eps)then
	ne_inner = ne_inner+1
	inner_ele(ne_inner)=ie
     

     else if (norm1 .gt. eps)then
	ne_outer = ne_outer+1
	outer_ele(ne_outer)=ie
     end if

100 continue
  end do!!!!!!!!!!differentiate inner and outer elements

  do ie=1,ne_outer
     do inl=1,nen
	node=ien(inl,outer_ele(ie))
	I_var(node)=0.0
     end do
  end do !set the indicator of the nodes of the outer elements to be 0

end subroutine set_element_index












