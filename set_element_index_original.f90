!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!differetiate inner outer elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine set_element_index(I_var,ne_inter,ne_inner,ne_outer,inter_ele,inner_ele,outer_ele,ien,I_fluid_center)
  use fluid_variables,only:nen,nn,ne

  integer ie,je,icount,inl,node
  integer ne_inter,ne_inner,ne_outer
  integer inter_ele(ne),inner_ele(ne),outer_ele(ne)
  integer ien(nen,ne)
  real(8) I_var(nn)
  real(8) flag(nen)
  real(8) norm1
  real(8) eps
  real(8) I_fluid_center(ne)
  integer temp

  I_fluid_center(:)=0.5
  eps = 1.0e-3
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
write(*,*)'ne_inter_ini=',ne_inter
!write(*,*)inter_ele(1:ne_inter)
temp=ne_inter
  do ie=1,ne_outer
     do inl=1,nen
	node=ien(inl,outer_ele(ie))
	I_var(node)=0.0
     end do
  end do

  do ie=1,ne_inner
     do inl=1,nen
	node=ien(inl,inner_ele(ie))
	if(abs(I_var(node)).le.1.0e-6) then
	   ne_inter=ne_inter+1
	   inter_ele(ne_inter)=inner_ele(ie)
	end if
	I_var(node)=1.0
     end do
     I_fluid_center(inner_ele(ie))=1.0
  end do ! set the indicator of the nodes of the inner elements to be 1

  do ie=1,ne_outer
     do inl=1,nen
	node=ien(inl,outer_ele(ie))
	if (abs(I_var(node)-1.0).le.1.0e-6) then
	   ne_inter=ne_inter+1
	   inter_ele(ne_inter)=outer_ele(ie)
	end if
	I_var(node)=0.0
     end do
     I_fluid_center(outer_ele(ie))=0.0
  end do !set the indicator of the nodes of the outer elements to be 0
write(*,*)'ne_inter_after=',ne_inter
!write(*,*)inter_ele(1:ne_inter)
!ne_inter=temp
end subroutine set_element_index












