!=====================================================
!differentiate inner and outer elements
!====================================================

subroutine set_element_index(infdomain_den,ne_inter,ne_inner,ne_outer,inter_ele,inner_ele,outer_ele,I_fluid_center,I_fluid_den,ien_den,ien,inter_ele_den,ne_inter_den)

  use fluid_variables, only:nen,ne,nn,nsd
  use denmesh_variables, only:nn_den,ne_den,nen_den
  integer infdomain_den(ne)
  integer ne_inter,ne_inner,ne_outer
  integer inter_ele(ne),inner_ele(ne),outer_ele(ne)
  integer ne_inner_tmp,ne_outer_tmp
  integer inner_ele_tmp(ne),outer_ele_tmp(ne)
  integer ien(nen,nn)
  real(8) I_fluid_center(ne), I_fluid_den(nn_den)
  integer ien_den(nen_den,nn_den)
  integer ne_inter_den,inter_ele_den(ne_den)
  real(8) eps
  integer I_var(nn)

  integer ie,je,icount,inl,node,flag
  integer flag_den(ne_den)
!write(*,*)'ne_inter=',ne_inter
!write(*,*)'I_fluid_center=',I_fluid_center(1:ne)
!stop
  eps=1.0e-4
  ne_inner=0
  ne_outer=0
  I_var(:)=0
  do ie=1,ne_den
     do je=1,ne_inter_den
	if(ie==inter_ele_den(je))then
	  flag_den(ie)=0
	  goto 100
	end if
     end do
     flag=0
     do inl=1,nen_den
        node=ien_den(inl,ie)
	if(abs(I_fluid_den(node)-1.0).gt.eps) then
	  flag=1
	end if
     end do

     if(flag==0) then
	flag_den(ie)=1
     else
	flag_den(ie)=-1
     end if
100 continue
  end do

  do ie=1,ne
     if(flag_den(infdomain_den(ie))==1) then
	I_fluid_center(ie)=1.0
!	ne_inner=ne_inner+1
!	inner_ele(ne_inner)=i
     else if(flag_den(infdomain_den(ie))==0) then
	I_fluid_center(ie)=0.5
!	do inl=1,nen
!	   node=ien(inl,i)
!	   I_var(node)=1
!	end do
     else
	I_fluid_center(ie)=0.0
!	ne_outer=ne_outer+1
!	outer_ele(ne_outer)=i
     end if
     do je=1,ne_inter
	if(ie==inter_ele(je)) then
	  do inl=1,nen
	     node=ien(inl,ie)
	     I_var(node)=1
	  end do
	  goto 200
	end if
    end do
    if(abs(I_fluid_center(ie)-1.0).lt.eps) then
	ne_inner=ne_inner+1
	inner_ele(ne_inner)=ie
    end if
    if(abs(I_fluid_center(ie)).lt.eps) then
	ne_outer=ne_outer+1
	outer_ele(ne_outer)=ie
    end if

200 continue

  end do

write(*,*)'ne_inter=',ne_inter
  ne_inter_temp=ne_inter
  ne_inner_tmp=0
  do ie=1,ne_inner
     flag=0
     do inl=1,nen
	node=ien(inl,inner_ele(ie))
	if(I_var(node)==1) then
	  flag=1
	end if
     end do
     if(flag==1) then
	ne_inter=ne_inter+1
	inter_ele(ne_inter)=inner_ele(ie)
     else if(flag==0) then
	ne_inner_tmp=ne_inner_tmp+1
	inner_ele_tmp(ne_inner_tmp)=inner_ele(ie)
     end if
  end do
  ne_outer_tmp=0
  do ie=1,ne_outer
     flag=0
     do inl=1,nen
	node=ien(inl,outer_ele(ie))
	if(I_var(node)==1)then
	  flag=1
	end if
     end do
     if(flag==1) then
	ne_inter=ne_inter+1
	inter_ele(ne_inter)=outer_ele(ie)
     else if(flag==0) then
	ne_outer_tmp=ne_outer_tmp+1
	outer_ele_tmp(ne_outer_tmp)=outer_ele(ie)
     end if
  end do
write(*,*)'ne_inter_after=',ne_inter
!  ne_inner=ne_inner_tmp
!  ne_outer=ne_outer_tmp
!  inner_ele(:)=inner_ele_tmp(:)
!  outer_ele(:)=outer_ele_tmp(:)
!do ie=1,ne_inter
!   I_fluid_center(inter_ele(ie))=0.0
!end do
end subroutine set_element_index








