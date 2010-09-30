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
  integer signf,ncount
  integer nn_flag(nn),corr_ele(ne),ne_corr,flag_inner,flag_outer

  I_fluid_center(:)=0.5
  eps = 1.0e-6
  ne_inner = 0
  ne_outer = 0
  nn_flag(1:nn)=0
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
  I_var(1:nn)=0.0
  do ie=1,ne_inter
     do inl=1,nen
        node=ien(inl,inter_ele(ie))
	I_var(node)=1.0
     end do 
  end do   ! set the nodes of interfacial elements to be 1.0

  ne_corr=0
  do ie=1,ne_inner
     signf=0
     ncount=0
     do inl=1,nen
        node=ien(inl,inner_ele(ie))
        if(abs(I_var(node)-1.0).lt.1.0e-6)then
	  ncount=ncount+1
	  signf=1
	end if
     end do
     if(signf==1) then
       ne_inter=ne_inter+1
       inter_ele(ne_inter)=inner_ele(ie)
     end if
     if(ncount==nen) then
!	I_fluid_center(inner_ele(ie))=0.5
	write(*,*)'correct I_fluid_center for inner element'
!when all the nodes of an inner element are actually the nodes of interfacial elements,set the center indicator to be 0.5
       ne_corr=ne_corr+1
       corr_ele(ne_corr)=inner_ele(ie)
     else
        I_fluid_center(inner_ele(ie))=1.0
        do inl=1,nen
	   node=ien(inl,inner_ele(ie))
	   nn_flag(node)=2
	end do
     end if
  end do     

  do ie=1,ne_outer
     signf=0
     do inl=1,nen
	node=ien(inl,outer_ele(ie))
	nn_flag(node)=1
	if(abs(I_var(node)-1.0).lt.1.0e-6)then
	  signf=1
	end if
     end do
     if (signf==1)then
	ne_inter=ne_inter+1
	inter_ele(ne_inter)=outer_ele(ie) 
     end if
     I_fluid_center(outer_ele(ie))=0.0
  end do
!===set the elements next to the interfacial elements to be the interfacial elements, and set the element center indicator of inner and outer element to be 1 and 0.

  do ie=1,ne_corr
     flag_inner=0
     flag_outer=0
     do inl=1,nen
	node=ien(inl,corr_ele(ie))
	if(nn_flag(node)==1)then
	  flag_outer=1
	end if
	if(nn_flag(node)==2)then
	  flag_inner=1
	end if
     end do
     if((flag_inner==1).and.(flag_outer==1))then
	I_fluid_center(corr_ele(ie))=0.5
	write(*,*)'correct I_fluid_center to 0.5'
     else if((flag_inner==1).and.(flag_outer==0))then
	I_fluid_center(corr_ele(ie))=1.0
	write(*,*)'correct I_fluid_center to 1.0'
     else if((flag_inner==0).and.(flag_outer==1))then
	I_fluid_center(corr_ele(ie))=0.0
	write(*,*)'correct I_fluid_Center to 0.0'
     else
	I_fluid_center(corr_ele(ie))=0.0
	write(*,*)'mess up element'
     end if
  end do
write(*,*)'ne_inter_after=',ne_inter
!write(*,*)inter_ele(1:ne_inter)
!ne_inter=temp
!stop
end subroutine set_element_index












