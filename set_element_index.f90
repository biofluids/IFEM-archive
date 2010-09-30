!=====================================================
!differentiate inner and outer elements
!====================================================

subroutine set_element_index(hg,x_center,Ic_inter,ne_inter,ne_inner,ne_outer,inter_ele,inner_ele,outer_ele,I_fluid_center,ien,ne_inter_temp)

  use fluid_variables, only:nen,ne,nn,nsd

  integer ne_inter,ne_inner,ne_outer
  integer inter_ele(ne),inner_ele(ne),outer_ele(ne)
  integer ne_inner_tmp,ne_outer_tmp
  integer inner_ele_tmp(ne),outer_ele_tmp(ne)
  integer ien(nen,nn)
  real(8) I_fluid_center(ne)
  real(8) eps
  integer I_var(nn)
  integer ne_inter_temp
  real(8) Ic_inter
  real(8) hs,sp,dx(nsd)
  real(8) x_center(nsd,ne),hg(ne)

  integer ie,je,icount,inl,node,flag
!write(*,*)'ne_inter=',ne_inter
!write(*,*)'I_fluid_center=',I_fluid_center(1:ne)
!stop
  eps=1.0e-2
  ne_inner=0
  ne_outer=0
  I_var(:)=0
  do ie=1,ne
     do je=1,ne_inter
	if(ie==inter_ele(je))then
	  I_fluid_center(ie)=0.0
	  do inl=1,nen
	     node=ien(inl,ie)
	     I_var(node)=1
	  end do
	  goto 100
	end if
     end do

     if(abs(I_fluid_center(ie)-1.0).le.eps)then
	ne_inner=ne_inner+1
	inner_ele(ne_inner)=ie
	I_fluid_center(ie)=1.0
     else 
	ne_outer=ne_outer+1
	outer_ele(ne_outer)=ie
	I_fluid_center(ie)=-1.0
     end if

100 continue
  end do
  do ie=1,ne_inter
     hs=hg(inter_ele(ie))
     do je=1,ne
	if(je.ne.inter_ele(ie)) then
	  dx(1:nsd)=abs(x_center(1:nsd,inter_ele(ie))-x_center(1:nsd,je))
	  call B_Spline(dx,hs,nsd,Sp)
	  I_fluid_center(inter_ele(ie))=I_fluid_center(inter_ele(ie))+I_fluid_center(je)*Sp
	end if
     end do
!     write(*,*)'I_fluid_center_interele=',inter_ele(ie),I_fluid_center(inter_ele(ie))
  end do
!  stop
do ie=1,ne_inter
  I_fluid_center(inter_ele(ie))=0.0
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
!	ne_inter=ne_inter+1
!	inter_ele(ne_inter)=inner_ele(ie)
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
!	ne_inter=ne_inter+1
!	inter_ele(ne_inter)=outer_ele(ie)
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








