!=============================================
!regen points in problematic elements
!=============================================

subroutine regen_in_ele(x,x_inter,x_center,Ic_inter,I_fluid_center,ien, &
			infdomain_inter,hg,norm_inter,curv_inter, &
			regen_point_flag,regen_ele,nn_regen_ele)

  use fluid_variables, only:nsd,ne,nn,nen
  use interface_variables

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,ne)
  integer infdomain_inter(maxmatrix)
  real(8) hg(ne),ien(nen,ne)
  real(8) norm_inter(nsd,maxmatrix)
  real(8) curv_inter(maxmatrix)
  real(8) Ic_inter, I_fluid_center(ne)
  integer regen_point_flag(maxmatrix)
  integer regen_ele(maxmatrix),nn_regen_ele

  real(8) corr_Ip(maxmatrix)
  integer i,j,icount,jcount,node,ie,inl,isd
  integer flag
  integer nn_inter_ele_regen
  real(8) x_inter_ele_regen(nsd,maxmatrix)
  
  integer nn_inter_new
  real(8) x_inter_new(nsd,maxmatrix)
  integer regen_point_flag_new(maxmatrix)


  call get_indicator(x_inter,x,hg,infdomain_inter,I_fluid_center,x_center,corr_Ip,Ic_inter)
  call get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,infdomain_inter, &
				norm_inter,curv_inter,hg)
!if(nn_regen_ele==0) then
!  goto 100
!end if

write(*,*)'begin element_wise regeneration'
write(*,*)'=============================='

do isd=1,5
  do i=1,nn_inter
     if(abs(curv_inter(i)).gt.curv_bound) then
	flag=0
	do j=1,nn_regen_ele
	   if(regen_point_flag(i)==regen_ele(j)) then
	     flag=1
	   end if
	end do
	if(flag==0) then
	   nn_regen_ele=nn_regen_ele+1
	   regen_ele(nn_regen_ele)=regen_point_flag(i)
	end if
     end if
  end do   ! if curv>max, find the element index and push it the regen_ele
write(*,*)'nn-regen_ele=',nn_regen_ele
if(nn_regen_ele==0) then
  goto 100
end if
  nn_inter_new=0
  do i=1,nn_regen_ele
     call points_regen_in_ele(x,x_inter,x_center,x_inter_ele_regen, &
			Ic_inter,nn_inter_ele_regen, &
			I_fluid_center,corr_Ip,hg,infdomain_inter,ien, &
			regen_ele(i))
	
	nn_inter_new=nn_inter_new+nn_inter_ele_regen
	x_inter_new(1:nsd,(nn_inter_new-nn_inter_ele_regen+1):nn_inter_new)= &
			x_inter_ele_regen(1:nsd,1:nn_inter_ele_regen)
!write(*,*)'x_inter_ele_regen',x_inter_ele_regen(1:nsd,1:nn_inter_ele_regen)
	regen_point_flag_new((nn_inter_new-nn_inter_ele_regen+1):nn_inter_new)= &
			regen_ele(i)
  end do
!write(*,*)'x_inter',x_inter(1:nsd,1:nn_inter)  
  
  do i=1,nn_inter
	flag=0
	do j=1,nn_regen_ele
	   if(regen_point_flag(i)==regen_ele(j)) then
		flag=1
	   end if
	end do
	if(flag==0) then
	   nn_inter_new=nn_inter_new+1
	   x_inter_new(1:nsd,nn_inter_new)=x_inter(1:nsd,i)
	   regen_point_flag_new(nn_inter_new)=regen_point_flag(i)
	end if
  end do

  nn_inter=nn_inter_new
  x_inter(1:nsd,1:nn_inter)=x_inter_new(1:nsd,1:nn_inter_new)
write(*,*)'nn_inter=',nn_inter
!write(*,*)'x_inter',x_inter(1:nsd,1:nn_inter)

  regen_point_flag(1:nn_inter)=regen_point_flag_new(1:nn_inter_new)


    call search_inf_inter(x_inter,x,nn,nn_inter,nsd,ne,nen,ien,infdomain_inter)

  call get_indicator(x_inter,x,hg,infdomain_inter,I_fluid_center,x_center,corr_Ip,Ic_inter)
  call get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,infdomain_inter, &
                                norm_inter,curv_inter,hg)


  nn_regen_ele=0
  regen_ele(:)=0

end do
100 continue
write(*,*)'============================================================'
end subroutine


     

  


  
