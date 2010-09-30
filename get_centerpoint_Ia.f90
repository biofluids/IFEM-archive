

subroutine get_centerpoint_Ia(x,I_fluid,x_center,I_fluid_center,corr_Ip,x_inter,hg,infdomain_inter,ne_inter_temp,inter_ele)

  use interface_variables
  use fluid_variables,only:nsd,nn,ne

  real(8) x_center(nsd,ne)
  real(8) I_fluid_center(ne)
  real(8) corr_Ip(maxmatrix)
  real(8) x_inter(nsd,maxmatrix)
  real(8) hg(ne)
  integer infdomain_inter(maxmatrix)
  real(8) I_fluid(nn)
  real(8) x(nsd,nn)
  integer ne_inter_temp
  integer inter_ele(ne)

  integer i,j,icount,node,inl,flag
  real(8) dx(nsd)
  real(8) hs,Sp,temp(ne_inter_temp)
  real(8) residual_center

  flag=1  
if(flag==1) then
  do i=1,ne_inter_temp
     temp(i)=0.0
     do j=1,nn_inter
	hs=hg(infdomain_inter(j))
	dx(:)=abs(x_center(:,inter_ele(i))-x_inter(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	temp(i)=temp(i)+corr_Ip(j)*Sp
     end do
    
     I_fluid_center(inter_ele(i))=I_fluid_center(inter_ele(i))+temp(i)
!     write(*,*)'I_fluid_center',inter_ele(i),'=',I_fluid_center(inter_ele(i))
  end do
!  write(*,*)'I_fluid_center=',I_fluid_center(:)
!  call getnorm(temp,temp,ne_inter_temp,residual_center)
!  residual_center=sqrt(residual_center)
!write(*,*)'residual_center=',residual_center
else if(flag==2) then
  do i=1,ne
     temp(i)=0.0
     hs=hg(infdomain_inter(i))
     do j=1,nn
	dx(:)=abs(x_center(:,i)-x(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	if(I_fluid(j).gt.1.0)then
	  I_fluid(j)=1.0
	else if(I_fluid(j).lt.0.0)then
	  I_fluid(j)=0.0
	end if
	temp(i)=temp(i)+Sp*I_fluid(j)
     end do
  end do
  call getnorm(temp(:)-I_fluid_center(:),temp(:)-I_fluid_center(:),ne,residual_center)
  I_fluid_center(:)=temp(:)
    residual_center=sqrt(residual_center)
write(*,*)'residual_center=',residual_center

end if
end subroutine get_centerpoint_Ia








  
