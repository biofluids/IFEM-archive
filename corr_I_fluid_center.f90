!=====

!correct I_fluid center for interfacial elements

!====

subroutine corr_I_fluid_center(I_fluid_center,ne_inter,x_inter,x_center,inter_ele,hg,corr_Ip)

  use interface_variables
  use fluid_variables,only:nsd,ne,nn

  real(8) I_fluid_center(ne)
  integer ne_inter
  real(8) x_inter(nsd,maxmatrix),x_center(nsd,ne)
  integer inter_ele(ne)
  real(8) hg(ne)
  real(8) corr_Ip(maxmatrix)

  integer i,j,ie,inl,node,icount,je
  real(8) I_fluid_center_corr(ne_inter)
  real(8) dx(nsd),Sp,hs

  I_fluid_center_corr(:)= 0.0

  do ie=1,ne_inter
     hs=hg(inter_ele(ie))
     do je=1,ne
	dx(1:nsd)=abs(x_center(1:nsd,inter_ele(ie))-x_center(1:nsd,je))
	call B_Spline(dx,hs,nsd,Sp)
	I_fluid_center_corr(ie)=I_fluid_center_corr(ie)+I_fluid_center(je)*Sp
     end do
  end do

  do ie=1,ne_inter
     hs=hg(inter_ele(ie))
     do icount=1,nn_inter
	dx(1:nsd)=abs(x_center(1:nsd,inter_ele(ie))-x_inter(1:nsd,icount))
	call B_Spline(dx,hs,nsd,Sp)
	I_fluid_center_corr(ie)=I_fluid_center_corr(ie)+corr_Ip(icount)*Sp
     end do
  end do

  do ie=1,ne_inter
     if(I_fluid_center_corr(ie).gt.0.0) then
!        I_fluid_center(inter_ele(ie))=1.0
     else
!	I_fluid_center(inter_ele(ie))=-1.0
     end if
     I_fluid_center(inter_ele(ie))=I_fluid_center_corr(ie)
  end do
  
end subroutine corr_I_fluid_center

