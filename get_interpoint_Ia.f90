!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Get initial indicator for interface points!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_interpoint_Ia(xp,xg,Ip,Ig,hg,infdomain_sub,I_fluid_center,x_center)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd

  real(8) xp(nsd,maxmatrix)          !coordinates of interface points
  real(8) xg(nsd,nn)                 !coordinates of fluid points
  real(8) Ip(maxmatrix)              !indicator of interface points
  real(8) Ig(nn)                     !indicator of fluid points
  real(8) hg(ne)                     !spacing
  integer infdomain_sub(maxmatrix)  !the elements that each inter points belongs to
  real(8) I_fluid_center(ne)
  real(8) x_center(nsd,ne)

  integer i,j
  real(8) dx(nsd)
  real(8) hs
  real(8) Sp
  real(8) length
!  I_fluid_center(:)=1.0
 
  Ip(1:nn_inter) = 0.0
  do i=1,nn_inter
     do j=1,ne
        hs=hg(j)
!	if (abs(Ig(j)).gt.1.0e-6) then
	   dx(1:nsd)=abs(xp(:,i)-x_center(:,j))
	   call B_Spline(dx,hs,nsd,Sp)
!	   call B_Spline_inter(xp(:,i),x_center(:,j),nsd,hs,'00',Sp)
!	   length=sqrt((xp(1,i)-x_center(1,j))**2.0+(xp(2,i)-x_center(2,j))**2.0)
!	   call B_Spline_0order(length,hs,Sp)
	   Ip(i)=Ip(i)+I_fluid_center(j)*Sp
!	end if
     end do
  end do
  Ig(1:nn)=0.0
  do i=1,nn
     do j=1,ne
	hs=hg(j)
	dx(1:nsd)=abs(xg(:,i)-x_center(:,j))
	call B_Spline(dx,hs,nsd,Sp)
!	call B_Spline_inter(xg(:,i),x_center(:,j),nsd,hs,'00',Sp)
!           length=sqrt((xg(1,i)-x_center(1,j))**2.0+(xg(2,i)-x_center(2,j))**2.0)
!           call B_Spline_0order(length,hs,Sp)
	Ig(i)=Ig(i)+I_fluid_center(j)*Sp
     end do
  end do
!write(*,*)'Ip(1:nn_inter)=',Ip(1:nn_inter)
!write(*,*)'Ig(1:nn)=',Ig(1:nn)
!stop
end subroutine get_interpoint_Ia









