!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Get initial indicator for interface points!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_interpoint_Ia(xp,xg,Ip,Ig,hg,infdomain_sub)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd

  real(8) xp(nsd,maxmatrix)          !coordinates of interface points
  real(8) xg(nsd,nn)                 !coordinates of fluid points
  real(8) Ip(maxmatrix)              !indicator of interface points
  real(8) Ig(nn)                     !indicator of fluid points
  real(8) hg(ne)                     !spacing
  integer infdomain_sub(maxmatrix)  !the elements that each inter points belongs to

  integer i,j
  real(8) dx(nsd)
  real(8) hs
  real(8) Sp

  Ip(1:nn_inter) = 0.0
  do i=1,nn_inter
     hs=hg(infdomain_sub(i))
     do j=1,nn
!	if (Ig(j).gt.1.0e-6) then
	   dx(1:nsd)=abs(xp(:,i)-xg(:,j))
	   call B_Spline(dx,hs,nsd,Sp)
	   Ip(i)=Ip(i)+Ig(j)*Sp
!	end if
     end do
  end do
end subroutine get_interpoint_Ia









