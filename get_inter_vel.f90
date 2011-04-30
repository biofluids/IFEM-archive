!===============================
!get interface velocity
!=================================

subroutine get_inter_vel(x,x_inter,infdomain,vel_fluid,vel_inter,hg)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix)
  integer infdomain(maxmatrix)
  real(8) vel_fluid(nsd,nn),vel_inter(nsd,maxmatrix)
  real(8) hg(ne)

  integer i,j
  real(8) dx(nsd),Sp,hs
  
  vel_inter(:,:)=0.0

  do i=1,nn_inter
     hs=hg(infdomain(i))
     do j=1,nn
	dx(:)=abs(x_inter(:,i)-x(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	vel_inter(:,i)=vel_inter(:,i)+vel_fluid(:,j)*Sp
     end do
  end do

end subroutine get_inter_vel
