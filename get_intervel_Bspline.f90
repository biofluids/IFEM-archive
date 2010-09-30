!=============================================================

!get the velocity for interfacial points using Bspline function

!=============================================================

subroutine get_intervel_Bspline(x,x_inter,infdomain_inter,vel_fluid,vel_inter,hg)

  use fluid_variables,only:nsd,nn,ne,nen
  use interface_variables

  real(8) x(nsd,nn)                   !coordinates of fluid nodes
  real(8) x_inter(nsd,maxmatrix)      !coordinates of interfacial points
  integer infdomain_inter(maxmatrix)
  real(8) vel_fluid(nsd,nn)           !velocity of fluid nodes
  real(8) vel_inter(nsd,maxmatrix)    !velocity of interfacial points
  real(8) hg(ne)

  integer i,j,node,isd,inl
  real(8) dx(nsd),Sp,hs
  
  vel_inter(:,:)=0.0

  do i=1,nn_inter
     hs=hg(infdomain_inter(i))
     do j=1,nn
	dx(:)=abs(x_inter(:,i)-x(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	vel_inter(:,i)=vel_inter(:,i)+vel_fluid(:,j)*Sp
     end do
  end do
     
end subroutine get_intervel_Bspline
