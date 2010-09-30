!============================================

!get the velocity for interface points using fem interpolation

!===========================================

subroutine get_intervel_fem(x,x_inter_regen,nn_inter_regen,infdomain_inter,vel_fluid,vel_inter,ien)

  use fluid_variables,only:nsd,nn,ne,nen
  use interface_variables

  real(8) x(nsd,nn)                    !coordinates of fluid nodes
  real(8) x_inter_regen(nsd,maxmatrix)       !coordinates of interfacial points
  integer infdomain_inter(maxmatrix)   !the element that each point belongs to
  real(8) vel_fluid(nsd,nn)            !velocity of fluid nodes
  real(8) vel_inter(nsd,maxmatrix)     !velocity of interfacial points
  integer ien(nen,ne)
  integer nn_inter_regen

  integer i,j,node,isd,inl,ie
  real(8) x_loc(nsd,nen)               !coordinates of local element
  real(8) x_inter_loc(nsd)             !coordinate of an interfacial points
  real(8) sh(nen)                      !shape function

  vel_inter(:,:) = 0.0
  do i=1,nn_inter_regen
     x_inter_loc(:)=x_inter_regen(:,i)
     do inl=1,nen
	node=ien(inl,infdomain_inter(i))
	x_loc(:,inl)=x(:,node)
     end do

     call sh_exchange(x_inter_loc,x_loc,nsd,nen,sh)

     do inl=1,nen
	node=ien(inl,infdomain_inter(i))
	vel_inter(:,i)=vel_inter(:,i)+vel_fluid(:,node)*sh(inl)
     end do
  end do
!write(*,*)'vel_inter=',sqrt(vel_inter(1,1:nn_inter)**2+vel_inter(2,1:nn_inter)**2)
end subroutine get_intervel_fem

