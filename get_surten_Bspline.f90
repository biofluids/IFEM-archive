!=============================================

!get surface tension for fluid nodes

!=============================================

subroutine get_surten_Bspline(x,x_inter,norm_inter,curv_inter,arc_inter,sur_fluid,hg,infdomain_inter,I_fluid,Ic_inter)

  use fluid_variables, only:nsd,nn,ne,nen,den_liq
  use interface_variables

  real(8) x(nsd,nn)                 !coordinates for fluid nodes
  real(8) x_inter(nsd,maxmatrix)    !corrdinates for interfacial points
  real(8) norm_inter(nsd,maxmatrix) !normal for interfacial points
  real(8) curv_inter(maxmatrix)     !curvature for interfacial points
  real(8) arc_inter(maxmatrix)      !arclength for interfacial points
  real(8) sur_fluid(nsd,nn)             !surface tension for fluid nodes
  real(8) I_fluid(nn)
  real(8) Ic_inter
  real(8) hg(ne)
  integer infdomain_inter(maxmatrix)

  integer i,j,isd
  real(8) hs,Sp,dx(nsd)
  real(8) delta
  real(8) den_p,den_f
  real(8) temp

  sur_fluid(:,:) = 0.0
  den_p=den_liq+(den_inter-den_liq)*Ic_inter
  do i=1,nn

     do j=1,nn_inter
	hs=hg(infdomain_inter(j))
	dx(:)=abs(x(:,i)-x_inter(:,j))
!	call getnorm(dx,dx,nsd,temp)
!	temp=sqrt(temp)
!	call B_Spline_0order(temp,hs,Sp)
!	delta=Sp/hs

	delta=1.0
	do isd=1,nsd
	   temp=dx(isd)
	   call B_Spline_0order(temp,hs,Sp)
	   delta=delta*Sp
	end do     !get delta function
!	delta=delta/(hs**nsd)
!	delta=delta/hs**2
	sur_fluid(:,i)=sur_fluid(:,i)+arc_inter(j)*sur_tension*curv_inter(j)*norm_inter(:,j)*delta
     end do
     den_f=den_liq+(den_inter-den_liq)*I_fluid(i)
     sur_fluid(:,i)=den_f/den_p*sur_fluid(:,i)
  end do
!     sur_fluid(:,:)=-sur_fluid(:,:)
!write(*,*)'sur_fluid=',sur_fluid(:,:)
!stop
end subroutine get_surten_Bspline




















     
