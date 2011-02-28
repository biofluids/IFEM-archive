!========================================
!get normal considering weight summation
!========================================

subroutine get_weight_derivative(x,xp,x_center,infdomain,hg,&
		     sum_0d,sum_1d,sum_2d)

  use interface_variables
  use fluid_variables, only:nsd,ne,nn
!  use mpi_variables

  real(8) x(nsd),xp(nsd,maxmatrix), x_center(nsd,ne)
  integer infdomain(maxmatrix)
  real(8) hg(ne)
  real(8) sum_0d(2),sum_1d(2,nsd),sum_2d(2,3*(nsd-1))
!Sp_sum_0d=S1+S2+...
!Sp_sum_1d=S1x+S2x+...,S1y+S2y+...,S1z+S2z+...
!Sp_sum_2d=S1xx+S2xx+...,S1yy+S2yy+...,(S1zz+S2zz+...);
!          S1xy+S2xy+...,(S1xz+S2xz+...,S1yz+S2yz+...)

  integer i,j,isd,jsd
  real(8) hs,Sp(nsd),temp,dx,Spp(nsd),temp_0

!  do i=1,nn_inter
     sum_0d(:)=0.0
     sum_1d(:,:)=0.0
     sum_2d(:,:)=0.0
     do isd=1,nsd  ! loop over isd
!=========================================================================
!used for calculate the summation of the 1st order derivative of the weight
	do j=1,ne  ! loop over center nodes
	   hs=hg(j)
	   do jsd=1,nsd  ! loop over jsd,if jsd=isd, do the derivative
	      dx=x(jsd)-x_center(jsd,j)
	      if(jsd==isd) then
		call B_Spline_1order(dx,hs,Sp(jsd))
		call B_Spline_2order(dx,hs,Spp(jsd))
	      else
		call B_Spline_0order(dx,hs,Sp(jsd))
		Spp(jsd)=Sp(jsd)
	      end if
	   end do ! get 1st&2nd order derivative for isd
	   temp=1.0
	   do jsd=1,nsd
	      temp=temp*Sp(jsd)
	   end do
	   sum_1d(1,isd)=sum_1d(1,isd)+temp  !S1x+S2x+...
	   temp=1.0
	   do jsd=1,nsd
	      temp=temp*Spp(jsd)
	   end do
	   sum_2d(1,isd)=sum_2d(1,isd)+temp  !S1xx+S2xx+...    isd==x
	end do
	do j=1,nn_inter ! loop over inter points
	   hs=hg(infdomain(j))
	   do jsd=1,nsd
	      dx=x(jsd)-xp(jsd,j)
	      if(jsd==isd) then
		call B_Spline_1order(dx,hs,Sp(jsd))
		call B_Spline_2order(dx,hs,Spp(jsd))
	      else
		call B_Spline_0order(dx,hs,Sp(jsd))
		Spp(jsd)=Sp(jsd)
	      end if
	   end do
	   temp=1.0
	   do jsd=1,nsd
	      temp=temp*Sp(jsd)
	   end do
!	   sum_1d(2,isd)=sum_1d(2,isd)+temp
	   temp=1.0
	   do jsd=1,nsd
	      temp=temp*Spp(jsd)
	   end do
!	   sum_2d(2,isd)=sum_2d(2,isd)+temp
	end do
! calculate the summation of the 1st&2nd order derivative of the weight for each dof
! for each inter point     
    end do ! end of loop over isd
!=====================================================
  if(nsd==2) then
	do j=1,ne
	   hs=hg(j)
	   do isd=1,nsd
		dx=x(isd)-x_center(isd,j)
		call B_Spline_1order(dx,hs,Sp(isd))
		call B_Spline_0order(dx,hs,Spp(isd))
	   end do
	   sum_2d(1,3)=sum_2d(1,3)+Sp(1)*Sp(2)
	   sum_0d(1)=sum_0d(1)+Spp(1)*Spp(2)
	end do
	do j=1,nn_inter
	   hs=hg(infdomain(j))
	   do isd=1,nsd
		dx=x(isd)-xp(isd,j)
		call B_Spline_1order(dx,hs,Sp(isd))
		call B_Spline_0order(dx,hs,Spp(isd))
	   end do
!	   sum_2d(2,3)=sum_2d(2,3)+Sp(1)*Sp(2)
!	   sum_0d(2)=sum_0d(2)+Spp(1)*Spp(2)
	end do
   else
!	if(myid==0) then
	  write(*,*)'no 3d for weight derivative'
	  stop
!	end if
  end if
!used for calcualte S1xy+S2xy+...
!===================================================

end subroutine get_weight_derivative







