!================================================
!used to calculate the derivative of the indicator for a point
!=================================================

subroutine get_indicator_derivative(x,xp,x_center,infdomain,hg,I_fluid_center,corr_Ip, &
				  sum_0d,sum_1d,sum_2d, &
				  II,dI,ddI)

  use interface_variables
  use fluid_variables,only:nsd,ne,nn
!  use mpi_variables

  real(8) x(nsd),xp(nsd,maxmatrix),x_center(nsd,ne)
  integer infdomain(maxmatrix)
  real(8) hg(ne)
  real(8) I_fluid_center(ne),corr_Ip(maxmatrix)
  real(8) sum_0d(2),sum_1d(2,nsd),sum_2d(2,3*(nsd-1))
  real(8) II,dI(nsd),ddI(3*(nsd-1))

  integer i,j,isd,jsd
  real(8) hs,S(nsd),Sp(nsd),Spp(nsd),temp,dx
  real(8) W0,W1(nsd),W2(3*(nsd-1))
  integer flag, jcount

!used to calculate the 1st&2nd derivative of indicator
  II=0.0
  dI(:)=0.0
  ddI(:)=0.0
  do flag=1,2

     if(flag==1) then
	jcount=ne
     else if(flag==2) then
	jcount=nn_inter
     end if

  do j=1,jcount
	if(flag==1) then
	   hs=hg(j)
	else if(flag==2) then
	   hs=hg(infdomain(j))
	end if
     do isd=1,nsd
	do jsd=1,nsd  !loop over jsd,if jsd=isd, do the derivative

	  if(flag==1) then
	   dx=x(jsd)-x_center(jsd,j)
	  else if(flag==2) then
	   dx=x(jsd)-xp(jsd,j)
	  end if

	   if(jsd==isd) then
		call B_Spline_1order(dx,hs,Sp(jsd))
		call B_Spline_2order(dx,hs,Spp(jsd))
	   else
		call B_Spline_0order(dx,hs,Sp(jsd))
		Spp(jsd)=Sp(jsd)
	   end if
	end do

	W1(isd)=1.0
	W2(isd)=1.0
	do jsd=1,nsd
	   W1(isd)=W1(isd)*Sp(jsd)
	   W2(isd)=W2(isd)*Spp(jsd)
	end do   !calculate S1x, S1xx
     end do

     if(nsd==2) then
	do isd=1,nsd
	   if(flag==1) then
		dx=x(isd)-x_center(isd,j)
	   else if(flag==2) then
		dx=x(isd)-xp(isd,j)
	   end if
	   call B_Spline_1order(dx,hs,Sp(isd))
	   call B_Spline_0order(dx,hs,Spp(isd))
	end do

	W0=Spp(1)*Spp(2)
	W2(3)=Sp(1)*Sp(2)
     end if

     if(flag==1) then  !for center points
	II=II+W0/sum_0d(flag)*I_fluid_center(j)

	dI(1:nsd)=dI(1:nsd)+(W1(1:nsd)/sum_0d(flag)-W0*sum_1d(flag,1:nsd)/sum_0d(flag)**2)*I_fluid_center(j)

	ddI(1:nsd)=ddI(1:nsd)+((W2(1:nsd)*sum_0d(flag)-W1(1:nsd)*sum_1d(flag,1:nsd))/sum_0d(flag)**2- &
		       (W1(1:nsd)*sum_1d(flag,1:nsd)+W0*sum_2d(flag,1:nsd))/sum_0d(flag)**2+ &
		       2.0*W0*sum_1d(flag,1:nsd)**2/sum_0d(flag)**3)*I_fluid_center(j)

	ddI(3)=ddI(3)+( (W2(3)*sum_0d(flag)-W1(1)*sum_1d(flag,2))/sum_0d(flag)**2- &
		        (W1(2)*sum_1d(flag,1)+W0*sum_2d(flag,3))/sum_0d(flag)**2+ &
		        2.0*W0*sum_1d(flag,1)*sum_1d(flag,2)/sum_0d(flag)**3)*I_fluid_center(j)

!	II=II+W0*I_fluid_center(j)
!        dI(1:nsd)=dI(1:nsd)+W1(1:nsd)*I_fluid_center(j)
!        ddI(1:nsd)=ddI(1:nsd)+W2(1:nsd)*I_fluid_center(j)
!        ddI(3)=ddI(3)+W2(3)*I_fluid_center(j)


     else if(flag==2) then ! for inter points

!        II=II+W0/sum_0d(flag)*corr_Ip(j)
!
!        dI(1:nsd)=dI(1:nsd)+(W1(1:nsd)/sum_0d(flag)-W0*sum_1d(flag,1:nsd)/sum_0d(flag)**2)*corr_Ip(j)
!
!        ddI(1:nsd)=ddI(1:nsd)+((W2(1:nsd)*sum_0d(flag)-W1(1:nsd)*sum_1d(flag,1:nsd))/sum_0d(flag)**2- &
!                       (W1(1:nsd)*sum_1d(flag,1:nsd)+W0*sum_2d(flag,1:nsd))/sum_0d(flag)**2+ &
!                       2.0*W0*sum_1d(flag,1:nsd)**2/sum_0d(flag)**3)*corr_Ip(j)
!
!        ddI(3)=ddI(3)+( (W2(3)*sum_0d(flag)-W1(1)*sum_1d(flag,2))/sum_0d(flag)**2- &
!                        (W1(2)*sum_1d(flag,1)+W0*sum_2d(flag,3))/sum_0d(flag)**2+ &
!                        2.0*W0*sum_1d(flag,1)*sum_1d(flag,2)/sum_0d(flag)**3)*corr_Ip(j)



	II=II+W0*corr_Ip(j)
	dI(1:nsd)=dI(1:nsd)+W1(1:nsd)*corr_Ip(j)
	ddI(1:nsd)=ddI(1:nsd)+W2(1:nsd)*corr_Ip(j)
	ddI(3)=ddI(3)+W2(3)*corr_Ip(j)

     end if
   end do! end of j loop
  end do ! end of flag loop

end subroutine get_indicator_derivative
























