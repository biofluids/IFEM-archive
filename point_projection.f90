!=====================================================================
!use newton iteration to project candidate points to the interface
!=====================================================================

subroutine point_projection(Ic_inter,xloc,x_center,I_fluid_center,x_inter,corr_Ip,hs,I_can,err_p,hg,infdomain_inter)

  use fluid_variables,only:nsd,nn,ne,nen
  use interface_variables

  real(8) Ic_inter              !constant interfaceial indicator
  real(8) xloc(nsd)             !coordinates for candidate point
  real(8) x_center(nsd,ne)          !coordinates for element center
  real(8) I_fluid_center(ne)    !initial indicator for element center
  real(8) x_inter(nsd,maxmatrix)!coordinates for interfacial points
  real(8) corr_Ip(maxmatrix)    !correction term
  real(8) hs                    !spacing
  real(8) I_can                 !indicator for candidate point
  real(8) err_p,err_p1(2)
  real(8) hg(ne)
  integer infdomain_inter(maxmatrix)

  integer i,j,isd,jsd,iit,nit
  real(8) dI(nsd)               !dI/dx,dI/dy
  real(8) dx,Sp(nsd),temp,ddx(nsd)      
  real(8) delta(nsd)            !delta x and delta y

  nit = 1
  err_p=1.0
!  do iit=1,nit      !begin newton interation
  do while((nit.le.5).and.(err_p.gt.1.0e-8))
!=========get the derivative of the indicator for the candidate point========
     dI(:)=0.0
     do i=1,ne      !loop over fluid centers
        hs=hg(i)
	do isd=1,nsd    ! loop over isd
	   do jsd=1,nsd ! loop over jsd,if jsd=isd, do the derivative
              dx=xloc(jsd)-x_center(jsd,i)
	      if(jsd==isd) then
		 call B_Spline_1order(dx,hs,Sp(jsd))
	      else
		 call B_Spline_0order(dx,hs,Sp(jsd))
	      end if
	    end do  ! get 1st derivative of Bspline function for isd
	    temp=1.0
	    do jsd=1,nsd
	       temp=temp*Sp(jsd)
	    end do
	    dI(isd)=dI(isd)+I_fluid_center(i)*temp
	end do  ! end of loop over isd
      end do    ! end of loop over fluid centers

     do i=1,nn_inter  ! loop over interfacial points
        hs=hg(infdomain_inter(i))/1.0
	do isd=1,nsd
	   do jsd=1,nsd
	      dx=xloc(jsd)-x_inter(jsd,i)
	      if(jsd==isd) then
		 call B_Spline_1order(dx,hs,Sp(jsd))
	      else
		 call B_Spline_0order(dx,hs,Sp(jsd))
	      end if
	   end do
	   temp=1.0
	   do jsd=1,nsd
	      temp=temp*Sp(jsd)
	   end do
	   dI(isd)=dI(isd)+corr_Ip(i)*temp
	end do
     end do   ! end of loop over interfacial points
!=============================================================================
!=====get delta x and delta y=================================================

     if (nsd==2) then
	temp=(dI(1)**2+dI(2)**2)/dI(1)
	delta(1)=(Ic_inter-I_can)/temp
	delta(2)=delta(1)*dI(2)/dI(1)
     else if(nsd==3) then
	write(*,*)'no 3D right now'
	stop
     end if
!===============================================================================
     xloc(1:nsd)=xloc(1:nsd)+delta(1:nsd)
!==============================================================================
! write(*,*)'I_can-Ic_inter=',I_can-Ic_inter
!     err_p=abs(I_can-Ic_inter)
     I_can=0.0
     do i=1,ne
        hs=hg(i)
	ddx(:)=abs(xloc(:)-x_center(:,i))
	call B_Spline(ddx,hs,nsd,temp)
	I_can=I_can+I_fluid_center(i)*temp
     end do
     do i=1,nn_inter
	hs=hg(infdomain_inter(i))/1.0
	ddx(:)=abs(xloc(:)-x_inter(:,i))
	call B_Spline(ddx,hs,nsd,temp)
	I_can=I_can+corr_Ip(i)*temp
     end do
 err_p1(1)=abs(I_can-Ic_inter)
 err_p1(2)=sqrt(delta(1)**2+delta(2)**2)
 err_p=maxval(err_p1)
!write(*,*)'delta=',delta
nit=nit+1
  end do ! end of loop for newton iteration
!write(*,*)'I_can=',I_can
end subroutine point_projection


