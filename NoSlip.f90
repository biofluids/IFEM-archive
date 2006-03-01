subroutine Noslip(nn_solid,x_solids,x_fluid,d_fluid,solid_force)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Mickael 03/02/06
  ! No-slip boundary conditions
  ! ATTENTION UNIFORM MESH
  
  use fluid_variables
  use run_variables, only: dt,tt
  use delta_nonuniform
  implicit none

  !real(8),allocatable :: shrknode(:,:)      !...shape function for each node, contains the weights

  integer :: nn_solid !,ndelta

  real(8) :: x_solids(nsd,nn_solid)
  real(8) :: x_fluid(nsd,nn),x(nsd)
  real(8) :: d_fluid(ndf,nn)
  real(8),dimension(1:nsd,1:nn_solid) :: solid_force   !...fluid structure interaction force


 !...local variables
  integer :: i,j,k,inf(nn),ninf,infJ(nn),ninfJ
  real(8) :: deltaX,Radius,normal(nsd,nn_solid),distX(nn_solid),distY(nn_solid)
  real(8) :: Uk(ndf,nn_solid),Press(4,nn_solid),Vf(2,nn_solid),Vfluid(2,4,nn_solid)

  real(8) :: Sol_Vel(ndf,nn_solid),Fluid_Vel(ndf,nn)

  real(8) :: xk(5,nn_solid),yk(5,nn_solid)
  real(8) :: dVxdx(nn_solid), dVxdy(nn_solid), dVydx(nn_solid),dVydy(nn_solid)
  real(8) :: dPdx(nn_solid),dPdy(nn_solid)
  real(8) :: d2Vxdx2(nn_solid),d2Vxdy2(nn_solid),d2Vydx2(nn_solid),d2Vydy2(nn_solid)

  real(8) :: Tempx1,Tempx2,Tempxk,Tempy3,Tempy4,Tempyk
  real(8) :: Temp2x1,Temp2x2,Temp2xk,Temp2y3,Temp2y4,Temp2yk

  real(8) :: Accel_Force(nsd,nn_solid), Press_Force(nsd,nn_solid)
  real(8) :: Iner_Force(nsd,nn_solid), Visc_Force(nsd,nn_solid)
  real(8) :: shp,b(4),bd(nsd,4)
  integer :: n,nnum,isd
  real(8) :: adist(nsd,nn)
  real(8) :: dwjp(nn),y(nsd),a(nsd),temp(nn)
  real(8) :: Magn_Accel(nn_solid),Magn_Iner(nn_solid),Magn_Press(nn_solid),Magn_Visc(nn_solid),Magn_Total(nn_solid)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Output force terms for point A (node96) and B (node71)
  open(unit = 500,file = "AcceForceA.plt",status='unknown')
  open(unit = 501,file = "AcceForceB.plt",status='unknown')
  open(unit = 502,file = "InerForceA.plt",status='unknown')
  open(unit = 503,file = "InerForceB.plt",status='unknown')
  open(unit = 504,file = "PressForceA.plt",status='unknown')
  open(unit = 505,file = "PressForceB.plt",status='unknown')
  open(unit = 506,file = "ViscForceA.plt",status='unknown')
  open(unit = 507,file = "ViscForceB.plt",status='unknown')
  open(unit = 508,file = "TotalForceA.plt",status='unknown')
  open(unit = 509,file = "TotalForceB.plt",status='unknown')
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  deltaX= abs(x_fluid(1,1)-x_fluid(1,2))

  !
  ! nodes along interface for 2D disk only:
  ! 2-5, 27-56, 65-82, 92-103

 
  Radius=0.25 ! Radius Disk

  Uk(3,1:nn_solid)=0.0
  Vf(:,:)=0.0
  dwjp(1:nn)=5.0d-3  ! Not sure about this value. TO BE CHECKED

  do k=1,nn_solid  
  	! Boundary nodes
	inf(1:nn_solid)=0

	if ((k.ge.2 .and. k.le.5) .or. (k.ge.27 .and. k.le.56) .or. &
				(k.ge.65 .and. k.le.82) .or. (k.ge.92 .and. k.le.103) ) then

		solid_force(1:nsd,k)=0.0
		
		! Velcoity at node k (no-slip)
		Uk(1,k)=0.0
		Uk(2,k)=0.0


		! Position boundary nodes and their 4 influence nodes
		xk(5,k)=x_solids(1,k)
		yk(5,k)=x_solids(2,k)

		if (yk(5,k) .lt. 1.0) then
											
			if (xk(5,k) .lt. 1.0) then
														
				xk(1,k)=x_solids(1,k)-deltaX
				yk(1,k)=x_solids(2,k)
				xk(2,k)=x_solids(1,k)-2.0*deltaX
				yk(2,k)=x_solids(2,k)
				xk(3,k)=x_solids(1,k)
				yk(3,k)=x_solids(2,k)-deltaX
				xk(4,k)=x_solids(1,k)
				yk(4,k)=x_solids(2,k)-2.0*deltaX
			
			else
				xk(1,k)=x_solids(1,k)+deltaX
				yk(1,k)=x_solids(2,k)
				xk(2,k)=x_solids(1,k)+2.0*deltaX
				yk(2,k)=x_solids(2,k)
				xk(3,k)=x_solids(1,k)
				yk(3,k)=x_solids(2,k)-deltaX
				xk(4,k)=x_solids(1,k)
				yk(4,k)=x_solids(2,k)-2.0*deltaX
			endif

		else
														
			if (xk(5,k) .lt. 1.0) then
				xk(1,k)=x_solids(1,k)-deltaX
				yk(1,k)=x_solids(2,k)
				xk(2,k)=x_solids(1,k)-2.0*deltaX
				yk(2,k)=x_solids(2,k)
				xk(3,k)=x_solids(1,k)
				yk(3,k)=x_solids(2,k)+deltaX
				xk(4,k)=x_solids(1,k)
				yk(4,k)=x_solids(2,k)+2.0*deltaX				
	
			else
				xk(1,k)=x_solids(1,k)+deltaX
				yk(1,k)=x_solids(2,k)
				xk(2,k)=x_solids(1,k)+2.0*deltaX
				yk(2,k)=x_solids(2,k)
				xk(3,k)=x_solids(1,k)
				yk(3,k)=x_solids(2,k)+deltaX
				xk(4,k)=x_solids(1,k)
				yk(4,k)=x_solids(2,k)+2.0*deltaX
								
			endif

		endif


		Sol_Vel(1:ndf,1:4)=0.0
		
		do j=1,4	! 4 nodes around each boundary nodes
			ninfJ=0
			temp(j)=0.0
			do i=1,nn

				! Fluid nodes inside influence domain
				if ((x_fluid(1,i) .le. ((xk(j,k)+2.0*deltaX))) .and. &
						(x_fluid(1,i) .ge. ((xk(j,k)-2.0*deltaX))))  then

					if ((x_fluid(2,i) .le. ((yk(j,k)+2.0*deltaX))) .and. &
							(x_fluid(2,i) .ge. ((yk(j,k)-2.0*deltaX)))) then
			
						! Remove nodes inside solid 
						if (((x_fluid(1,i)-1)**2+(x_fluid(2,i)-1)**2) .gt. (Radius**2)) then
			  				! get a list of influence nodes from the fluids grid
							ninfJ=ninfJ+1	!number of influence nodes
							infJ(ninfJ)=i
							
							! Store Velocity and pressure of influence nodes of nodes 1, 2, 3 and 4
							Fluid_Vel(1,ninfJ)=d_fluid(1,i)
							Fluid_Vel(2,ninfJ)=d_fluid(2,i)
							Fluid_Vel(3,ninfJ)=d_fluid(3,i)
					
						endif
					endif
				endif

			enddo


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Interpolate velocity and pressure for points j=1,4
			! Get solids point J coordinates 
			x(1)=xk(j,k)
			x(2)=yk(j,k)
		
			adist(1:nsd,1:nn)=deltaX ! TO BE CHECKED

			! calculate the correction function
			call correct2dl(b,bd,x,x_fluid,adist,dwjp,nn,infJ,ninfJ,maxconn)

			do n = 1, ninfJ
				nnum = infJ(n)

				do isd = 1,nsd
					y(isd) = x_fluid(isd,nnum)
					a(isd) = adist(isd,nnum)
				enddo

				call RKPMshape2d(shp,b,bd,x,y,a,dwjp(nnum))

				shrknode(n,j)=shp
	
				Sol_Vel(1:ndf,j) = Sol_Vel(1:ndf,j) + Fluid_Vel(1:ndf,n) * shrknode(n,j)

			enddo
	

		enddo



		!===========================
		! Identify nodes for pressue at point k				  
		! normal for a sphere
		if (yk(5,k) .lt. 1.0) then
														
			if (xk(5,k) .lt. 1.0) then
										
				normal(1,k)=-(-(xk(5,k)-1))/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
				normal(2,k)=(yk(5,k)-1)/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
							
			else
									
				normal(1,k)=-(-(xk(5,k)-1))/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
				normal(2,k)=(yk(5,k)-1)/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
									
			endif

		else
														
			if (xk(5,k) .lt. 1.0) then
				normal(1,k)=-(-(xk(5,k)-1))/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
				normal(2,k)=(yk(5,k)-1)/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
									
			else
								
				normal(1,k)=-(-(xk(5,k)-1))/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
				normal(2,k)=(yk(5,k)-1)/sqrt((-(xk(5,k)-1))**2+(yk(5,k)-1)**2)
								
			endif

		endif

		distX(k)=deltaX*normal(1,k)
		distY(k)=deltaX*normal(2,k)

		ninf=0					
		do i=1,nn
			! Influence nodes for pressure at boundary (node k)
			if (abs(x_fluid(1,i)-(xk(5,k)+distX(k))) .lt. deltaX) then
				if (abs(x_fluid(2,i)-(yk(5,k)+distY(k))) .lt. deltaX) then
					! get a list of influence nodes from the fluids grid
					ninf=ninf+1	!number of influence nodes
					inf(ninf)=i
					Press(ninf,k)=d_fluid(3,i)
					Vfluid(1,ninf,k)=d_fluid(1,i)
					Vfluid(2,ninf,k)=d_fluid(2,i)
				
				endif
			endif
		enddo
			
	
	

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Interpolate pressure at solid node k
		! Get solids point coordinate
		x(1)=xk(5,k)+distX(k)
		x(2)=yk(5,k)+distY(k)
		
		adist(1:nsd,1:nn)=deltaX ! TO BE CHECKED

		! calculate the correction function
		call correct2dl(b,bd,x,x_fluid,adist,dwjp,nn,inf,ninf,maxconn)
		  
		temp(k)=0.0d0

		do n = 1, ninf
			nnum = inf(n)
			do isd = 1,nsd
				y(isd) = x_fluid(isd,nnum)
				a(isd) = adist(isd,nnum)
			 enddo
		
			 
			call RKPMshape2d(shp,b,bd,x,y,a,dwjp(nnum))
	

			shrknode(n,k)=shp

			! Interpolate pressure at node k		
			Uk(3,k) = Uk(3,k) + Press(n,k)*shrknode(n,k)
			Vf(1,k) = Vf(1,k) + Vfluid(1,n,k)*shrknode(n,k)
			Vf(2,k) = Vf(2,k) + Vfluid(2,n,k)*shrknode(n,k)

			temp(k)=temp(k)+shp

		enddo


		!=================================================================================
		! Compute gradient and divergence at node k (2 velocity components and pressure)

		! Gradient

		Tempx1=(xk(5,k)-xk(2,k))/((xk(1,k)-xk(2,k))*(xk(1,k)-xk(5,k)))
		Tempx2=(xk(5,k)-xk(1,k))/((xk(2,k)-xk(1,k))*(xk(2,k)-xk(5,k)))
		Tempxk=((xk(5,k)-xk(1,k))+(xk(5,k)-xk(2,k)))/((xk(5,k)-xk(1,k))*(xk(5,k)-xk(2,k)))

		Tempy3=(yk(5,k)-yk(4,k))/((yk(3,k)-yk(4,k))*(yk(3,k)-yk(5,k)))
		Tempy4=(yk(5,k)-yk(3,k))/((yk(4,k)-yk(3,k))*(yk(4,k)-yk(5,k)))
		Tempyk=((yk(5,k)-yk(3,k))+(yk(5,k)-yk(4,k)))/((yk(5,k)-yk(3,k))*(yk(5,k)-yk(4,k)))

		! dVx/dx
		dVxdx(k)=Tempx1*Sol_Vel(1,1)+Tempx2*Sol_Vel(1,2)+Tempxk*Uk(1,k)
		! dVx/dy
		dVxdy(k)=Tempy3*Sol_Vel(1,3)+Tempy4*Sol_Vel(1,4)+Tempyk*Uk(1,k)

		! dVy/dx
		dVydx(k)=Tempx1*Sol_Vel(2,1)+Tempx2*Sol_Vel(2,2)+Tempxk*Uk(2,k)
		! dVy/dy
		dVydy(k)=Tempy3*Sol_Vel(2,3)+Tempy4*Sol_Vel(2,4)+Tempyk*Uk(2,k)

		! dP/dx
		dPdx(k)=Tempx1*Sol_Vel(3,1) +Tempx2*Sol_Vel(3,2) +Tempxk*Uk(3,k)
		! dP/dy
		dPdy(k)=Tempy3*Sol_Vel(3,3)+Tempy4*Sol_Vel(3,4) +Tempyk*Uk(3,k)

		! Divergence

		Temp2x1=2/((xk(1,k)-xk(2,k))*(xk(1,k)-xk(5,k)))
		Temp2x2=2/((xk(2,k)-xk(1,k))*(xk(2,k)-xk(5,k)))
		Temp2xk=2/((xk(5,k)-xk(1,k))*(xk(5,k)-xk(2,k)))

		Temp2y3=2/((yk(3,k)-yk(4,k))*(yk(3,k)-yk(5,k)))
		Temp2y4=2/((yk(4,k)-yk(3,k))*(yk(4,k)-yk(5,k)))
		Temp2yk=2/((yk(5,k)-yk(3,k))*(yk(5,k)-yk(4,k)))
	
		! d2Vx/dx2
		d2Vxdx2(k)=Temp2x1*Sol_Vel(1,1)+Temp2x2*Sol_Vel(1,2)+Temp2xk*Uk(1,k)
		! d2Vx/dy2
		d2Vxdy2(k)=Temp2y3*Sol_Vel(1,3)+Temp2y4*Sol_Vel(1,4)+Temp2yk*Uk(1,k)

		! d2Vy/dx2
		d2Vydx2(k)=Temp2x1*Sol_Vel(2,1)+Temp2x2*Sol_Vel(2,2)+Temp2xk*Uk(2,k)
		! d2Vy/dy2
		d2Vydy2(k)=Temp2y3*Sol_Vel(2,3)+Temp2y4*Sol_Vel(2,4)+Temp2yk*Uk(2,k)

		! Acceleration force components and magnitude	
		Accel_Force(1,k)=-den_liq*Vf(1,k)/dt
		Accel_Force(2,k)=-den_liq*Vf(2,k)/dt
		Magn_Accel(k)=sqrt(Accel_Force(1,k)**2+Accel_Force(2,k)**2)
		
		! Inertia force components and magnitude
		Iner_Force(1,k)=den_liq*(dVxdx(k)+dVydy(k))*Uk(1,k)
		Iner_Force(2,k)=den_liq*(dVxdx(k)+dVydy(k))*Uk(2,k)
		Magn_Iner(k)=sqrt(Iner_Force(1,k)**2+Iner_Force(2,k)**2)

		! Viscous force components and magnitude
		Visc_Force(1,k)=-vis_liq*(d2Vxdx2(k)+d2Vxdy2(k))
		Visc_Force(2,k)=-vis_liq*(d2Vydx2(k)+d2Vydy2(k))
		Magn_Press(k)=sqrt(Visc_Force(1,k)**2+Visc_Force(2,k)**2)

		! Pressure force components and magnitude
		Press_Force(1,k)=dPdx(k)
		Press_Force(2,k)=dPdy(k)
		Magn_Visc(k)=sqrt(Press_Force(1,k)**2+Press_Force(2,k)**2)



		solid_force(1:nsd,k)= solid_force(1:nsd,k)+ &
					Accel_Force(1:nsd,k) + &
					Iner_Force(1:nsd,k) + &
					Press_Force(1:nsd,k)  + &
					Visc_Force(1:nsd,k)
		
		if (k==71) then
			write(501, 3000) tt, 	Accel_Force(1,k), 	Accel_Force(2,k), 	Magn_Accel(k)
			write(503, 3000) tt, 	Iner_Force(1,k), 	Iner_Force(2,k), 	Magn_Iner(k)	
			write(505, 3000) tt, 	Press_Force(1,k), 	Press_Force(2,k), 	Magn_Press(k)
			write(507, 3000) tt, 	Visc_Force(1,k), 	Visc_Force(2,k),	Magn_Visc(k)
			write(509, 3000) tt,    solid_Force(1,k),        solid_Force(2,k),      Magn_Total(k)
		elseif (k==96) then
			write(500, 3000) tt, 	Accel_Force(1,k), 	Accel_Force(2,k), 	Magn_Accel(k)
                        write(502, 3000) tt, 	Iner_Force(1,k), 	Iner_Force(2,k), 	Magn_Iner(k)
                        write(504, 3000) tt, 	Press_Force(1,k), 	Press_Force(2,k), 	Magn_Press(k)
                        write(506, 3000) tt, 	Visc_Force(1,k), 	Visc_Force(2,k),	Magn_Visc(k)
			write(508, 3000) tt,    solid_Force(1,k),        solid_Force(2,k),      Magn_Total(k)
		endif	
		3000 format(6e12.4,  2X,  6e12.4,  2X,  6e12.4,  2X,  6e12.4)


	endif

  enddo

 return
end subroutine NoSlip

