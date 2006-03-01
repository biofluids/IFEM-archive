	sh(0,1) = sq(0,1,iq) 
	sh(0,2) = sq(0,2,iq) 
	sh(0,3) = sq(0,3,iq) 
	sh(0,4) = sq(0,4,iq) 

	xr(1,1) = sq(1,1,iq) * x(1,1) + sq(1,2,iq) * x(1,2) &
		    + sq(1,3,iq) * x(1,3) + sq(1,4,iq) * x(1,4)
	xr(1,2) = sq(1,1,iq) * x(2,1) + sq(1,2,iq) * x(2,2) &
		    + sq(1,3,iq) * x(2,3) + sq(1,4,iq) * x(2,4)
	xr(2,1) = sq(2,1,iq) * x(1,1) + sq(2,2,iq) * x(1,2) &
		    + sq(2,3,iq) * x(1,3) + sq(2,4,iq) * x(1,4)
	xr(2,2) = sq(2,1,iq) * x(2,1) + sq(2,2,iq) * x(2,2) &
		    + sq(2,3,iq) * x(2,3) + sq(2,4,iq) * x(2,4)

!  jacobian
    cf(1,1) = + (xr(1,1)*xr(2,2) - xr(2,1)*xr(1,2))
      
    det = cf(1,1)

	sx(1,1) = xr(2,2)/det
	sx(1,2) = -xr(1,2)/det
	sx(2,1) = -xr(2,1)/det
	sx(2,2) = xr(1,1)/det

!  global first derivatives

    sh(1,1)= sq(1,1,iq) * sx(1,1) + sq(2,1,iq) * sx(1,2) 
    sh(1,2)= sq(1,2,iq) * sx(1,1) + sq(2,2,iq) * sx(1,2) 
    sh(1,3)= sq(1,3,iq) * sx(1,1) + sq(2,3,iq) * sx(1,2) 
	sh(1,4)= sq(1,4,iq) * sx(1,1) + sq(2,4,iq) * sx(1,2) 
    sh(2,1)= sq(1,1,iq) * sx(2,1) + sq(2,1,iq) * sx(2,2) 
    sh(2,2)= sq(1,2,iq) * sx(2,1) + sq(2,2,iq) * sx(2,2) 
    sh(2,3)= sq(1,3,iq) * sx(2,1) + sq(2,3,iq) * sx(2,2) 
	sh(2,4)= sq(1,4,iq) * sx(2,1) + sq(2,4,iq) * sx(2,2) 

 !write(*,*) 'shape=',sh(0,1),sh(0,2),sh(0,3),sh(0,4)
 !write(*,*) 'shape derivatives=', sh(1,1),sh(1,2),sh(1,3), sh(1,4)
 ! write(*,*) '2nd shape derivatives=', sh(2,1),sh(2,2),sh(2,3), sh(2,4)
 !stop