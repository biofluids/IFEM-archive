	sh(0,1) = sq2d(0,1,iq) 
	sh(0,2) = sq2d(0,2,iq) 
	sh(0,3) = sq2d(0,3,iq) 
	sh(0,4) = sq2d(0,4,iq) 

	xr(1,1) = sq2d(1,1,iq) * x(1,1) + sq2d(1,2,iq) * x(1,2) &
		    + sq2d(1,3,iq) * x(1,3) + sq2d(1,4,iq) * x(1,4)
	xr(1,2) = sq2d(1,1,iq) * x(2,1) + sq2d(1,2,iq) * x(2,2) &
		    + sq2d(1,3,iq) * x(2,3) + sq2d(1,4,iq) * x(2,4)
	xr(2,1) = sq2d(2,1,iq) * x(1,1) + sq2d(2,2,iq) * x(1,2) &
		    + sq2d(2,3,iq) * x(1,3) + sq2d(2,4,iq) * x(1,4)
	xr(2,2) = sq2d(2,1,iq) * x(2,1) + sq2d(2,2,iq) * x(2,2) &
		    + sq2d(2,3,iq) * x(2,3) + sq2d(2,4,iq) * x(2,4)

c  jacobian
    cf(1,1) = + (xr(1,1)*xr(2,2) - xr(2,1)*xr(1,2))
      
    det = cf(1,1)

	sx(1,1) = xr(2,2)/det
	sx(1,2) = -xr(1,2)/det
	sx(2,1) = -xr(2,1)/det
	sx(2,2) = xr(1,1)/det

c  global first derivatives

    sh(1,1)= sq2d(1,1,iq) * sx(1,1) + sq2d(2,1,iq) * sx(1,2) 
    sh(1,2)= sq2d(1,2,iq) * sx(1,1) + sq2d(2,2,iq) * sx(1,2) 
    sh(1,3)= sq2d(1,3,iq) * sx(1,1) + sq2d(2,3,iq) * sx(1,2) 
	sh(1,4)= sq2d(1,4,iq) * sx(1,1) + sq2d(2,4,iq) * sx(1,2) 
    sh(2,1)= sq2d(1,1,iq) * sx(2,1) + sq2d(2,1,iq) * sx(2,2) 
    sh(2,2)= sq2d(1,2,iq) * sx(2,1) + sq2d(2,2,iq) * sx(2,2) 
    sh(2,3)= sq2d(1,3,iq) * sx(2,1) + sq2d(2,3,iq) * sx(2,2) 
	sh(2,4)= sq2d(1,4,iq) * sx(2,1) + sq2d(2,4,iq) * sx(2,2) 

