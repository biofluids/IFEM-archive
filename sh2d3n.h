      sh(0,1) = sq(0,1,iq)
      sh(0,2) = sq(0,2,iq)
      sh(0,3) = sq(0,3,iq)

      xr(1,1) = x(1, 1) - x(1, 3) 
      xr(1,2) = x(2, 1) - x(2, 3) 
      xr(2,1) = x(1, 2) - x(1, 3) 
      xr(2,2) = x(2, 2) - x(2, 3) 
    
!  jacobian
      cf(1,1) = + (xr(1,1)*xr(2,2) - xr(2,1)*xr(1,2))
      
      det = cf(1,1)

!  global first derivatives

      sh(1,1)= xr(2,2)/det
      sh(1,2)= -xr(1,2)/det
      sh(1,3)= -xr(2,2)/det + xr(1,2)/det
      sh(2,1)= -xr(2,1)/det
      sh(2,2)= xr(1,1)/det
      sh(2,3)= -xr(1,1)/det + xr(2,1)/det