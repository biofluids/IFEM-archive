c  Copied from sh3d4n.h
      xr(1,1) = xn(1,1) - xn(1,4)
      xr(1,2) = xn(2,1) - xn(2,4)
      xr(1,3) = xn(3,1) - xn(3,4)
      xr(2,1) = xn(1,2) - xn(1,4)
      xr(2,2) = xn(2,2) - xn(2,4)
      xr(2,3) = xn(3,2) - xn(3,4)
      xr(3,1) = xn(1,3) - xn(1,4)
      xr(3,2) = xn(2,3) - xn(2,4)
      xr(3,3) = xn(3,3) - xn(3,4)
c  jacobian
      cf(1,1) = + (xr(2,2)*xr(3,3) - xr(3,2)*xr(2,3))
      cf(1,2) = - (xr(1,2)*xr(3,3) - xr(3,2)*xr(1,3))
      cf(1,3) = + (xr(1,2)*xr(2,3) - xr(2,2)*xr(1,3))
      det = ( xr(1,1) * cf(1,1)
     &     + xr(2,1) * cf(1,2)
     &     + xr(3,1) * cf(1,3) )
      vol = abs(det)/6.0

     
