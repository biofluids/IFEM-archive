      vol = 0

      do iq = 1,nquad

c  Copied from sh3d8n.h
        
        xr(1,1) =
     &       + sq(1,1,iq) * xn(1, 1) + sq(1,2,iq) * xn(1, 2)
     &       + sq(1,3,iq) * xn(1, 3) + sq(1,4,iq) * xn(1, 4)
     &       + sq(1,5,iq) * xn(1, 5) + sq(1,6,iq) * xn(1, 6)
     &       + sq(1,7,iq) * xn(1, 7) + sq(1,8,iq) * xn(1, 8)
        xr(1,2) =
     &       + sq(1,1,iq) * xn(2, 1) + sq(1,2,iq) * xn(2, 2)
     &       + sq(1,3,iq) * xn(2, 3) + sq(1,4,iq) * xn(2, 4)
     &       + sq(1,5,iq) * xn(2, 5) + sq(1,6,iq) * xn(2, 6)
     &       + sq(1,7,iq) * xn(2, 7) + sq(1,8,iq) * xn(2, 8)
        xr(1,3) =
     &       + sq(1,1,iq) * xn(3, 1) + sq(1,2,iq) * xn(3, 2)
     &       + sq(1,3,iq) * xn(3, 3) + sq(1,4,iq) * xn(3, 4)
     &       + sq(1,5,iq) * xn(3, 5) + sq(1,6,iq) * xn(3, 6)
     &       + sq(1,7,iq) * xn(3, 7) + sq(1,8,iq) * xn(3, 8)
        xr(2,1) =
     &       + sq(2,1,iq) * xn(1, 1) + sq(2,2,iq) * xn(1, 2)
     &       + sq(2,3,iq) * xn(1, 3) + sq(2,4,iq) * xn(1, 4)
     &       + sq(2,5,iq) * xn(1, 5) + sq(2,6,iq) * xn(1, 6)
     &       + sq(2,7,iq) * xn(1, 7) + sq(2,8,iq) * xn(1, 8)
        xr(2,2) =
     &       + sq(2,1,iq) * xn(2, 1) + sq(2,2,iq) * xn(2, 2)
     &       + sq(2,3,iq) * xn(2, 3) + sq(2,4,iq) * xn(2, 4)
     &       + sq(2,5,iq) * xn(2, 5) + sq(2,6,iq) * xn(2, 6)
     &       + sq(2,7,iq) * xn(2, 7) + sq(2,8,iq) * xn(2, 8)
        xr(2,3) =
     &       + sq(2,1,iq) * xn(3, 1) + sq(2,2,iq) * xn(3, 2)
     &       + sq(2,3,iq) * xn(3, 3) + sq(2,4,iq) * xn(3, 4)
     &       + sq(2,5,iq) * xn(3, 5) + sq(2,6,iq) * xn(3, 6)
     &       + sq(2,7,iq) * xn(3, 7) + sq(2,8,iq) * xn(3, 8)
        xr(3,1) =
     &       + sq(3,1,iq) * xn(1, 1) + sq(3,2,iq) * xn(1, 2)
     &       + sq(3,3,iq) * xn(1, 3) + sq(3,4,iq) * xn(1, 4)
     &       + sq(3,5,iq) * xn(1, 5) + sq(3,6,iq) * xn(1, 6)
     &       + sq(3,7,iq) * xn(1, 7) + sq(3,8,iq) * xn(1, 8)
        xr(3,2) =
     &       + sq(3,1,iq) * xn(2, 1) + sq(3,2,iq) * xn(2, 2)
     &       + sq(3,3,iq) * xn(2, 3) + sq(3,4,iq) * xn(2, 4)
     &       + sq(3,5,iq) * xn(2, 5) + sq(3,6,iq) * xn(2, 6)
     &       + sq(3,7,iq) * xn(2, 7) + sq(3,8,iq) * xn(2, 8)
        xr(3,3) =
     &       + sq(3,1,iq) * xn(3, 1) + sq(3,2,iq) * xn(3, 2)
     &       + sq(3,3,iq) * xn(3, 3) + sq(3,4,iq) * xn(3, 4)
     &       + sq(3,5,iq) * xn(3, 5) + sq(3,6,iq) * xn(3, 6)
     &       + sq(3,7,iq) * xn(3, 7) + sq(3,8,iq) * xn(3, 8)

c  jacobian
        cf(1,1) = + (xr(2,2)*xr(3,3) - xr(3,2)*xr(2,3))
        cf(1,2) = - (xr(1,2)*xr(3,3) - xr(3,2)*xr(1,3))
        cf(1,3) = + (xr(1,2)*xr(2,3) - xr(2,2)*xr(1,3))
        
        det = ( xr(1,1) * cf(1,1)
     &       + xr(2,1) * cf(1,2)
     &       + xr(3,1) * cf(1,3) )

        vol = vol + abs(det) * wq(iq);

      end do
