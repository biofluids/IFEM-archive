      subroutine window(aw,awdx,awdy,awdxx,awdxy,
     &                  awdyy,xx,yy,xj,yj,a1,a2)
c***********************************************************
c
c   This is a subroutine to compute 2-D cubic spline window
c   function and its first and second derivatives.
c   In this subroutine, the window function is constructed by 
c   Cartesian product.
c   The mathematical formulation of this code is based on:
c
c   ``Ten lectures on Wavelets'' by Ingrid Daubechies, 
c     page 79.
c
c   arguments:
c
c   i         a1:  dilation parameter in X-direction;
c
c   i         a2:  dilation parameter in Y-direction;
c
c   Note: there is a difference between dilation parameter and
c         the dilation coefficient.
c
c      o      aw:  2-D cubic spline window function:
c
c   Note:     Phi_r := (1/rho) Phi(X/rho).
c
c      o    awdx:  d/dx(aw);
c
c      o    awdy:  d/dy(aw);
c
c      o   awdxx:  d^2/dx^2(aw);
c
c      o   awdyy:  d^2/dy^2(aw);
c
c      o   awdxy:  d^2/dxdy(aw);
c
c   i         xx:  X-coordinate for arbitary point;
c
c   i         yy:  Y-coordinate for arbitary point;
c
c   i         xj:  X-coordinate for particle j;
c
c   i         yj:  Y-coordinate for particle j;
c
c
c    Author: Shaofan Li
c
c    Date: October, 1996
c   
c
c*************************************************************
      implicit double precision (a-h,o-z)
c
c.....Cartsian product window function
c
c.....normalize the argument
c
      x1  = (xj-xx)/a1
      x2  =  x1 * x1
      x3  =  x2 * x1
c
      y1  = (yj-yy)/a2
      y2  =  y1 * y1
      y3  =  y2 * y1
c
      two1 = -2.0d0
      one1 = -1.0d0
      zero =  0.0d0
      one2 =  1.0d0
      two2 =  2.0d0
c
c.....dfx := d(xr)/dx; dfy:= d(yr)/dy
c
      dx1  = -1.0d0/a1
      dx2  =  dx1 * dx1
      dy1  = -1.0d0/a2
      dy2  =  dy1 * dy1
c
      hv = 1.0d0/(a1*a2)
c
      if((x1.ge.two1).and.(x1.lt.one1)) then
        awx   = (1.0d0/6.0d0)*(2.0d0 + x1)**3.
        awxd  = dx1*0.5d0*(2.0d0 + x1)**2.0
        awxdd = dx2*(2.0d0 + x1)
      elseif ((x1.ge.one1).and.(x1.lt.zero)) then
        awx   =   2.0d0/3.0d0 - x2 - 0.50d0*x3
        awxd  = - dx1*(2.0d0*x1 + 1.50d0*x2)
        awxdd = - dx2*(2.0d0 + 3.0d0*x1)
      elseif ((x1.ge.zero).and.(x1.lt.one2)) then
        awx   =   2.0d0/3.0d0 - x2 + 0.50d0*x3
        awxd  = - dx1*(2.0d0*x1 - 1.50d0*x2)
        awxdd = - dx2*(2.0d0 - 3.0d0*x1 )
      elseif ((x1.ge.one2).and.(x1.le.two2)) then
        awx   =   (1.0d0/6.0d0)*(2.0d0 - x1)**3.
        awxd  = - dx1*0.5d0*(2.0d0 - x1)**2.0
        awxdd =   dx2*(2.0d0 - x1)
      else
        awx   = 0.00d0
        awxd  = 0.00d0
        awxdd = 0.00d0
      endif
c     
c
      if((y1.ge.two1).and.(y1.lt.one1)) then
        awy   = (1.0d0/6.0d0)*(2.0d0 + y1)**3.0
        awyd  = dy1*0.5d0*(2.0d0 + y1 )**2.
        awydd = dy2*(2.0d0 + y1 )
      elseif ((y1.ge.one1).and.(y1.lt.zero)) then
        awy   =   2.0d0/3.0d0 - y2 - 0.5d0*y3
        awyd  = - dy1*(2.0d0*y1 + 1.5d0*y2)
        awydd = - dy2*(2.0d0 + 3.0d0*y1 )
      elseif ((y1.ge.zero).and.(y1.lt.one2)) then
        awy   =   2.0d0/3.0d0 - y2 + 0.5d0*y3
        awyd  = - dy1*(2.0d0*y1 - 1.5d0*y2)
        awydd = - dy2*(2.0d0 - 3.0d0*y1 )
      elseif ((y1.ge.one2).and.(y1.le.two2)) then
        awy   =   (1.0d0/6.0d0)*(2.0d0 - y1)**3.0
        awyd  = - dy1*0.5d0*(2.0d0 - y1)**2.0
        awydd =   dy2*(2.0d0 - y1)
      else
        awy   = 0.00d0
        awyd  = 0.00d0 
        awydd = 0.00d0
      endif
c     
c
      aw    = awx*awy*hv
      awdx  = awxd*awy*hv
      awdy  = awx*awyd*hv
      awdxy = awxd*awyd*hv
      awdxx = awxdd*awy*hv
      awdyy = awx*awydd*hv
c
      return
      end
c
