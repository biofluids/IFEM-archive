      subroutine window3d(aw,awdx,awdy,awdz,
     &                    awdxx,awdyy,awdzz,
     &                    awdxy,awdyz,awdzx,
     &            xp,yp,zp,xj,yj,zj,a1,a2,a3)
c***********************************************************
c
c   This is a subroutine to compute 3-D cubic spline window
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
c   i         a3:  dilation parameter in Z-direction;
c
c   Note: there is a difference between dilation parameter and
c         the dilation coefficient.
c
c      o      aw:  3-D cubic spline window function:
c
c   Note:     Phi_r := (1/rho) Phi(X/rho).
c
c      o    awdx:  d/dx(aw);
c
c      o    awdy:  d/dy(aw);
c
c      o    awdz:  d/dz(aw);
c
c      o   awdxx:  d^2/dx^2(aw);
c
c      o   awdyy:  d^2/dy^2(aw);
c
c      o   awdzz:  d^2/dz^2(aw);
c
c      o   awdxy:  d^2/dxdy(aw);
c
c      o   awdyz:  d^2/dydz(aw);
c
c      o   awdzx:  d^2/dzdx(aw);
c
c   i         xp:  X-coordinate for arbitary point;
c
c   i         yp:  Y-coordinate for arbitary point;
c
c   i         zp:  Z-coordinate for arbitary point;
c
c   i         xj:  X-coordinate for particle j;
c
c   i         yj:  Y-coordinate for particle j;
c
c   i         zj:  Y-coordinate for particle j;
c
c.....cartsian product window function
c
c    Author: Shaofan Li
c
c    Date: October, 1996
c   
c
c*************************************************************
c
      implicit none
c
c.....global variables
c
      real*8 aw,awdx,awdy,awdz
      real*8 awdxx,awdyy,awdzz
      real*8 awdxy,awdyz,awdzx
      real*8 a1,a2,a3
      real*8 xj,yj,zj,xp,yp,zp
c
c.....Local variables
c
      real*8 hv,zero
      real*8 x1,x2,x3,y1,y2,y3,z1,z2,z3
      real*8 dx1,dx2,dy1,dy2,dz1,dz2
      real*8 onen,onep,twon,twop
      real*8 awx,awxd,awxdd
      real*8 awy,awyd,awydd
      real*8 awz,awzd,awzdd
c
c.....normalize the argument
c
      x1  = (xj-xp)/a1
      x2  =  x1*x1
      x3  =  x2*x1
c
      y1  = (yj-yp)/a2
      y2  =  y1 * y1
      y3  =  y2 * y1
c
      z1  = (zj-zp)/a3
      z2  =  z1 * z1
      z3  =  z2 * z1
c
      twon = -2.00
      onen = -1.00
      zero =  0.00
      onep =  1.00
      twop =  2.00
c
c.....dfx := d(xr)/dx; dfy:= d(yr)/dy
c
      dx1  = -1.00/a1
      dx2 = dx1 * dx1
      dy1  = -1.00/a2
      dy2 = dy1 * dy1
      dz1  = -1.00/a3
      dz2 = dz1 * dz1
c
      hv = 1.0/(a1*a2*a3)
c
      if((x1.ge.twon).and.(x1.lt.onen)) then
        awx   = (1.0/6.0)*(2.0 + x1)**3.
        awxd  = dx1*0.5*(2.0 + x1)**2.
        awxdd = dx2*(2.0 + x1)
      elseif ((x1.ge.onen).and.(x1.lt.zero)) then
        awx   = 2.0/3.0 - x2 - 0.5*x3
        awxd  = - dx1*(2.0*x1 + 1.5*x2)
        awxdd = - dx2*(2.0 + 3.0*x1 )
      elseif ((x1.ge.zero).and.(x1.lt.onep)) then
        awx   = 2.0/3.0 - x2 + 0.5*x3
        awxd  = - dx1*(2.0*x1 - 1.5*x2)
        awxdd = - dx2*(2.0 - 3.0*x1 )
      elseif ((x1.ge.onep).and.(x1.le.twop)) then
        awx   = (1.0/6.0)*(2.0 - x1)**3.
        awxd  = - dx1*0.5*(2.0 - x1)**2.
        awxdd = dx2*(2.0 - x1)
      else
        awx   = 0.00
        awxd  = 0.00 
        awxdd = 0.00
      endif
c     
c
      if((y1.ge.twon).and.(y1.lt.onen)) then
        awy   = (1.0/6.0)*(2.0 + y1)**3.
        awyd  = dy1*0.5*(2.0 + y1 )**2.
        awydd = dy2*(2.0 + y1 )
      elseif ((y1.ge.onen).and.(y1.lt.zero)) then
        awy   = 2.0/3.0 - y2 - 0.5*y3
        awyd  = - dy1*(2.0*y1 + 1.5*y2)
        awydd = - dy2*(2.0 + 3.0*y1 )
      elseif ((y1.ge.zero).and.(y1.lt.onep)) then
        awy   = 2.0/3.0 - y2 + 0.5*y3
        awyd  = - dy1*(2.0*y1 - 1.5*y2)
        awydd = - dy2*(2.0 - 3.0*y1 )
      elseif ((y1.ge.onep).and.(y1.le.twop)) then
        awy   = (1.0/6.0)*(2.0 - y1)**3.
        awyd  = - dy1*0.5*(2.0 - y1)**2.
        awydd = dy2*(2.0 - y1)
      else
        awy   = 0.00
        awyd  = 0.00 
        awydd = 0.00
      endif
c     
c
      if((z1.ge.twon).and.(z1.lt.onen)) then
        awz   = (1.0/6.0)*(2.0 + z1)**3.
        awzd  = dz1*0.5*(2.0 + z1 )**2.
        awzdd = dz2*(2.0 + z1 )
      elseif ((z1.ge.onen).and.(z1.lt.zero)) then
        awz   = 2.0/3.0 - z2 - 0.5*z3
        awzd  = - dz1*(2.0*z1 + 1.5*z2)
        awzdd = - dz2*(2.0 + 3.0*z1 )
      elseif ((z1.ge.zero).and.(z1.lt.onep)) then
        awz   = 2.0/3.0 - z2 + 0.5*z3
        awzd  = - dz1*(2.0*z1 - 1.5*z2)
        awzdd = - dz2*(2.0 - 3.0*z1 )
      elseif ((z1.ge.onep).and.(z1.le.twop)) then
        awz   = (1.0/6.0)*(2.0 - z1)**3.
        awzd  = - dz1*0.5*(2.0 - z1)**2.
        awzdd = dz2*(2.0 - z1)
      else
        awz   = 0.00
        awzd  = 0.00 
        awzdd = 0.00
      endif
c     
c
      aw    = awx*awy*awz*hv
      awdx  = awxd*awy*awz*hv
      awdy  = awx*awyd*awz*hv
      awdz  = awx*awy*awzd*hv
      awdxy = awxd*awyd*awz*hv
      awdyz = awx*awyd*awzd*hv
      awdzx = awxd*awy*awzd*hv
      awdxx = awxdd*awy*awz*hv
      awdyy = awx*awydd*awz*hv
      awdzz = awx*awy*awzdd*hv
c
      return
      end
c
