!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine RKPMshape3d(shp,b,bd,cpt,cjt,haj,wjt)
!
!
!     This subroutine is to calculate the shape functions and 
!     their derivatives for the moving least square method
!     ---- reproducing kernel method.
!     This code is only offering shape functions and its' 1st 
!     order derivatives for 3-D case.
!
!     That is :  N_{jt}(xpt,ypt)
!     jt: the integer index for shape function, i.e.
!         we are computing N_{jt}; 1<jt<nep
!     
!     Date: June, 1994
!
!     -arguments:
!
!  i      cjt(3):  the jt th particle's global coordinations; 
!                  (jt is the integer index for shape function)
!                  i.e.  N_{jt},  1 < jt < nep
!
!                      such as x(jt):= cjt(1)
!                              y(jt):= cjt(2)
!                              z(jt):= cjt(3)
!
!  i      cpt(3): the point at where the shape function is
!                     evaluated.;
!                 x(pt):= cpt(1)
!                 y(pt):= cpt(2)
!                 z(pt):= cpt(3)
!
!
!  i      haj(3): the radius of particle's compact support.
!                 the radius at X direction: haj(1)
!                 the radius at Y direction: haj(2)
!                 the radius at Z direction: haj(3)
!
!  i      wjt: integration weight at particle jt;
!
!     o   coref: correct function at pt. cpt and particle cjt
!         C(cjt, cpt);
!
!
!     o   shp: shape functions and their 1st derivative
!                       at point cpt;
!
!     o   shpd(1): the X-derivative of shape function
!                          N_{jt},x  at cpt;  
!
!     o   shpd(2): the Y-derivative of shape function
!                          N_{jt},y at cpt;  
!    
!     o   shpd(3): the Z-derivative of shape function
!                          N_{jt},z at cpt;  
!
!     The subroutine is only designed for solving 3-D locally
!     linear interpolation moving field problem. i.e.
!
!    l   p0   := (1,0,0,0); 
!
!    l   p(4) := (1,x,y,z);
!
!    There two different subroutines are called by this subroutine:
!
!    window2dr.f ----- the cartesian product type window function;
!  
!    window2dc.f ----- the circular disc type window function
!
!**************************************************************
!
      implicit none
      real(8) :: cpt(3),cjt(3),shp,shpd(3)
      real(8) :: p(4),pdx(4),pdy(4),pdz(4)
      real(8) :: b(4),bd(3,4),haj(3)
      real(8) :: coref,cdx,cdy,cdz
      real(8) :: xp,yp,zp,xj,yj,zj
      real(8) :: a1,a2,a3
      real(8) :: aw,awdx,awdy,awdz
      real(8) :: awdxx,awdyy,awdzz
      real(8) :: awdxy,awdyz,awdzx
      real(8) :: wjt,xx,yy,zz
      integer :: i,nn
!
!.....find the polynominal p
!
      nn = 4
      xp = cpt(1)
      yp = cpt(2)
      zp = cpt(3)

      xj = cjt(1)
      yj = cjt(2)
      zj = cjt(3)
!
!.....weight
!

!.....assign value to working array
!
      a1 = haj(1)
      a2 = haj(2)
      a3 = haj(3)

      p(1) = 1.00
      p(2) = (xj - xp)/a1
      p(3) = (yj - yp)/a2
      p(4) = (zj - zp)/a3

      pdx(1) =   0.00
      pdx(2) = - 1.00/a1
      pdx(3) =   0.00
      pdx(4) =   0.00

      pdy(1) =   0.00
      pdy(2) =   0.00
      pdy(3) = - 1.00/a2
      pdy(4) =   0.00

      pdz(1) =   0.00
      pdz(2) =   0.00
      pdz(3) =   0.00
      pdz(4) = - 1.00/a3
!
!.....compute correct function  C
!
      coref = 0.00
      do i = 1, nn
       coref = coref + p(i)*b(i)
      enddo


      cdx = 0.00
      cdy = 0.00
      cdz = 0.00
      do i = 1, nn
       cdx = cdx + pdx(i)*b(i) + p(i)*bd(1,i)
       cdy = cdy + pdy(i)*b(i) + p(i)*bd(2,i)
       cdz = cdz + pdz(i)*b(i) + p(i)*bd(3,i)
      enddo

!
!.....check whether or not the point is within the compact support
!
      xx = abs(xj - xp)/a1
      yy = abs(yj - yp)/a2
      zz = abs(zj - zp)/a3

      if (xx .ge. 2.0 .or. yy .ge. 2.0 .or. zz .ge. 2.0) then
        shp      = 0.00
        shpd(1)  = 0.00
        shpd(2)  = 0.00
        shpd(3)  = 0.00
      else
          call window3d(aw,awdx,awdy,awdz,  &
              awdxx,awdyy,awdzz,            &
              awdxy,awdyz,awdzx,            &
              xp,yp,zp,xj,yj,zj,a1,a2,a3)
!
!.....calculate the shape functions
!
          shp   = coref*aw*wjt
!
!.....calculate the 1st derivatives of the shape function
!
          shpd(1) = (cdx*aw+coref*awdx)*wjt
          shpd(2) = (cdy*aw+coref*awdy)*wjt
          shpd(3) = (cdz*aw+coref*awdz)*wjt

      endif


  999 continue

      return
      end

!***********************************************************
      subroutine window3d(aw,awdx,awdy,awdz,  &
                          awdxx,awdyy,awdzz,  &
                          awdxy,awdyz,awdzx,  &
                  xp,yp,zp,xj,yj,zj,a1,a2,a3)
!***********************************************************
!
!   This is a subroutine to compute 3-D cubic spline window
!   function and its first and second derivatives.
!   In this subroutine, the window function is constructed by 
!   Cartesian product.
!   The mathematical formulation of this code is based on:
!
!   ``Ten lectures on Wavelets'' by Ingrid Daubechies, 
!     page 79.
!
!   arguments:
!
!   i         a1:  dilation parameter in X-direction;
!
!   i         a2:  dilation parameter in Y-direction;
!
!   i         a3:  dilation parameter in Z-direction;
!
!   Note: there is a difference between dilation parameter and
!         the dilation coefficient.
!
!      o      aw:  3-D cubic spline window function:
!
!   Note:     Phi_r := (1/rho) Phi(X/rho).
!
!      o    awdx:  d/dx(aw);
!
!      o    awdy:  d/dy(aw);
!
!      o    awdz:  d/dz(aw);
!
!      o   awdxx:  d^2/dx^2(aw);
!
!      o   awdyy:  d^2/dy^2(aw);
!
!      o   awdzz:  d^2/dz^2(aw);
!
!      o   awdxy:  d^2/dxdy(aw);
!
!      o   awdyz:  d^2/dydz(aw);
!
!      o   awdzx:  d^2/dzdx(aw);
!
!   i         xp:  X-coordinate for arbitary point;
!
!   i         yp:  Y-coordinate for arbitary point;
!
!   i         zp:  Z-coordinate for arbitary point;
!
!   i         xj:  X-coordinate for particle j;
!
!   i         yj:  Y-coordinate for particle j;
!
!   i         zj:  Y-coordinate for particle j;
!
!.....cartsian product window function
!
!    Author: Shaofan Li
!
!    Date: October, 1996
!   
!
!*************************************************************
!
      implicit none
!
!.....global variables
!
      real(8) aw,awdx,awdy,awdz
      real(8) awdxx,awdyy,awdzz
      real(8) awdxy,awdyz,awdzx
      real(8) a1,a2,a3
      real(8) xj,yj,zj,xp,yp,zp
!
!.....Local variables
!
      real(8) hv,zero
      real(8) x1,x2,x3,y1,y2,y3,z1,z2,z3
      real(8) dx1,dx2,dy1,dy2,dz1,dz2
      real(8) onen,onep,twon,twop
      real(8) awx,awxd,awxdd
      real(8) awy,awyd,awydd
      real(8) awz,awzd,awzdd
!
!.....normalize the argument
!
      x1  = (xj-xp)/a1
      x2  =  x1*x1
      x3  =  x2*x1

      y1  = (yj-yp)/a2
      y2  =  y1 * y1
      y3  =  y2 * y1

      z1  = (zj-zp)/a3
      z2  =  z1 * z1
      z3  =  z2 * z1

      twon = -2.00
      onen = -1.00
      zero =  0.00
      onep =  1.00
      twop =  2.00
!
!.....dfx := d(xr)/dx; dfy:= d(yr)/dy
!
      dx1  = -1.00/a1
      dx2 = dx1 * dx1
      dy1  = -1.00/a2
      dy2 = dy1 * dy1
      dz1  = -1.00/a3
      dz2 = dz1 * dz1

      hv = 1.0/(a1*a2*a3)

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

      return
      end
