      subroutine RKPMshape2d(shp,b,bd,cpt,cjt,haj,wjt)
!
!**************************************************************
!
!     This subroutine is to calculate the shape functions and 
!     their derivatives for the moving least square method
!     ---- reproducing kernel method.
!     This code is only offering shape functions and its' 1st 
!     order derivatives for 2-D case.
!
!     That is :  N_{jt}(xpt,ypt)
!     jt: the integer index for shape function, i.e.
!         we are computing N_{jt}; 1<jt<nep
!     
!     Date: June, 1994
!
!     -arguments:
!
!  i      cjt(2):  the jt th particle's global coordinations; 
!                  (jt is the integer index for shape function)
!                  i.e.  N_{jt},  1 < jt < nep
!
!                      such as x(jt):= cjt(1)
!                              y(jt):= cjt(2)
!
!
!  i      cpt(2): the point at where the shape function is
!                     evaluated.;
!                 x(pt):= cpt(1)
!                 y(pt):= cpt(2)
!
!
!  i      haj(2): the radius of particle's compact support.
!                 the radius at X direction: haj(1)
!                 the radius at Y direction: haj(2)
!
!
!  i      wjt: integration weight at particle jt;
!
!
!     o   coref: correct function at pt. cpt and particle cjt
!         C(cjt, cpt);
!
!
!     o   shp: shape functions and their 1st derivative
!                       at point cpt;
!
!     o   shp(0)    : shape function N_{jt} at cpt ;
!
!     o   shp(1): the X-derivative of shape function
!                          N_{jt},x  at cpt;  
!
!     o   shp(2): the Y-derivative of shape function
!                          N_{jt},y at cpt;  
!    
!
!     The subroutine is only designed for solving 2-D locally
!     linear interpolation moving field problem. i.e.
!
!    l   p0   := (1,0,0); 
!
!    l   p(3) := (1,x,y);
!
!    There two different subroutines are called by this subroutine:
!
!    window2dr.f ----- the cartesian product type window function;
!  
!    window2dc.f ----- the circular disc type window function
!
!**************************************************************
      implicit none
!
      real*8 cpt(2),shpd(2)
      real*8 cjt(2),p(3),pdx(3),pdy(3)
      real*8 b(*),bd(2,*),haj(2)
      real*8 coref,cdx,cdy,ha1,ha2,xx,yy,xj,yj,xp,yp,shp,wjt
      real*8 aw,awdx,awdy,awdxx,awdxy,awdyy
      integer i
!
!.....find the polynominal p
!
      xp = cpt(1)
      yp = cpt(2)
!
      xj = cjt(1)
      yj = cjt(2)
!
!.....assign value to working array
!
      ha1    = haj(1)
      ha2    = haj(2)
!
      p(1)   = 1.0d0
      p(2)   = (xj - xp)/ha1
      p(3)   = (yj - yp)/ha2
!
      pdx(1) =   0.0d0
      pdx(2) = - 1.0d0/ha1
      pdx(3) =   0.0d0
!
      pdy(1) =   0.0d0
      pdy(2) =   0.0d0
      pdy(3) = - 1.0d0/ha2
!
!.....compute correct function  C
!
      coref = 0.0d0
      do i = 1,3
		coref = coref + p(i)*b(i)
      enddo
!
      cdx = 0.0d0
      cdy = 0.0d0
      do i = 1,3
		 cdx = cdx + pdx(i)*b(i) + p(i)*bd(1,i)
		 cdy = cdy + pdy(i)*b(i) + p(i)*bd(2,i)
      enddo

!
!.....check whether or not the point is within the compact support
!
      xx = dabs(xj - xp)/ha1
      yy = dabs(yj - yp)/ha2

      if (xx .gt. 2.0 .or. yy .gt. 2.0) then
		  shp      = 0.00
		  shpd(1)  = 0.00
		  shpd(2)  = 0.00
      else
          call window(aw,awdx,awdy,awdxx,awdxy,awdyy,xp,yp,xj,yj,ha1,ha2)
!
!.....calculate the shape functions
!
          shp     = coref*aw*wjt
!
!.....calculate the 1st derivatives of the shape function
!
          shpd(1) = (cdx*aw+coref*awdx)*wjt
          shpd(2) = (cdy*aw+coref*awdy)*wjt
!
      end if
!
  999 continue
!
      return
      end
!



      subroutine window(aw,awdx,awdy,awdxx,awdxy,awdyy,xx,yy,xj,yj,a1,a2)
!***********************************************************
!
!   This is a subroutine to compute 2-D cubic spline window
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
!   Note: there is a difference between dilation parameter and
!         the dilation coefficient.
!
!      o      aw:  2-D cubic spline window function:
!
!   Note:     Phi_r := (1/rho) Phi(X/rho).
!
!      o    awdx:  d/dx(aw);
!
!      o    awdy:  d/dy(aw);
!
!      o   awdxx:  d^2/dx^2(aw);
!
!      o   awdyy:  d^2/dy^2(aw);
!
!      o   awdxy:  d^2/dxdy(aw);
!
!   i         xx:  X-coordinate for arbitary point;
!
!   i         yy:  Y-coordinate for arbitary point;
!
!   i         xj:  X-coordinate for particle j;
!
!   i         yj:  Y-coordinate for particle j;
!
!
!    Author: Shaofan Li
!
!    Date: October, 1996
!   
!
!*************************************************************
      implicit double precision (a-h,o-z)
!
!.....Cartsian product window function
!
!.....normalize the argument
!
      x1  = (xj-xx)/a1
      x2  =  x1 * x1
      x3  =  x2 * x1
!
      y1  = (yj-yy)/a2
      y2  =  y1 * y1
      y3  =  y2 * y1
!
      two1 = -2.0d0
      one1 = -1.0d0
      zero =  0.0d0
      one2 =  1.0d0
      two2 =  2.0d0
!
!.....dfx := d(xr)/dx; dfy:= d(yr)/dy
!
      dx1  = -1.0d0/a1
      dx2  =  dx1 * dx1
      dy1  = -1.0d0/a2
      dy2  =  dy1 * dy1
!
      hv = 1.0d0/(a1*a2)
!
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

      aw    = awx*awy*hv
      awdx  = awxd*awy*hv
      awdy  = awx*awyd*hv
      awdxy = awxd*awyd*hv
      awdxx = awxdd*awy*hv
      awdyy = awx*awydd*hv

      return
      end

