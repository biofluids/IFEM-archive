      subroutine RKPMshape_La(shp,shpd,b,bd,cpt,cjt,haj,wjt)
c
c**************************************************************
c
c     This subroutine is to calculate the shape functions and 
C     their derivatives for the moving least square method
c     ---- reproducing kernel method.
c     This code is only offering shape functions and its' 1st 
c     order derivatives for 2-D case.
c
c     That is :  N_{jt}(xpt,ypt)
c     jt: the integer index for shape function, i.e.
c         we are computing N_{jt}; 1<jt<nep
c     
c     Date: June, 1994
c
c     -arguments:
c
c  i      cjt(2):  the jt th particle's global coordinations; 
c                  (jt is the integer index for shape function)
c                  i.e.  N_{jt},  1 < jt < nep
c
c                      such as x(jt):= cjt(1)
c                              y(jt):= cjt(2)
c
c
c  i      cpt(2): the point at where the shape function is
c                     evaluated.;
c                 x(pt):= cpt(1)
c                 y(pt):= cpt(2)
c
c
c  i      haj(2): the radius of particle's compact support.
c                 the radius at X direction: haj(1)
c                 the radius at Y direction: haj(2)
c
c
c  i      wjt: integration weight at particle jt;
c
c
c     o   coref: correct function at pt. cpt and particle cjt
c         C(cjt, cpt);
c
c
c     o   shp: shape functions and their 1st derivative
c                       at point cpt;
c
c     o   shp(0)    : shape function N_{jt} at cpt ;
c
c     o   shp(1): the X-derivative of shape function
c                          N_{jt},x  at cpt;  
c
c     o   shp(2): the Y-derivative of shape function
c                          N_{jt},y at cpt;  
c    
c
c     The subroutine is only designed for solving 2-D locally
c     linear interpolation moving field problem. i.e.
c
c    l   p0   := (1,0,0); 
c
c    l   p(3) := (1,x,y);
c
c    There two different subroutines are called by this subroutine:
c
c    window2dr.f ----- the cartesian product type window function;
c  
c    window2dc.f ----- the circular disc type window function
c
c**************************************************************
      implicit none
c
      real*8 cpt(2),shpd(2)
      real*8 cjt(2),p(3),pdx(3),pdy(3)
      real*8 b(*),bd(2,*),haj(2)
      real*8 coref,cdx,cdy,ha1,ha2,xx,yy,xj,yj,
     &       xp,yp,shp,wjt,
     &       aw,awdx,awdy,awdxx,awdxy,awdyy
      integer i
c
c.....find the polynominal p
c
      xp = cpt(1)
      yp = cpt(2)
c
      xj = cjt(1)
      yj = cjt(2)
c
c.....assign value to working array
c
      ha1    = haj(1)
      ha2    = haj(2)
c
      p(1)   = 1.0d0
      p(2)   = (xj - xp)/ha1
      p(3)   = (yj - yp)/ha2
c
      pdx(1) =   0.0d0
      pdx(2) = - 1.0d0/ha1
      pdx(3) =   0.0d0
c
      pdy(1) =   0.0d0
      pdy(2) =   0.0d0
      pdy(3) = - 1.0d0/ha2
c
c.....compute correct function  C
c
      coref = 0.0d0
      do i = 1,3
	 coref = coref + p(i)*b(i)
      enddo
c
      cdx = 0.0d0
      cdy = 0.0d0
      do i = 1,3
	 cdx = cdx + pdx(i)*b(i) + p(i)*bd(1,i)
	 cdy = cdy + pdy(i)*b(i) + p(i)*bd(2,i)
      enddo
c
c
c.....check whether or not the point is within the compact support
c
      xx = dabs(xj - xp)/ha1
      yy = dabs(yj - yp)/ha2

      if (xx .gt. 2.0 .or. yy .gt. 2.0) then
	  shp      = 0.00
	  shpd(1)  = 0.00
	  shpd(2)  = 0.00
      else
          call window(aw,awdx,awdy,
     &           awdxx,awdxy,awdyy,
     &         xp,yp,xj,yj,ha1,ha2)
c
c.....calculate the shape functions
c
          shp     = coref*aw*wjt
c
c.....calculate the 1st derivatives of the shape function
c
          shpd(1) = (cdx*aw+coref*awdx)*wjt
          shpd(2) = (cdy*aw+coref*awdy)*wjt
c
      end if
c
  999 continue
c
      return
      end
c
