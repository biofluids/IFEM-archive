      subroutine RKPMshape3dtl(shp,shpd,b,bd,cpt,cjt,haj,wjt)
c
c**************************************************************
c
c     This subroutine is to calculate the shape functions and 
C     their derivatives for the moving least square method
c     ---- reproducing kernel method.
c     This code is only offering shape functions and its' 1st 
c     order derivatives for 3-D case.
c
c     That is :  N_{jt}(xpt,ypt)
c     jt: the integer index for shape function, i.e.
c         we are computing N_{jt}; 1<jt<nep
c     
c     Date: June, 1994
c
c     -arguments:
c
c  i      cjt(3):  the jt th particle's global coordinations; 
c                  (jt is the integer index for shape function)
c                  i.e.  N_{jt},  1 < jt < nep
c
c                      such as x(jt):= cjt(1)
c                              y(jt):= cjt(2)
c                              z(jt):= cjt(3)
c
c
c  i      cpt(3): the point at where the shape function is
c                     evaluated.;
c                 x(pt):= cpt(1)
c                 y(pt):= cpt(2)
c                 z(pt):= cpt(3)
c
c
c  i      haj(3): the radius of particle's compact support.
c                 the radius at X direction: haj(1)
c                 the radius at Y direction: haj(2)
c                 the radius at Z direction: haj(3)
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
c     o   shpd(1): the X-derivative of shape function
c                          N_{jt},x  at cpt;  
c
c     o   shpd(2): the Y-derivative of shape function
c                          N_{jt},y at cpt;  
c    
c     o   shpd(3): the Z-derivative of shape function
c                          N_{jt},z at cpt;  
c
c     The subroutine is only designed for solving 3-D locally
c     linear interpolation moving field problem. i.e.
c
c    l   p0   := (1,0,0,0,0,0,0); 
c
c    l   p(*) := (1,x,y,z,xy,yz,zx);
c
c    There two different subroutines are called by this subroutine:
c
c    window2dr.f ----- the cartesian product type window function;
c  
c    window2dc.f ----- the circular disc type window function
c
c**************************************************************
c
      implicit none
c
      real*8 cpt(3),cjt(3),shp,shpd(3)
      real*8 p(7),pdx(7),pdy(7),pdz(7)
      real*8 b(*),bd(3,*),haj(3)
      real*8 coref,cdx,cdy,cdz
      real*8 xp,yp,zp,xj,yj,zj
      real*8 a1,a2,a3
      real*8 aw,awdx,awdy,awdz
      real*8 awdxx,awdyy,awdzz
      real*8 awdxy,awdyz,awdzx
      real*8 wjt,xx,yy,zz
      integer i,nn
c
c.....find the polynominal p
c
      nn = 7
      xp = cpt(1)
      yp = cpt(2)
      zp = cpt(3)
c
      xj = cjt(1)
      yj = cjt(2)
      zj = cjt(3)
c
c.....weight
c
      
c.....assign value to working array
c
      a1 = haj(1)
      a2 = haj(2)
      a3 = haj(3)
c
      p(1) = 1.00
      p(2) = (xj - xp)/a1
      p(3) = (yj - yp)/a2
      p(4) = (zj - zp)/a3
      p(5) = (xj - xp)*(yj - yp)/(a1*a2)
      p(6) = (yj - yp)*(zj - zp)/(a2*a3)
      p(7) = (zj - zp)*(xj - xp)/(a3*a1)
c
      pdx(1) =   0.00
      pdx(2) = - 1.00/a1
      pdx(3) =   0.00
      pdx(4) =   0.00
      pdx(5) = - (yj - yp)/(a1*a2)
      pdx(6) =   0.00
      pdx(7) = - (zj - zp)/(a1*a3)
c
      pdy(1) =   0.00
      pdy(2) =   0.00
      pdy(3) = - 1.00/a2
      pdy(4) =   0.00
      pdy(5) = - (xj - xp)/(a1*a2)
      pdy(6) = - (zj - zp)/(a2*a3)
      pdy(7) =   0.00
c
      pdz(1) =   0.00
      pdz(2) =   0.00
      pdz(3) =   0.00
      pdz(4) = - 1.00/a3
      pdz(5) =   0.00
      pdz(6) = - (yj - yp)/(a2*a3)
      pdz(7) = - (xj - xp)/(a1*a3)
c
c.....compute correct function  C
c
      coref = 0.00
      do i = 1, nn
	 coref = coref + p(i)*b(i)
      enddo
c
c
      cdx = 0.00
      cdy = 0.00
      cdz = 0.00
      do i = 1, nn
	 cdx = cdx + pdx(i)*b(i) + p(i)*bd(1,i)
	 cdy = cdy + pdy(i)*b(i) + p(i)*bd(2,i)
	 cdz = cdz + pdz(i)*b(i) + p(i)*bd(3,i)
      enddo
c
c
c.....check whether or not the point is within the compact support
c
      xx = dabs(xj - xp)/a1
      yy = dabs(yj - yp)/a2
      zz = dabs(zj - zp)/a3

      if (xx .ge. 2.0 .or. yy .ge. 2.0 
     &    .or. zz .ge. 2.0) then
	  shp      = 0.00
	  shpd(1)  = 0.00
	  shpd(2)  = 0.00
	  shpd(3)  = 0.00
      else
          call window3d(aw,awdx,awdy,awdz,
     &         awdxx,awdyy,awdzz,
     &         awdxy,awdyz,awdzx,
     &         xp,yp,zp,xj,yj,zj,a1,a2,a3)
c
c.....calculate the shape functions
c
          shp   = coref*aw*wjt
c
c.....calculate the 1st derivatives of the shape function
c
          shpd(1) = (cdx*aw+coref*awdx)*wjt
          shpd(2) = (cdy*aw+coref*awdy)*wjt
          shpd(3) = (cdz*aw+coref*awdz)*wjt
c
      end if
c
c
  999 continue
c
      return
      end
c
