      subroutine RKPMshape_qd(shp,shpd,shpdd,
     &                    b,bd,bdd,cpt,cjt,anode,wjt)
c
c**************************************************************
c
c     This subroutine is to calculate the shape functions and 
c     their derivatives for moving least square reproducing kernel
c     shape function.
c     This code is only offering shape functions and its' 1st 
c     order derivatives for 2-D case.
c
c     That is :  N_{jt}(xpt,ypt)
c     jt: the integer index for shape function, i.e.
c         we are computing N_{jt}; 1<jt<nep
c     
c     Date: October, 1996
c
c     Author: Shaofan Li
c
c     Copyright @ Northwestern University
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
c  i      dcjt(2): the dilation parameter vector of the jth particle: 
c                      the dilation parameter at X direction: dcjt(1)
c                      the dilation parameter at Y direction: dcjt(2)
c
c         for irregular mesh, this matter should be re-discussed;
c
c  i      cpt(2): the point at where the shape function is
c                     evaluated.;
c                 x(pt):= cpt(1)
c                 y(pt):= cpt(2)
c
c  i      wjt: integration weight at jt;
c
c     o   shp: shape functions and their 1st derivative
c                       at point cpt;
c
c     o   shp    : shape function N_{jt} at cpt ;
c
c     o   shpd(1): the X-derivative of shape function
c                          N_{jt},x  at cpt;  
c
c     o   shpd(2): the Y-derivative of shape function
c                          N_{jt},y at cpt;  
c
c     o   shpdd(1): N_{jt},xx  at cpt;  
c
c     o   shpdd(2): N_{jt},xy  at cpt;  
c
c     o   shpdd(3): N_{jt},yy  at cpt;  
c
c    
c    l   a1      : dilation parameter; scalar;
c                  for uniformal mesh, we choose a1:= dcjt(1,ip)
c                                                a2:= dcjt(2,ip)
c
c     The subroutine is only designed for solving 2-D locally
c     linear interpolation moving field problem. i.e.
c
c    l   p0   := (1,0,0,0,0,0); 
c
c    l   p(6) := (1,x,y,x*x,x*y,y*y);
c
c    There two different subroutines are called by this subroutine:
c
c    cubic2d.f ----- the cubic spline function.
c  
c    quinwt2.f ----- the quintic spline function.
c
c**************************************************************
      implicit double precision (a-h,o-z)
      dimension cpt(2),cjt(2),anode(2)
      dimension p(6),pdx(6),pdy(6)
      dimension pdxx(6),pdxy(6),pdyy(6)
      dimension b(6),bd(2,6),bdd(3,6)
      dimension shpd(2),shpdd(3)
c
c.....find the polynominal p
c
      xpt = cpt(1)
      ypt = cpt(2)
c
      xjt = cjt(1)
      yjt = cjt(2)
c
c.....assign value to working array
c
      a1 =  anode(1)
      a2 =  anode(2)
c
      p(1) = 1.0d0
      p(2) = (xjt - xpt)/a1
      p(3) = (yjt - ypt)/a2
      p(4) = p(2) * p(2)
      p(5) = p(2) * p(3)
      p(6) = p(3) * p(3)
c
      dx   =  -1.0/a1
      dy   =  -1.0/a2
c      
      pdx(1) = 0.0
      pdx(2) = dx
      pdx(3) = 0.0
      pdx(4) = 2.0*dx*p(2)
      pdx(5) = dx*p(3)
      pdx(6) = 0.0
c
      pdy(1) = 0.0
      pdy(2) = 0.0
      pdy(3) = dy
      pdy(4) = 0.0
      pdy(5) = dy*p(2)
      pdy(6) = 2.0*dy*p(3)
c
      pdxx(1) = 0.0
      pdxx(2) = 0.0
      pdxx(3) = 0.0
      pdxx(4) = 2.0*dx*dx
      pdxx(5) = 0.0
      pdxx(6) = 0.0
c
      pdxy(1) = 0.0
      pdxy(2) = 0.0
      pdxy(3) = 0.0
      pdxy(4) = 0.0
      pdxy(5) = dx*dy
      pdxy(6) = 0.0
c
      pdyy(1) = 0.0
      pdyy(2) = 0.0
      pdyy(3) = 0.0
      pdyy(4) = 0.0
      pdyy(5) = 0.0
      pdyy(6) = 2.0*dy*dy
c
c
c.....compute correct function  C
c
      corf   = 0.00d0
      cdx    = 0.00d0
      cdy    = 0.00d0
c
      do i = 1,6
	 corf  = corf + p(i)*b(i)
	 cdx   = cdx + pdx(i)*b(i) + p(i)*bd(1,i)
	 cdy   = cdy + pdy(i)*b(i) + p(i)*bd(2,i)
      enddo
c
c
c....compute the second derivatives
c
      cdxx = 0.0d0
      cdxy = 0.0d0
      cdyy = 0.0d0
c
      do i = 1,6
         cdxx = cdxx + pdxx(i)*b(i) 
     &        + 2.0*pdx(i)*bd(1,i) + p(i)*bdd(1,i)
         cdxy = cdxy + pdxy(i)*b(i) + pdx(i)*bd(2,i) 
     &        + pdy(i)*bd(1,i) + p(i)*bdd(2,i)
         cdyy = cdyy + pdyy(i)*b(i)
     &        + 2.0*pdy(i)*bd(2,i) + p(i)*bdd(3,i)
      enddo
c
c.....check whether or not the point is within the compact support
c
      xx = dabs(xjt - xpt)/a1
      yy = dabs(yjt - ypt)/a2

      if (xx .gt. 2.0 .or. yy .gt. 2.0) then
	  shp      = 0.0d0
	  shpd(1)  = 0.0d0
	  shpd(2)  = 0.0d0
	  shpdd(1) = 0.0d0
	  shpdd(2) = 0.0d0
	  shpdd(3) = 0.0d0
      else
	  call window(aw,awdx,awdy,awdxx,awdxy,awdyy,
     &                xpt,ypt,xjt,yjt,a1,a2)
c
c.....calculate the shape functions
c
          dwj   = wjt
	  aw    = aw * dwj
	  awdx  = awdx * dwj
	  awdy  = awdy * dwj
	  awdxx = awdxx * dwj
	  awdxy = awdxy * dwj
	  awdyy = awdyy * dwj
c
          shp   = corf * aw
c
c.....calculate the 1st derivatives of the shape function
c
          shpd(1) = cdx*aw + corf*awdx
          shpd(2) = cdy*aw + corf*awdy
c
c.....calculate the 2nd derivatives of the shape function
c
          shpdd(1) = cdxx*aw + 2.0*cdx*awdx + corf*awdxx 
          shpdd(2) = cdxy*aw + cdx*awdy + cdy*awdx + corf*awdxy 
          shpdd(3) = cdyy*aw + 2.0*cdy*awdy + corf*awdyy 
c
      endif
c
  999 continue
c
      return
      end
c
