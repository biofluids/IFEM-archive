      subroutine RKPMshape3d(shp,shpd,b,bd,cpt,cjt,haj,wjt,iInter)
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
      real*8 b(*),bd(3,*),haj(3)
      real*8 wjt
      integer iInter
c
      if (iInter .eq. 1) then
         call RKPMshape3dl(shp,shpd,b,bd,cpt,cjt,haj,wjt)
      elseif(iInter .eq. 11) then
         call RKPMshape3dtl(shp,shpd,b,bd,cpt,cjt,haj,wjt)
      endif
c
c
      return
      end
c
