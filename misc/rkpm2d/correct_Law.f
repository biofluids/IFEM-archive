      subroutine correct_Law(b,bd,wx,wxd,wy,wyd,
     &                       cpt,cjp,dcjp,dwjp,nep)
c**************************************************************
c
c     This subroutine is to calculate the b vector and its
c     derivatives,for the moving least square reproducing kernel 
c     interpolation shape function.
c
c     This code is only offering b vector and its 1st derivatives
c     for the 2-D case (RGT1a).
c
c     The interpolation function is choosen as
c     
c     p = ( 1 , x , y )
c
c     Date: June, 1994
c
c     -arguments:
c
c  i      nep: the numbers of particles, therefore, the numbers
c              of elements too.
c
c  i      cjp(2,nep):  the particle's global coordinations;
c                      such as x(nep):= cjp(1,nep)
c                              y(nep):= cjp(2,nep)
c
c  i      dcjp(2,nep): the array that stores all the dilation
c                     vectors of the particles .
c                     The dilation parameter at X direction: dcjp(1,nep)
c                     The dilation parameter at Y direction: dcjp(2,nep)
c
c  i      dwjp(nep): the integration weight at each particle point
c
c
c  i      cpt(2): the point at where the shape function is
c                     evaluated.;
c                 xp:= cpt(1)
c                 yp:= cpt(2)
c
c    
c   l     ha1 : dilation parameter; scalar;
c
c   l     ha2 : dilation parameter; scalar;
c
c               we input the adjust constants ax, ay, and we construct
c               the dilation parameter in such a way: 
c               (for non-uniform mesh )
c
c                 ha1 = dcjp(1,j)
c                 ha2 = dcjp(2,j)
c
c    l   gm(3,3) : the basic moment matrix;
c
c    l   gminv(3,3): the inverse of the gm, i.e. gm^{-1}
c
c    l   gmdx(3,3): the derivatives of gm: gm_x;
c
c    l   gmdy(3,3): the derivatives of gm: gm_y;
c
c    l   adx(3,3):  M^(-1) M_x
c
c    l   ady(3,3):  M^(-1) M_y
c
c
c  o     b(3)   :   the b vector, also the first row of gminv
c
c  o     bd(2,3):   the 1st derivatives of the b vector
c
c                   bd(1,3) :=  b_x
c                   bd(2,3) :=  b_y
c
c     The subroutine is only designed to generate 2-D
c     moving least square shape function based on linear polynomial,
c
c    l   p(3) := (1,x,y);
c
c     The mathematical formulation is as follows
c
c     b   = (1/det)* [ a11, - a12, a13 ];
c
c     b_x = - M^{-1} M_x b;
c
c     b_y = - M^{-1} M_y b;
c
c
c    There two different subroutines could be called by this subroutine:
c
c    cubic2d.f ----- the cubic spline window function;
c  
c    quintw2.f ----- the fifth order spline window function.
c
c**************************************************************
      implicit double precision (a-h,o-z)
c
      dimension b(*),bd(2,*),cpt(2)
      dimension wx(*),wxd(2,*),wy(*),wyd(2,*)
      dimension cjp(2,nep),dcjp(2,nep),dwjp(nep)
      dimension gminv(3,3),gmdx(3,3),gmdy(3,3)
      dimension adx(3,3),ady(3,3)
c
c.....set the initial value for moment and its derivatives:
c
      zero = 0.00
c
      am00 = 0.00
      am10 = 0.00
      am01 = 0.00
      am20 = 0.00
      am11 = 0.00
      am02 = 0.00
c
c
      am00dx = 0.00
      am10dx = 0.00
      am01dx = 0.00
      am20dx = 0.00
      am11dx = 0.00
      am02dx = 0.00
c
      am00dy = 0.00
      am10dy = 0.00
      am01dy = 0.00
      am20dy = 0.00
      am11dy = 0.00
      am02dy = 0.00
c
c
c.....set the initial value for all array: 
c
      do i = 1, 3
	 do  j = 1, 3
	    gminv(i,j)  = 0.00
	    gmdx(i,j)   = 0.00
	    gmdy(i,j)   = 0.00
         end do
	 b(i)     = 0.00
	 bd(1,i)  = 0.00
	 bd(2,i)  = 0.00
c
	 wx(i)    = 0.00
	 wxd(1,i) = 0.00
	 wxd(2,i) = 0.00
c
	 wy(i)    = 0.00
	 wyd(1,i) = 0.00
	 wyd(2,i) = 0.00
c
      end do  
c
c.....input the value cpt 
c
      xp = cpt(1)
      yp = cpt(2)
c
c.....main loop: calculate moment by Trapezodial rule
c
      do 30 j = 1, nep
c
c........define intermediate variable
c
c........ha1 and ha2
c
	 ha1 =  dcjp(1,j)
	 ha2 =  dcjp(2,j)
c
	 dsj = dwjp(j)

	 xj  = cjp(1,j)
	 yj  = cjp(2,j)
         dx  = -1.0/ha1
	 dy  = -1.0/ha2
c
	 r10 = (xj - xp)/ha1
	 r01 = (yj - yp)/ha2
	 r20 = r10*r10
         r11 = r10*r01
         r02 = r01*r01
c
         xx = dabs(r10)
	 yy = dabs(r01)
c
	 if((xx.ge.2.0) .or. (yy.ge.2.0)) go to 30
c
         call window(aw,awdx,awdy,awdxx,awdxy,
     &               awdyy,xp,yp,xj,yj,ha1,ha2)
c
         aw   = aw*dsj
c
         am00 = am00  +  aw
         am10 = am10  +  r10*aw
         am01 = am01  +  r01*aw
         am20 = am20  +  r20*aw
         am11 = am11  +  r11*aw
         am02 = am02  +  r02*aw
c
c
	 awdx = awdx*dsj
	 awdy = awdy*dsj
c
         am00dx  = am00dx  + awdx
	 am10dx  = am10dx  + dx*aw + r10*awdx
	 am01dx  = am01dx  + r01*awdx
	 am20dx  = am20dx  + 2.0*dx*r10*aw + r20*awdx
	 am11dx  = am11dx  + dx*r01*aw + r11*awdx
	 am02dx  = am02dx  + r02*awdx
c
         am00dy  = am00dy  + awdy
	 am10dy  = am10dy  + r10*awdy
	 am01dy  = am01dy  + dy*aw + r01*awdy
	 am20dy  = am20dy  + r20*awdy
	 am11dy  = am11dy  + r10*dy*aw + r11*awdy
	 am02dy  = am02dy  + 2.0*dy*r01*aw + r02*awdy
c
  30  continue
c
c
c.....end of the main loop
c
c.....assemble the M matrix (gm(i,j))
c
c      gm(1,1) = am00
c      gm(1,2) = am10
c      gm(2,1) = am10
c      gm(1,3) = am01
c      gm(3,1) = am01
c      gm(2,2) = am20
c      gm(2,3) = am11
c      gm(3,2) = am11
c      gm(3,3) = am02
c
c
c.....assemble the cofactor matrices ( a(i,j))
c
      a11 =  am20*am02  - am11*am11
      a12 =  am10*am02  - am11*am01
      a21 =  a12
      a13 =  am10*am11  - am20*am01
      a31 =  a13
      a22 =  am00*am02  - am01*am01
      a23 =  am00*am11  - am10*am01
      a32 =  a23
      a33 =  am00*am20  - am10*am10
c
c.....calculate the determinat det
c
      det = am00*a11 - am10*a12 + am01*a13 
c
c.....test convergence criteria
c
      zero = 0.0d0
      if(det .le. zero) then
	print *, 'det =', det
	print *, 'STOP! the determinat det < 0 '
	print *, am00, am10, am01 
	print *, am10, am20, am11 
	print *, am01, am11, am02 
	stop
      else
      end if
c
c.....assemble the gminv(i,j)
c
      cdet = 1.0/det
      gminv(1,1) =  a11*cdet
      gminv(1,2) = -a12*cdet
      gminv(2,1) = -a12*cdet
      gminv(1,3) =  a13*cdet
      gminv(3,1) =  a13*cdet
      gminv(2,2) =  a22*cdet
      gminv(2,3) = -a23*cdet
      gminv(3,2) = -a23*cdet
      gminv(3,3) =  a33*cdet
c
c
c.....calculate the derivative of gm
c
      gmdx(1,1) = am00dx
      gmdx(1,2) = am10dx
      gmdx(2,1) = am10dx
      gmdx(1,3) = am01dx
      gmdx(3,1) = am01dx
      gmdx(2,2) = am20dx
      gmdx(2,3) = am11dx
      gmdx(3,2) = am11dx
      gmdx(3,3) = am02dx
c
      gmdy(1,1) = am00dy
      gmdy(1,2) = am10dy
      gmdy(2,1) = am10dy
      gmdy(1,3) = am01dy
      gmdy(3,1) = am01dy
      gmdy(2,2) = am20dy
      gmdy(2,3) = am11dy
      gmdy(3,2) = am11dy
      gmdy(3,3) = am02dy
c
c.....assign the value for b vector
c
      do i = 1,3
         b(i)  = gminv(1,i)
         wx(i) = gminv(2,i)
         wy(i) = gminv(3,i)
      enddo
c
c.....find the value for bd(2,3)
c
c.....compute M^(-1)dM/dx and M^(-1)dM/dy
c
      do i = 1,3
	 do j = 1,3
	    adx(i,j) = 0.0d0    
	    ady(i,j) = 0.0d0    
	    do k = 1,3
	       adx(i,j) = adx(i,j) + gminv(i,k)*gmdx(k,j)
	       ady(i,j) = ady(i,j) + gminv(i,k)*gmdy(k,j)
            enddo
         enddo
      enddo
c
c
      do i = 1,3
         do j = 1,3
	    bd(1,i) = bd(1,i) - adx(i,j)*b(j)
	    bd(2,i) = bd(2,i) - ady(i,j)*b(j)
c
	    wxd(1,i) = wxd(1,i) - adx(i,j)*wx(j)
	    wxd(2,i) = wxd(2,i) - ady(i,j)*wx(j)
c
	    wyd(1,i) = wyd(1,i) - adx(i,j)*wy(j)
	    wyd(2,i) = wyd(2,i) - ady(i,j)*wy(j)
         end do
      end do
c
c
c
  999 continue
c
      return
      end
c
