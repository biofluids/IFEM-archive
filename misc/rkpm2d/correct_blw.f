      subroutine correct_blw(b,bd,bdd,
     &           wx,wxd,wxdd,wy,wyd,wydd,
     &           wxy,wxyd,wxydd,
     &           cpt,cjp,anode,dwjp,nep)
c
c**************************************************************
c
c     This subroutine is to calculate the b vector and its
c     derivatives,for the moving least square reproducing kernel
c     shape functions.
c
c     This code is only offering b vector and its derivatives
c     up to second order for the 2-D case.
c     The interpolation function is choosen as
c     
c     p= ( 1, x , y, x*y ) .
c
c
c     Author: Shaofan Li
c
c     Date: Octorber, 1996
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
c  i      anode(2,nep): the array that contains particle's
c                      dilation vector, i.e.
c                          a1(ip) =  anode(1,ip)
c                          a2(ip) =  anode(2,ip)
c
c  i      dwjp(nep): the wight at each particle point
c
c         The only difference between dcjp and dwjp is at the
c         boundary.
c
c  i      cpt(2): the point at where the shape function is
c                     evaluated.;
c                 x(pt):= cpt(1)
c                 y(pt):= cpt(2)
c
c    
c   Note:
c   There is differences between dilation parameter and dilation
c   coeffecients, i.e.
c
c   l     a1, a2:  dilation parameter; scalar;
c                  for uniformal mesh, we choose a1:=ax*dx
c                                                a2:=ay*dy
c         in which, the dilation coefficients are, roughly, 1 <ax,ay < 3
c
c                 For non-uniform mesh, 
c                 we input the adjust constants ax, ay, and 
c                 the dilation parameters are constructed in
c                 an abstract way, i.e. 
c
c                 a1 =  dcjp(1,j)
c                 a2 =  dcjp(2,j)  .
c
c    l   gm(4,4) : the basic moment matrix;
c
c    l   gminv(4,4): the inverse of the gm, i.e. gm^{-1}
c
c    l   gmdx(4,4): the derivatives of gm: d/dx (gm);
c
c    l   gmdy(4,4): the derivatives of gm: d/dy (gm);
c
c    l   gmdxx(4,4): d^2/dx^2 (gm)
c
c    l   gmdxy(4,4): d^2/dxy (gm)
c
c    l   gmdyy(4,4): d^2/dy^2 (gm)
c
c    l   adx(4,4):  M^(-1) d/dx M
c
c    l   ady(4,4):  M^(-1) d/dy M
c
c    l   adxx(4,4):  local array to compute b,xx
c
c        adxx := M^{-1} M_xx
c
c    l   adxy(4,4):  local array to compute b,xy
c
c        adxy := M^{-1} M_xy
c
c    l   adyy(4,4):  local array to compute b,yy
c
c        adyy = M^{-1} M_yy
c
c  o     b(4)   :   the b vector, also the first row of gminv
c
c  o     bd(2,4):   the 1st derivatives of the b vector
c
c                   bd(1,:) := d/dx b
c                   bd(2,:) := d/dy b
c
c  o     bdd(3,4):  the 2nd derivatives of the b vector
c
c                   bdd(1,:) := d^2/dxx b
c                   bdd(2,:) := d^2/dxy b
c                   bdd(3,:) := d^2/dyy b
c
c     The subroutine is only designed for solving 2-D locally
c     linear interpolation moving field problem. i.e.
c
c    l   p0   := (1,0,0,0); 
c
c    l   p(4) := (1,x,y,x*y);
c
c    There two different subroutines are called by this subroutine:
c
c    cubic2d.f ----- the cubic spline window function;
c  
c    quintw2.f ----- the quintic spline function.
c
c**************************************************************
      implicit double precision (a-h,o-z)
c
      dimension cjp(2,nep),dwjp(nep),anode(2,nep)
c
      dimension gminv(4,4),gmdx(4,4),gmdy(4,4)
      dimension gmdxx(4,4),gmdxy(4,4),gmdyy(4,4)
      dimension adx(4,4),ady(4,4)
      dimension adxx(4,4),adxy(4,4),adyy(4,4)
c
      dimension b(*),bd(2,*),bdd(3,*),cpt(2)
      dimension wx(*),wxd(2,*),wxdd(3,*)
      dimension wy(*),wyd(2,*),wydd(3,*)
      dimension wxy(*),wxyd(2,*),wxydd(3,*)
c
c.....set the initial value for moment and its derivatives:
c
      zero    = 0.00d0
c
      am00    = 0.00d0
      am10    = 0.00d0
      am01    = 0.00d0
      am11    = 0.00d0
c
      am20    = 0.00d0
      am02    = 0.00d0
      am21    = 0.00d0
      am12    = 0.00d0
      am22    = 0.00d0
c
      am00dx  = 0.00d0
      am10dx  = 0.00d0
      am01dx  = 0.00d0
      am11dx  = 0.00d0
c
      am20dx  = 0.00d0
      am02dx  = 0.00d0
      am21dx  = 0.00d0
      am12dx  = 0.00d0
      am22dx  = 0.00d0
c
      am00dy  = 0.00d0
      am10dy  = 0.00d0
      am01dy  = 0.00d0
      am11dy  = 0.00d0
c
      am20dy  = 0.00d0
      am02dy  = 0.00d0
      am21dy  = 0.00d0
      am12dy  = 0.00d0
      am22dy  = 0.00d0
c
      am00dxx = 0.00d0
      am10dxx = 0.00d0
      am01dxx = 0.00d0
      am11dxx = 0.00d0
c
      am20dxx = 0.00d0
      am02dxx = 0.00d0
      am21dxx = 0.00d0
      am12dxx = 0.00d0
      am22dxx = 0.00d0
c
c
      am00dxy = 0.00d0
      am10dxy = 0.00d0
      am01dxy = 0.00d0
      am11dxy = 0.00d0
c
      am20dxy = 0.00d0
      am02dxy = 0.00d0
      am21dxy = 0.00d0
      am12dxy = 0.00d0
      am22dxy = 0.00d0
c
c
      am00dyy = 0.00d0
      am10dyy = 0.00d0
      am01dyy = 0.00d0
      am11dyy = 0.00d0
c
      am20dyy = 0.00d0
      am02dyy = 0.00d0
      am21dyy = 0.00d0
      am12dyy = 0.00d0
      am22dyy = 0.00d0
c
c.....set the initial value for all array: 
c
      do i = 1, 4
	 do  j = 1, 4
	    gminv(i,j)  = 0.00d0
	    gmdx(i,j)   = 0.00d0
	    gmdy(i,j)   = 0.00d0
	    gmdxx(i,j)  = 0.00d0
	    gmdxy(i,j)  = 0.00d0
	    gmdyy(i,j)  = 0.00d0
c
	    adx(i,j)    = 0.0d0
	    ady(i,j)    = 0.0d0
	    adxx(i,j)   = 0.0d0
	    adxy(i,j)   = 0.0d0
	    adyy(i,j)   = 0.0d0
         end do
	 b(i)     = 0.00d0
	 bd(1,i)  = 0.00d0
	 bd(2,i)  = 0.00d0
	 bdd(1,i) = 0.00d0
	 bdd(2,i) = 0.00d0
	 bdd(3,i) = 0.00d0
c
	 wx(i)     = 0.00d0
	 wxd(1,i)  = 0.00d0
	 wxd(2,i)  = 0.00d0
	 wxdd(1,i) = 0.00d0
	 wxdd(2,i) = 0.00d0
	 wxdd(3,i) = 0.00d0
c
	 wy(i)     = 0.00d0
	 wyd(1,i)  = 0.00d0
	 wyd(2,i)  = 0.00d0
	 wydd(1,i) = 0.00d0
	 wydd(2,i) = 0.00d0
	 wydd(3,i) = 0.00d0
c
	 wxy(i)     = 0.00d0
	 wxyd(1,i)  = 0.00d0
	 wxyd(2,i)  = 0.00d0
	 wxydd(1,i) = 0.00d0
	 wxydd(2,i) = 0.00d0
	 wxydd(3,i) = 0.00d0
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
         a1  = anode(1,j)
         a2  = anode(2,j)
	 dsj = dwjp(j)
c
	 xj  = cjp(1,j)
	 yj  = cjp(2,j)
	 r10 = (xj - xp)/a1
	 r01 = (yj - yp)/a2
	 r11 = r10*r01
c
	 r20 = r10*r10
	 r02 = r01*r01
	 r21 = r20*r01
	 r12 = r10*r02
	 r22 = r20*r02
c
         dx  = - 1.0d0/a1
	 dy  = - 1.0d0/a2
c
         dxx = dx * dx
	 dxy = dx * dy
	 dyy = dy * dy
c
         xx = dabs(r10)
	 yy = dabs(r01)

	 if(xx.gt. 2.0d0 .or. yy.gt. 2.0d0 ) go to 30
c
c
         call window(aw,awdx,awdy,awdxx,awdxy,
     &               awdyy,xp,yp,xj,yj,a1,a2)
c
	 aw    = aw*dsj
c
         am00   = am00  + aw
         am10   = am10  + r10*aw
         am01   = am01  + r01*aw
         am11   = am11  + r11*aw
c
	 am20   = am20  + r20*aw
         am02   = am02  + r02*aw
         am21   = am21  + r21*aw
         am12   = am12  + r12*aw
         am22   = am22  + r22*aw
c
	 awdx = awdx*dsj
	 awdy = awdy*dsj
c
         am00dx  = am00dx + awdx
	 am10dx  = am10dx + dx*aw + r10*awdx
	 am01dx  = am01dx + r01*awdx
	 am11dx  = am11dx + dx*r01*aw + r11*awdx
c
	 am20dx  = am20dx + 2.0d0*dx*r10*aw + r20*awdx
	 am02dx  = am02dx + r02*awdx
	 am21dx  = am21dx + 2.0d0*dx*r11*aw + r21*awdx
         am12dx  = am12dx + dx*r02*aw + r12*awdx
	 am22dx  = am22dx + 2.0d0*dx*r12*aw + r22*awdx
c
c
         am00dy  = am00dy + awdy
	 am10dy  = am10dy + r10*awdy
	 am01dy  = am01dy + dy*aw + r01*awdy
	 am11dy  = am11dy + r10*dy*aw + r11*awdy
c
	 am20dy  = am20dy + r20*awdy
	 am02dy  = am02dy + 2.0d0*r01*dy*aw + r02*awdy
	 am21dy  = am21dy + dy*r20*aw + r21*awdy
         am12dy  = am12dy + 2.0d0*dy*r11*aw + r12*awdy
	 am22dy  = am22dy + 2.0d0*dy*r21*aw + r22*awdy
c
c........the second derivatives
c
	 awdxx = awdxx*dsj
	 awdxy = awdxy*dsj
	 awdyy = awdyy*dsj
c
	 am00dxx = am00dxx + awdxx
	 am10dxx = am10dxx + 2.0d0*dx*awdx + r10*awdxx
         am01dxx = am01dxx + r01*awdxx
	 am11dxx = am11dxx + 2.0d0*dx*r01*awdx + r11*awdxx
c
	 am20dxx = am20dxx + 2.0d0*dxx*aw + 4.0d0*dx*r10*awdx 
     &           + r20*awdxx
	 am02dxx = am02dxx + r02*awdxx
	 am21dxx = am21dxx + 2.0d0*dxx*r01*aw 
     &           + 4.0d0*dx*r11*awdx + r21*awdxx
         am12dxx = am12dxx + 2.0d0*dx*r02*awdx + r12*awdxx
	 am22dxx = am22dxx + 2.0d0*dxx*r02*aw 
     &           + 4.0d0*dx*r12*awdx + r22*awdxx
c
c.................
c
	 am00dxy = am00dxy + awdxy
	 am10dxy = am10dxy + dx*awdy + r10*awdxy
	 am01dxy = am01dxy + dy*awdx + r01*awdxy
	 am11dxy = am11dxy + dxy*aw + r10*dy*awdx + dx*r01*awdy 
     &           + r11*awdxy
c
	 am20dxy = am20dxy + 2.0d0*dx*r10*awdy + r20*awdxy
         am02dxy = am02dxy + 2.0d0*dy*r01*awdx + r02*awdxy
	 am21dxy = am21dxy + 2.0d0*dxy*r10*aw + r20*dy*awdx 
     &           + 2.0d0*dx*r11*awdy + r21*awdxy
	 am12dxy = am12dxy + 2.0d0*dxy*r01*aw + r02*dx*awdy 
     &           + 2.0d0*dy*r11*awdx + r12*awdxy
	 am22dxy = am22dxy + 4.0d0*dxy*r11*aw + 2.0d0*r21*dy*awdx
     &           + 2.0d0*r12*dx*awdy + r22*awdxy
c
c.................
c
	 am00dyy = am00dyy + awdyy
	 am10dyy = am10dyy + r10*awdyy
         am01dyy = am01dyy + 2.0d0*dy*awdy + r01*awdyy
	 am11dyy = am11dyy + 2.0d0*r10*dy*awdy + r11*awdyy
c
	 am20dyy = am20dyy + r20*awdyy
	 am02dyy = am02dyy + 2.0d0*dyy*aw + 4.0d0*dy*r01*awdy 
     &           + r02*awdyy
	 am21dyy = am21dyy + 2.0d0*r20*dy*awdy + r21*awdyy
	 am12dyy = am12dyy + 2.0d0*r10*dyy*aw 
     &           + 4.0d0*r11*dy*awdy + r12*awdyy
	 am22dyy = am22dyy + 2.0d0*r20*dyy*aw 
     &           + 4.0d0*r21*dy*awdy +r22*awdyy
c
  30  continue
c
c.....end of the main loop
c
c.....assemble the M matrix (gm(i,j))
c
      f11 = am20*(am02*am22 - am12*am12)
     &    - am11*(am11*am22 - am12*am21)
     &    + am21*(am11*am12 - am02*am21)
c
      f12 = am10*(am02*am22 - am12*am12)
     &    - am11*(am01*am22 - am12*am11)
     &    + am21*(am01*am12 - am02*am11)
c
      f13 = am10*(am11*am22 - am12*am21)
     &    - am20*(am01*am22 - am12*am11)
     &    + am21*(am01*am21 - am11*am11)
c
      f14 = am10*(am11*am12 - am02*am21)
     &    - am20*(am01*am12 - am02*am11)
     &    + am11*(am01*am21 - am11*am11)
c
      f21 = f12
c
      f22 = am00*(am02*am22 - am12*am12)
     &    - am01*(am01*am22 - am12*am11)
     &    + am11*(am01*am12 - am02*am11)
c
      f23 = am00*(am11*am22 - am12*am21)
     &    - am10*(am01*am22 - am12*am11)
     &    + am11*(am01*am21 - am11*am11)
c
      f24 = am00*(am11*am12 - am02*am21)
     &    - am10*(am01*am12 - am02*am11)
     &    + am01*(am01*am21 - am11*am11)
c
      f31 = f13
c
      f32 = f23
c
      f33 = am00*(am20*am22 - am21*am21)
     &    - am10*(am10*am22 - am21*am11)
     &    + am11*(am10*am21 - am20*am11)
c
      f34 = am00*(am20*am12 - am11*am21)
     &    - am10*(am10*am12 - am11*am11)
     &    + am01*(am10*am21 - am20*am11)
c
      f41  = f14
c
      f42  = f24
c
      f43  = f34
c
      f44  = am00*(am20*am02 - am11*am11)
     &     - am10*(am10*am02 - am11*am01)
     &     + am01*(am10*am11 - am20*am01)
c
c.....find the inverse
c
      det = am00*f11- am10*f12 + am01*f13 - am11*f14
c
c.....test convergence criteria
c
      if(det .le. zero) then
	print *, 'det =', det
	print *, 'singular pt.', xp, yp
	print *, 'STOP! the determinat det < 0 '
	stop
      else
      end if
c
      cdet = 1.0d0/det
c
      gminv(1,1) =   f11 * cdet
      gminv(1,2) = - f12 * cdet
      gminv(1,3) =   f13 * cdet
      gminv(1,4) = - f14 * cdet
c
      gminv(2,1) = - f21 * cdet
      gminv(2,2) =   f22 * cdet
      gminv(2,3) = - f23 * cdet
      gminv(2,4) =   f24 * cdet
c
      gminv(3,1) =   f31 * cdet
      gminv(3,2) = - f32 * cdet
      gminv(3,3) =   f33 * cdet
      gminv(3,4) = - f34 * cdet
c
      gminv(4,1) = - f41 * cdet
      gminv(4,2) =   f42 * cdet
      gminv(4,3) = - f43 * cdet
      gminv(4,4) =   f44 * cdet
c
c.....calculate the derivative of gm
c
      gmdx(1,1) = am00dx
      gmdx(1,2) = am10dx
      gmdx(1,3) = am01dx
      gmdx(1,4) = am11dx
c
      gmdx(2,1) = am10dx
      gmdx(2,2) = am20dx
      gmdx(2,3) = am11dx
      gmdx(2,4) = am21dx
c
      gmdx(3,1) = am01dx
      gmdx(3,2) = am11dx
      gmdx(3,3) = am02dx
      gmdx(3,4) = am12dx
c
      gmdx(4,1) = am11dx
      gmdx(4,2) = am21dx
      gmdx(4,3) = am12dx
      gmdx(4,4) = am22dx
c
      gmdy(1,1) = am00dy
      gmdy(1,2) = am10dy
      gmdy(1,3) = am01dy
      gmdy(1,4) = am11dy
c
      gmdy(2,1) = am10dy
      gmdy(2,2) = am20dy
      gmdy(2,3) = am11dy
      gmdy(2,4) = am21dy
c
      gmdy(3,1) = am01dy
      gmdy(3,2) = am11dy
      gmdy(3,3) = am02dy
      gmdy(3,4) = am12dy
c
      gmdy(4,1) = am11dy
      gmdy(4,2) = am21dy
      gmdy(4,3) = am12dy
      gmdy(4,4) = am22dy
c
c.....calculate the second derivative of gm
c 
      gmdxx(1,1) = am00dxx
      gmdxx(1,2) = am10dxx
      gmdxx(1,3) = am01dxx
      gmdxx(1,4) = am11dxx
c
      gmdxx(2,1) = am10dxx
      gmdxx(2,2) = am20dxx
      gmdxx(2,3) = am11dxx
      gmdxx(2,4) = am21dxx
c
      gmdxx(3,1) = am01dxx
      gmdxx(3,2) = am11dxx
      gmdxx(3,3) = am02dxx
      gmdxx(3,4) = am12dxx
c
      gmdxx(4,1) = am11dxx
      gmdxx(4,2) = am21dxx
      gmdxx(4,3) = am12dxx
      gmdxx(4,4) = am22dxx
c
c 
      gmdxy(1,1) = am00dxy
      gmdxy(1,2) = am10dxy
      gmdxy(1,3) = am01dxy
      gmdxy(1,4) = am11dxy
c
      gmdxy(2,1) = am10dxy
      gmdxy(2,2) = am20dxy
      gmdxy(2,3) = am11dxy
      gmdxy(2,4) = am21dxy
c
      gmdxy(3,1) = am01dxy
      gmdxy(3,2) = am11dxy
      gmdxy(3,3) = am02dxy
      gmdxy(3,4) = am12dxy
c
      gmdxy(4,1) = am11dxy
      gmdxy(4,2) = am21dxy
      gmdxy(4,3) = am12dxy
      gmdxy(4,4) = am22dxy
c
c
      gmdyy(1,1) = am00dyy
      gmdyy(1,2) = am10dyy
      gmdyy(1,3) = am01dyy
      gmdyy(1,4) = am11dyy
c
      gmdyy(2,1) = am10dyy
      gmdyy(2,2) = am20dyy
      gmdyy(2,3) = am11dyy
      gmdyy(2,4) = am21dyy
c
      gmdyy(3,1) = am01dyy
      gmdyy(3,2) = am11dyy
      gmdyy(3,3) = am02dyy
      gmdyy(3,4) = am12dyy
c
      gmdyy(4,1) = am11dyy
      gmdyy(4,2) = am21dyy
      gmdyy(4,3) = am12dyy
      gmdyy(4,4) = am22dyy
c
c.....assign the value for b vector
c
      do i = 1, 4
         b(i)   = gminv(1,i)
         wx(i)  = gminv(2,i)
         wy(i)  = gminv(3,i)
         wxy(i) = gminv(4,i)
      enddo
c
c.....find the value for bd(2,3)
c
c.....compute M^(-1)dM/dx and M^(-1)dM/dy
c
      do i = 1,4
	 do j = 1,4
	    do k = 1,4
	       adx(i,j)  = adx(i,j)  + gminv(i,k)*gmdx(k,j)
	       ady(i,j)  = ady(i,j)  + gminv(i,k)*gmdy(k,j)
               adxx(i,j) = adxx(i,j) + gminv(i,k)*gmdxx(k,j)
               adxy(i,j) = adxy(i,j) + gminv(i,k)*gmdxy(k,j)
               adyy(i,j) = adyy(i,j) + gminv(i,k)*gmdyy(k,j)
            enddo
         enddo
      enddo
c
c
      do i = 1,4
         do j = 1,4
	    bd(1,i) = bd(1,i) - adx(i,j)*b(j)
	    bd(2,i) = bd(2,i) - ady(i,j)*b(j)
c
	    wxd(1,i) = wxd(1,i) - adx(i,j)*wx(j)
	    wxd(2,i) = wxd(2,i) - ady(i,j)*wx(j)
c	    
	    wyd(1,i) = wyd(1,i) - adx(i,j)*wy(j)
	    wyd(2,i) = wyd(2,i) - ady(i,j)*wy(j)
c
	    wxyd(1,i) = wxyd(1,i) - adx(i,j)*wxy(j)
	    wxyd(2,i) = wxyd(2,i) - ady(i,j)*wxy(j)
c
         enddo
      enddo
c
c
c....compute the 2nd derivatives bdd(3,3)
c
c
      do i = 1,4
	 do j = 1,4
	    bdd(1,i) = bdd(1,i) - 2.0d0*adx(i,j)*bd(1,j)
     &               - adxx(i,j)*b(j)
	    bdd(2,i) = bdd(2,i) - adx(i,j)*bd(2,j)
     &               - ady(i,j)*bd(1,j) - adxy(i,j)*b(j)
            bdd(3,i) = bdd(3,i) - 2.0d0*ady(i,j)*bd(2,j)
     &               - adyy(i,j)*b(j)
c
	    wxdd(1,i) = wxdd(1,i) - 2.0d0*adx(i,j)*wxd(1,j)
     &               - adxx(i,j)*wx(j)
	    wxdd(2,i) = wxdd(2,i) - adx(i,j)*wxd(2,j)
     &               - ady(i,j)*wxd(1,j) - adxy(i,j)*wx(j)
            wxdd(3,i) = wxdd(3,i) - 2.0d0*ady(i,j)*wxd(2,j)
     &               - adyy(i,j)*wx(j)
c
	    wydd(1,i) = wydd(1,i) - 2.0d0*adx(i,j)*wyd(1,j)
     &               - adxx(i,j)*wy(j)
	    wydd(2,i) = wydd(2,i) - adx(i,j)*wyd(2,j)
     &               - ady(i,j)*wyd(1,j) - adxy(i,j)*wy(j)
            wydd(3,i) = wydd(3,i) - 2.0d0*ady(i,j)*wyd(2,j)
     &               - adyy(i,j)*wy(j)
c
c
	    wxydd(1,i) = wxydd(1,i) - 2.0d0*adx(i,j)*wxyd(1,j)
     &               - adxx(i,j)*wxy(j)
	    wxydd(2,i) = wxydd(2,i) - adx(i,j)*wxyd(2,j)
     &               - ady(i,j)*wxyd(1,j) - adxy(i,j)*wxy(j)
            wxydd(3,i) = wxydd(3,i) - 2.0d0*ady(i,j)*wxyd(2,j)
     &               - adyy(i,j)*wxy(j)
c
         enddo
      enddo
c
c
  999 continue
c
      return
      end
c
