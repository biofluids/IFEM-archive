      subroutine correct_qdw(b,bd,bdd,wx,wxd,wxdd,wy,wyd,wydd,
     &           wxx,wxxd,wxxdd,wxy,wxyd,wxydd,
     &           wyy,wyyd,wyydd,cpt,cjp,anode,dwjp,nep)
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
c     p= ( 1, x , y, x*x, x*y, y*y ) .
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
c   Note:
c   There is differences between dilation parameter and dilation
c   coeffecients, i.e.
c
c   l     a1, a2:  dilation parameter; scalar;
c                  for uniformal mesh, we choose a1:=ax*dx
c                                                a2:=ay*dy
c         in which, the dilation coefficients are, roughly, 1 < ax,ay < 3
c
c                 For non-uniform mesh, 
c                 we input the adjust constants ax, ay, and 
c                 the dilation parameters are constructed in
c                 an abstract way, i.e. 
c
c                 a1 =  dcjp(1,j)
c                 a2 =  dcjp(2,j)  .
c
c    l   gm(6,6) : the basic moment matrix;
c
c    l   gminv(6,6): the inverse of the gm, i.e. gm^{-1}
c
c    l   gmdx(6,6): the derivatives of gm: d/dx (gm);
c
c    l   gmdy(6,6): the derivatives of gm: d/dy (gm);
c
c    l   gmdxx(6,6): d^2/dx^2 (gm)
c
c    l   gmdxy(6,6): d^2/dxy (gm)
c
c    l   gmdyy(6,6): d^2/dy^2 (gm)
c
c    l   adx(6,6):  M^(-1) d/dx M
c
c    l   ady(6,6):  M^(-1) d/dy M
c
c    l   adxx(6,6):  local array to compute b,xx
c
c        adxx := M^{-1} M_xx
c
c    l   adxy(6,6):  local array to compute b,xy
c
c        adxy := M^{-1} M_xy
c
c    l   adyy(6,6):  local array to compute b,yy
c
c        adyy = M^{-1} M_yy
c
c  o     b(6)   :   the b vector, also the first row of gminv
c
c  o     bd(2,6):   the 1st derivatives of the b vector
c
c                   bd(1,:) := d/dx b
c                   bd(2,:) := d/dy b
c
c  o     bdd(3,6):  the 2nd derivatives of the b vector
c
c                   bdd(1,:) := d^2/dxx b
c                   bdd(2,:) := d^2/dxy b
c                   bdd(3,:) := d^2/dyy b
c
c     The subroutine is only designed for solving 2-D locally
c     linear interpolation moving field problem. i.e.
c
c    l   p0 := (1,0,0,0,0,0); 
c
c    l   p(3) := (1,x,y,x*x,x*y,y*y);
c
c    There two different subroutines are called by this subroutine:
c
c    window:
c
c    cubic2d.f ----- the cubic spline window function;
c  
c    quintw2.f ----- the quintic spline function.
c
c**************************************************************
      implicit double precision (a-h,o-z)
c
      dimension cjp(2,nep),cpt(2),anode(2,nep),dwjp(nep)
      dimension gm(6,6),gminv(6,6),gmdx(6,6),gmdy(6,6)
      dimension gmdxx(6,6),gmdxy(6,6),gmdyy(6,6)
      dimension adx(6,6),ady(6,6)
      dimension adxx(6,6),adxy(6,6),adyy(6,6)
      dimension b(6),bd(2,6),bdd(3,6)
      dimension wx(6),wxd(2,6),wxdd(3,6)
      dimension wy(6),wyd(2,6),wydd(3,6)
      dimension wxx(6),wxxd(2,6),wxxdd(3,6)
      dimension wxy(6),wxyd(2,6),wxydd(3,6)
      dimension wyy(6),wyyd(2,6),wyydd(3,6)
c
c.....set the initial value for moment and its derivatives:
c
      zero   = 0.00
c
      am00    = 0.00
      am10    = 0.00
      am01    = 0.00
      am20    = 0.00
      am11    = 0.00
      am02    = 0.00
c
      am21    = 0.00
      am12    = 0.00
      am30    = 0.00
      am03    = 0.00
      am22    = 0.00
      am31    = 0.00
      am13    = 0.00
      am40    = 0.00
      am04    = 0.00
c
      am00dx  = 0.00
      am10dx  = 0.00
      am01dx  = 0.00
      am20dx  = 0.00
      am11dx  = 0.00
      am02dx  = 0.00
c
      am21dx  = 0.00
      am12dx  = 0.00
      am30dx  = 0.00
      am03dx  = 0.00
      am22dx  = 0.00
      am31dx  = 0.00
      am13dx  = 0.00
      am40dx  = 0.00
      am04dx  = 0.00
c
      am00dy  = 0.00
      am10dy  = 0.00
      am01dy  = 0.00
      am20dy  = 0.00
      am11dy  = 0.00
      am02dy  = 0.00
c
      am21dy  = 0.00
      am12dy  = 0.00
      am30dy  = 0.00
      am03dy  = 0.00
      am22dy  = 0.00
      am31dy  = 0.00
      am13dy  = 0.00
      am40dy  = 0.00
      am04dy  = 0.00
c
      am00dxx = 0.00
      am10dxx = 0.00
      am01dxx = 0.00
      am20dxx = 0.00
      am11dxx = 0.00
      am02dxx = 0.00
c
      am21dxx = 0.00
      am12dxx = 0.00
      am30dxx = 0.00
      am03dxx = 0.00
      am22dxx = 0.00
      am31dxx = 0.00
      am13dxx = 0.00
      am40dxx = 0.00
      am04dxx = 0.00
c
c
      am00dxy = 0.00
      am10dxy = 0.00
      am01dxy = 0.00
      am20dxy = 0.00
      am11dxy = 0.00
      am02dxy = 0.00
c
      am21dxy = 0.00
      am12dxy = 0.00
      am30dxy = 0.00
      am03dxy = 0.00
      am22dxy = 0.00
      am31dxy = 0.00
      am13dxy = 0.00
      am40dxy = 0.00
      am04dxy = 0.00
c
c
      am00dyy = 0.00
      am10dyy = 0.00
      am01dyy = 0.00
      am20dyy = 0.00
      am11dyy = 0.00
      am02dyy = 0.00
c
      am21dyy = 0.00
      am12dyy = 0.00
      am30dyy = 0.00
      am03dyy = 0.00
      am22dyy = 0.00
      am31dyy = 0.00
      am13dyy = 0.00
      am40dyy = 0.00
      am04dyy = 0.00
c
c.....set the initial value for all array: 
c
      do i = 1, 6
	 do  j = 1, 6
	    gm(i,j)     = 0.00
	    gminv(i,j)  = 0.00
	    gmdx(i,j)   = 0.00
	    gmdy(i,j)   = 0.00
	    gmdxx(i,j)  = 0.00
	    gmdxy(i,j)  = 0.00
	    gmdyy(i,j)  = 0.00
         end do
	 b(i)     = 0.00
	 bd(1,i)  = 0.00
	 bd(2,i)  = 0.00
	 bdd(1,i) = 0.00
	 bdd(2,i) = 0.00
	 bdd(3,i) = 0.00
      end do  
c
c.....input the value cpt 
c
      xpt = cpt(1)
      ypt = cpt(2)
c
c.....main loop: calculate moment by Trapezodial rule
c
      do 30 j = 1, nep
c
c........define intermediate variable
c
c........a1 and a2
c
	 a1  = anode(1,j)
	 a2  = anode(2,j)
c
	 dsj = dwjp(j)

	 xj  = cjp(1,j)
	 yj  = cjp(2,j)
	 r10 = (xj - xpt)/a1
	 r01 = (yj - ypt)/a2
	 r20 = r10*r10
	 r11 = r10*r01
	 r02 = r01*r01
c
	 r21 = r20*r01
	 r12 = r10*r02
	 r30 = r20*r10
	 r03 = r01*r02
c
	 r40 = r30*r10
	 r31 = r30*r01
	 r22 = r20*r02
	 r13 = r10*r03
	 r04 = r01*r03
c
         xx = dabs(r10)
	 yy = dabs(r01)

	 if(xx.gt.2.0.or.yy.gt.2.0) go to 30
c
c The length 2.0 is only cubic spline window function.
c
         call window(aw,awdx,awdy,awdxx,awdxy,
     &               awdyy,xpt,ypt,xj,yj,a1,a2)
c
c
	 aw    = aw*dsj
c
         am00   = am00  + aw
         am10   = am10  + r10*aw
         am01   = am01  + r01*aw
	 am20   = am20  + r20*aw
         am11   = am11  + r11*aw
         am02   = am02  + r02*aw
c
         am21   = am21  + r21*aw
         am12   = am12  + r12*aw
         am30   = am30  + r30*aw
         am03   = am03  + r03*aw
         am22   = am22  + r22*aw
         am31   = am31  + r31*aw
         am13   = am13  + r13*aw
         am40   = am40  + r40*aw
         am04   = am04  + r04*aw
c
	 awdx = awdx*dsj
	 awdy = awdy*dsj
	 dx   = -1.00/a1
	 dy   = -1.00/a2
c
         am00dx  = am00dx + awdx
	 am10dx  = am10dx + dx*aw + r10*awdx
	 am01dx  = am01dx + r01*awdx
	 am20dx  = am20dx + 2.0*dx*r10*aw + r20*awdx
	 am11dx  = am11dx + dx*r01*aw + r11*awdx
	 am02dx  = am02dx + r02*awdx
c
	 am21dx  = am21dx + 2.0*dx*r11*aw + r21*awdx
         am12dx  = am12dx + dx*r02*aw + r12*awdx
         am30dx  = am30dx + 3.0*dx*r20*aw + r30*awdx
	 am03dx  = am03dx + r03*awdx
	 am22dx  = am22dx + 2.0*dx*r12*aw + r22*awdx
         am31dx  = am31dx + 3.0*dx*r21*aw + r31*awdx
	 am13dx  = am13dx + dx*r03*aw + r13*awdx
	 am40dx  = am40dx + 4.0*dx*r30*aw + r40*awdx
         am04dx  = am04dx + r04*awdx
c
         am00dy  = am00dy  + awdy
	 am10dy  = am10dy  + r10*awdy
	 am01dy  = am01dy  + dy*aw + r01*awdy
	 am20dy  = am20dy  + r20*awdy
	 am11dy  = am11dy  + r10*dy*aw + r11*awdy
	 am02dy  = am02dy  + 2.0*r01*dy*aw + r02*awdy
c
	 am21dy  = am21dy + dy*r20*aw + r21*awdy
         am12dy  = am12dy + 2.0*dy*r11*aw + r12*awdy
         am30dy  = am30dy + r30*awdy
	 am03dy  = am03dy + 3.0*dy*r02*aw + r03*awdy
	 am22dy  = am22dy + 2.0*dy*r21*aw + r22*awdy
         am31dy  = am31dy + dy*r30*aw + r31*awdy
	 am13dy  = am13dy + 3.0*dy*r12*aw + r13*awdy
	 am40dy  = am40dy + r40*awdy
         am04dy  = am04dy + 4.0*dy*r03*aw + r04*awdy
c
c........the second derivatives
c
	 awdxx = awdxx*dsj
	 awdxy = awdxy*dsj
	 awdyy = awdyy*dsj
         dxx   = dx*dx
	 dxy   = dx*dy
	 dyy   = dy*dy
c
c
	 am00dxx = am00dxx + awdxx
	 am10dxx = am10dxx + 2.0*dx*awdx + r10*awdxx
         am01dxx = am01dxx + r01*awdxx
	 am20dxx = am20dxx + 2.0*dxx*aw + 4.0d0*dx*r10*awdx 
     &           + r20*awdxx
	 am11dxx = am11dxx + 2.0*dx*r01*awdx + r11*awdxx
	 am02dxx = am02dxx + r02*awdxx
c
	 am21dxx = am21dxx + 2.0*dxx*r01*aw +4.0*dx*r11*awdx 
     &           + r21*awdxx
         am12dxx = am12dxx + 2.0*dx*r02*awdx + r12*awdxx
         am30dxx = am30dxx + 6.0*dxx*r10*aw + 6.0*dx*r20*awdx 
     &           + r30*awdxx
	 am03dxx = am03dxx + r03*awdxx
	 am22dxx = am22dxx + 2.0*dxx*r02*aw + 4.0*dx*r12*awdx 
     &           + r22*awdxx
	 am31dxx = am31dxx + 6.0*dxx*r11*aw + 6.0*dx*r21*awdx 
     &           + r31*awdxx
         am13dxx = am13dxx + 2.0d0*dx*r03*awdx + r13*awdxx
	 am40dxx = am40dxx + 12.0*dxx*r20*aw + 8.0*dx*r30*awdx 
     &           + r40*awdxx
	 am04dxx = am04dxx + r04*awdxx
c
	 am00dxy = am00dxy + awdxy
	 am10dxy = am10dxy + dx*awdy + r10*awdxy
	 am01dxy = am01dxy + dy*awdx + r01*awdxy
	 am20dxy = am20dxy + 2.0*dx*r10*awdy + r20*awdxy
	 am11dxy = am11dxy + dxy*aw + r10*dy*awdx + dx*r01*awdy 
     &           + r11*awdxy
         am02dxy = am02dxy + 2.0*dy*r01*awdx + r02*awdxy
c
	 am21dxy = am21dxy + 2.0*dxy*r10*aw + r20*dy*awdx 
     &           + 2.0*dx*r11*awdy + r21*awdxy
	 am12dxy = am12dxy + 2.0*dxy*r01*aw + r02*dx*awdy 
     &           + 2.0*dy*r11*awdx + r12*awdxy
         am30dxy = am30dxy + 3.0*r20*dx*awdy + r30*awdxy
	 am03dxy = am03dxy + 3.0*r02*dy*awdx + r03*awdxy
	 am22dxy = am22dxy + 4.0*dxy*r11*aw + 2.0*r21*dy*awdx
     &           + 2.0*r12*dx*awdy + r22*awdxy
	 am31dxy = am31dxy + 3.0*r20*dxy*aw + 3.0*dx*r21*awdy
     &           + r30*dy*awdx + r31*awdxy
	 am13dxy = am13dxy + 3.0*r02*dxy*aw + 3.0*dy*r12*awdx
     &           + r03*dx*awdy + r13*awdxy
	 am40dxy = am40dxy + 4.0*r30*dx*awdy + r40*awdxy
	 am04dxy = am04dxy + 4.0*r03*dy*awdx + r04*awdxy
c
c
	 am00dyy = am00dyy + awdyy
	 am10dyy = am10dyy + r10*awdyy
         am01dyy = am01dyy + 2.0*dy*awdy + r01*awdyy
	 am20dyy = am20dyy + r20*awdyy
	 am11dyy = am11dyy + 2.0*r10*dy*awdy + r11*awdyy
	 am02dyy = am02dyy + 2.0*dyy*aw + 4.0*dy*r01*awdy 
     &           + r02*awdyy
c
	 am21dyy = am21dyy + 2.0*r20*dy*awdy + r21*awdyy
	 am12dyy = am12dyy + 2.0*r10*dyy*aw + 4.0*r11*dy*awdy 
     &           + r12*awdyy
	 am30dyy = am30dyy + r30*awdyy 
	 am03dyy = am03dyy + 6.0*dyy*r01*aw + 6.0*dy*r02*awdy 
     &           + r03*awdyy
	 am22dyy = am22dyy + 2.0*r20*dyy*aw 
     &           + 4.0*r21*dy*awdy +r22*awdyy
         am31dyy = am31dyy + 2.0*r30*dy*awdy + r31*awdyy
	 am13dyy = am13dyy + 6.0*dyy*r11*aw + 6.0*r12*dy*awdy 
     &           +  r13*awdyy
	 am40dyy = am40dyy + r40*awdyy
	 am04dyy = am04dyy + 12.0*dyy*r02*aw + 8.0*dy*r03*awdy 
     &           + r04*awdyy

  30  continue
c
c
c.....end of the main loop
c
c.....assemble the M matrix (gm(i,j))
c
      gm(1,1) = am00
      gm(1,2) = am10
      gm(1,3) = am01
      gm(1,4) = am20
      gm(1,5) = am11
      gm(1,6) = am02
c
      gm(2,1) = am10
      gm(2,2) = am20
      gm(2,3) = am11
      gm(2,4) = am30
      gm(2,5) = am21
      gm(2,6) = am12
c
      gm(3,1) = am01
      gm(3,2) = am11
      gm(3,3) = am02
      gm(3,4) = am21
      gm(3,5) = am12
      gm(3,6) = am03
c
      gm(4,1) = am20
      gm(4,2) = am30
      gm(4,3) = am21
      gm(4,4) = am40
      gm(4,5) = am31
      gm(4,6) = am22
c
      gm(5,1) = am11
      gm(5,2) = am21
      gm(5,3) = am12
      gm(5,4) = am31
      gm(5,5) = am22
      gm(5,6) = am13
c
      gm(6,1) = am02
      gm(6,2) = am12
      gm(6,3) = am03
      gm(6,4) = am22
      gm(6,5) = am13
      gm(6,6) = am04
c
c
c
c.....find the inverse
c
       n  = 6
       np = 6
       call gjinv(gm,gminv,n,np,flag)
c
c      do i = 1, 6
c         print *, (gminv(i,j),j=1,6)
c      enddo
c      stop
c
      det = am00*gminv(1,1)- am10*gminv(1,2)+ am01*gminv(1,3)
     &    - am20*gminv(1,4)+ am11*gminv(1,5)- am02*gminv(1,6)
c
c.....test convergence criteria
c
      if(det .le. zero) then
	print *, 'det =', det
	print *, 'STOP! the determinat det < 0 '
	stop
      else
      end if
c
c
c.....calculate the derivative of gm
c
      gmdx(1,1) = am00dx
      gmdx(1,2) = am10dx
      gmdx(1,3) = am01dx
      gmdx(1,4) = am20dx
      gmdx(1,5) = am11dx
      gmdx(1,6) = am02dx
c
      gmdx(2,2) = am20dx
      gmdx(2,3) = am11dx
      gmdx(2,4) = am30dx
      gmdx(2,5) = am21dx
      gmdx(2,6) = am12dx
c
      gmdx(3,3) = am02dx
      gmdx(3,4) = am21dx
      gmdx(3,5) = am12dx
      gmdx(3,6) = am03dx
c
      gmdx(4,4) = am40dx
      gmdx(4,5) = am31dx
      gmdx(4,6) = am22dx
c
      gmdx(5,5) = am22dx
      gmdx(5,6) = am13dx
c
      gmdx(6,6) = am04dx
c
c
      gmdy(1,1) = am00dy
      gmdy(1,2) = am10dy
      gmdy(1,3) = am01dy
      gmdy(1,4) = am20dy
      gmdy(1,5) = am11dy
      gmdy(1,6) = am02dy
c
      gmdy(2,2) = am20dy
      gmdy(2,3) = am11dy
      gmdy(2,4) = am30dy
      gmdy(2,5) = am21dy
      gmdy(2,6) = am12dy
c
      gmdy(3,3) = am02dy
      gmdy(3,4) = am21dy
      gmdy(3,5) = am12dy
      gmdy(3,6) = am03dy
c
      gmdy(4,4) = am40dy
      gmdy(4,5) = am31dy
      gmdy(4,6) = am22dy
c
      gmdy(5,5) = am22dy
      gmdy(5,6) = am13dy
c
      gmdy(6,6) = am04dy
c
c.....calculate the second derivative of gm
c 
      gmdxx(1,1) = am00dxx
      gmdxx(1,2) = am10dxx
      gmdxx(1,3) = am01dxx
      gmdxx(1,4) = am20dxx
      gmdxx(1,5) = am11dxx
      gmdxx(1,6) = am02dxx
c
      gmdxx(2,2) = am20dxx
      gmdxx(2,3) = am11dxx
      gmdxx(2,4) = am30dxx
      gmdxx(2,5) = am21dxx
      gmdxx(2,6) = am12dxx
c
      gmdxx(3,3) = am02dxx
      gmdxx(3,4) = am21dxx
      gmdxx(3,5) = am12dxx
      gmdxx(3,6) = am03dxx
c
      gmdxx(4,4) = am40dxx
      gmdxx(4,5) = am31dxx
      gmdxx(4,6) = am22dxx
c
      gmdxx(5,5) = am22dxx
      gmdxx(5,6) = am13dxx
c
      gmdxx(6,6) = am04dxx
c
      gmdxy(1,1) = am00dxy
      gmdxy(1,2) = am10dxy
      gmdxy(1,3) = am01dxy
      gmdxy(1,4) = am20dxy
      gmdxy(1,5) = am11dxy
      gmdxy(1,6) = am02dxy
c
      gmdxy(2,2) = am20dxy
      gmdxy(2,3) = am11dxy
      gmdxy(2,4) = am30dxy
      gmdxy(2,5) = am21dxy
      gmdxy(2,6) = am12dxy
c
      gmdxy(3,3) = am02dxy
      gmdxy(3,4) = am21dxy
      gmdxy(3,5) = am12dxy
      gmdxy(3,6) = am03dxy
c
      gmdxy(4,4) = am40dxy
      gmdxy(4,5) = am31dxy
      gmdxy(4,6) = am22dxy
c
      gmdxy(5,5) = am22dxy
      gmdxy(5,6) = am13dxy
c
      gmdxy(6,6) = am04dxy
c
      gmdyy(1,1) = am00dyy
      gmdyy(1,2) = am10dyy
      gmdyy(1,3) = am01dyy
      gmdyy(1,4) = am20dyy
      gmdyy(1,5) = am11dyy
      gmdyy(1,6) = am02dyy
c
      gmdyy(2,2) = am20dyy
      gmdyy(2,3) = am11dyy
      gmdyy(2,4) = am30dyy
      gmdyy(2,5) = am21dyy
      gmdyy(2,6) = am12dyy
c
      gmdyy(3,3) = am02dyy
      gmdyy(3,4) = am21dyy
      gmdyy(3,5) = am12dyy
      gmdyy(3,6) = am03dyy
c
      gmdyy(4,4) = am40dyy
      gmdyy(4,5) = am31dyy
      gmdyy(4,6) = am22dyy
c
      gmdyy(5,5) = am22dyy
      gmdyy(5,6) = am13dyy
c
      gmdyy(6,6) = am04dyy
c
      do i = 1,6
	 do j = i,6
            gmdx(j,i)  = gmdx(i,j)
            gmdy(j,i)  = gmdy(i,j)
            gmdxx(j,i) = gmdxx(i,j)
            gmdxy(j,i) = gmdxy(i,j)
            gmdyy(j,i) = gmdyy(i,j)
	 enddo
      enddo
c
c.....assign the value for b vector
c
      do i = 1, 6
         b(i)   = gminv(1,i)
         wx(i)  = gminv(2,i)
         wy(i)  = gminv(3,i)
         wxx(i) = gminv(4,i)
         wxy(i) = gminv(5,i)
         wyy(i) = gminv(6,i)
      enddo
c
c.....find the value for bd(2,3)
c
c.....compute M^(-1)dM/dx and M^(-1)dM/dy
c
      do i = 1,6
	 do j = 1,6
	    adx(i,j)  = 0.0    
	    ady(i,j)  = 0.0    
	    adxx(i,j) = 0.0    
	    adxy(i,j) = 0.0    
            adyy(i,j) = 0.0    
	    do k = 1,6
	       adx(i,j) = adx(i,j) + gminv(i,k)*gmdx(k,j)
	       ady(i,j) = ady(i,j) + gminv(i,k)*gmdy(k,j)
	       adxx(i,j) = adxx(i,j) + gminv(i,k)*gmdxx(k,j)
	       adxy(i,j) = adxy(i,j) + gminv(i,k)*gmdxy(k,j)
               adyy(i,j) = adyy(i,j) + gminv(i,k)*gmdyy(k,j)
            enddo
         enddo
      enddo
c
c
      do i = 1,6
	 bd(1,i)  = 0.0
	 bd(2,i)  = 0.0
c
	 wxd(1,i) = 0.0
	 wxd(2,i) = 0.0
	 wyd(1,i) = 0.0
	 wyd(2,i) = 0.0
c
	 wxxd(1,i) = 0.0
	 wxxd(2,i) = 0.0
	 wxyd(1,i) = 0.0
	 wxyd(2,i) = 0.0
	 wyyd(1,i) = 0.0
	 wyyd(2,i) = 0.0
         do j = 1,6
	    bd(1,i) = bd(1,i) - adx(i,j)*b(j)
	    bd(2,i) = bd(2,i) - ady(i,j)*b(j)
c
	    wxd(1,i) = wxd(1,i) - adx(i,j)*wx(j)
	    wxd(2,i) = wxd(2,i) - ady(i,j)*wx(j)
c
	    wyd(1,i) = wyd(1,i) - adx(i,j)*wy(j)
	    wyd(2,i) = wyd(2,i) - ady(i,j)*wy(j)
c
	    wxxd(1,i) = wxxd(1,i) - adx(i,j)*wxx(j)
	    wxxd(2,i) = wxxd(2,i) - ady(i,j)*wxx(j)
c
	    wxyd(1,i) = wxyd(1,i) - adx(i,j)*wxy(j)
	    wxyd(2,i) = wxyd(2,i) - ady(i,j)*wxy(j)
c
	    wyyd(1,i) = wyyd(1,i) - adx(i,j)*wyy(j)
	    wyyd(2,i) = wyyd(2,i) - ady(i,j)*wyy(j)
         enddo
      enddo
c
c
c....compute the 2nd derivatives bdd(3,3)
c
c
      do i = 1,6
	 bdd(1,i) = 0.00
	 bdd(2,i) = 0.00
	 bdd(3,i) = 0.00
c
	 wxdd(1,i) = 0.0 
	 wxdd(2,i) = 0.0 
	 wxdd(3,i) = 0.0 
c
	 wydd(1,i) = 0.0 
	 wydd(2,i) = 0.0 
	 wydd(3,i) = 0.0 
c
	 wxxdd(1,i) = 0.0 
	 wxxdd(2,i) = 0.0 
	 wxxdd(3,i) = 0.0 
c
	 wxydd(1,i) = 0.0 
	 wxydd(2,i) = 0.0 
	 wxydd(3,i) = 0.0 
c
	 wyydd(1,i) = 0.0 
	 wyydd(2,i) = 0.0 
	 wyydd(3,i) = 0.0 
c
	 do j = 1,6
	    bdd(1,i) = bdd(1,i) - 2.0*adx(i,j)*bd(1,j)
     &               - adxx(i,j)*b(j)
	    bdd(2,i) = bdd(2,i) - adx(i,j)*bd(2,j)
     &               - ady(i,j)*bd(1,j) - adxy(i,j)*b(j)
            bdd(3,i) = bdd(3,i) - 2.0*ady(i,j)*bd(2,j)
     &               - adyy(i,j)*b(j)
c
	    wxdd(1,i) = wxdd(1,i) - 2.0*adx(i,j)*wxd(1,j)
     &               - adxx(i,j)*wx(j)
	    wxdd(2,i) = wxdd(2,i) - adx(i,j)*wxd(2,j)
     &               - ady(i,j)*wxd(1,j) - adxy(i,j)*wx(j)
            wxdd(3,i) = wxdd(3,i) - 2.0*ady(i,j)*wxd(2,j)
     &               - adyy(i,j)*wx(j)
c
	    wydd(1,i) = wydd(1,i) - 2.0*adx(i,j)*wyd(1,j)
     &               - adxx(i,j)*wy(j)
	    wydd(2,i) = wydd(2,i) - adx(i,j)*wyd(2,j)
     &               - ady(i,j)*wyd(1,j) - adxy(i,j)*wy(j)
            wydd(3,i) = wydd(3,i) - 2.0*ady(i,j)*wyd(2,j)
     &               - adyy(i,j)*wy(j)
c
	    wxxdd(1,i) = wxxdd(1,i) - 2.0*adx(i,j)*wxxd(1,j)
     &               - adxx(i,j)*wxx(j)
	    wxxdd(2,i) = wxxdd(2,i) - adx(i,j)*wxxd(2,j)
     &               - ady(i,j)*wxxd(1,j) - adxy(i,j)*wxx(j)
            wxxdd(3,i) = wxxdd(3,i) - 2.0*ady(i,j)*wxxd(2,j)
     &               - adyy(i,j)*wxx(j)
c
	    wxydd(1,i) = wxydd(1,i) - 2.0*adx(i,j)*wxyd(1,j)
     &               - adxx(i,j)*wxy(j)
	    wxydd(2,i) = wxydd(2,i) - adx(i,j)*wxyd(2,j)
     &               - ady(i,j)*wxyd(1,j) - adxy(i,j)*wxy(j)
            wxydd(3,i) = wxydd(3,i) - 2.0*ady(i,j)*wxyd(2,j)
     &               - adyy(i,j)*wxy(j)
c
	    wyydd(1,i) = wyydd(1,i) - 2.0*adx(i,j)*wyyd(1,j)
     &               - adxx(i,j)*wyy(j)
	    wyydd(2,i) = wyydd(2,i) - adx(i,j)*wyyd(2,j)
     &                 - ady(i,j)*wyyd(1,j) - adxy(i,j)*wyy(j)
            wyydd(3,i) = wyydd(3,i) - 2.0*ady(i,j)*wyyd(2,j)
     &                 - adyy(i,j)*wyy(j)
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
