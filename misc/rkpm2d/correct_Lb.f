      subroutine correct_Lb(b,bd,bdd,cpt,cjp,dcjp,dwjp,nep,
     &                      n_support,Lmap)
c
c**************************************************************
c
c     This subroutine is to calculate the b vector and its
c     derivatives,for the moving least square reproducing kernel 
c     interpolant shape function.
c
c     This code is only offering b vector and its derivatives
c     up to second order for the 2-D case.
c     The interpolation function is choosen as
c     
c     p= ( 1, x , y )
c
c     Author : Shaofan Li;
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
c  i      dcjp(2,nep): the increment's of particle's coordinates
c                          a1(ip) =  dcjp(1,ip)
c                          a2(ip) =  dcjp(2,ip)
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
c   l     a1, a2:  dilation parameter; scalar;
c                  for uniformal mesh, we choose a1:=ax*dx
c                                                a2:=ay*dy
c                  roughly, 1 <ax,ay < 3
c
c                 we input the adjust constants ax, ay, and we construct
c                 the dilation parameter in such a way: 
c                 (for non-uniform mesh )
c
c                 a1 =  dcjp(1,j)
c                 a2 =  dcjp(2,j)
c
c    l   gm(3,3) : the basic moment matrix;
c
c    l   gminv(3,3): the inverse of the gm, i.e. gm^{-1}
c
c    l   gmdx(3,3): the derivatives of gm: d/dx (gm);
c
c    l   gmdy(3,3): the derivatives of gm: d/dy (gm);
c
c    l   gmdxx(3,3): d^2/dx^2 (gm)
c
c    l   gmdxy(3,3): d^2/dxy (gm)
c
c    l   gmdyy(3,3): d^2/dy^2 (gm)
c
c    l   adx(3,3):  M^(-1) d/dx M
c
c    l   ady(3,3):  M^(-1) d/dy M
c
c    l   adxx(3,3):  local array to compute b,xx
c
c    l   adxy(3,3):  local array to compute b,xy
c
c    l   adyy(3,3):  local array to compute b,yy
c
c  o     b(3)   :   the b vector, also the first row of gminv
c
c  o     bd(2,3):   the 1st derivatives of the b vector
c
c                   bd(1,3) := d/dx b
c                   bd(2,3) := d/dy b
c
c  o     bdd(3,3):  the 2nd derivatives of the b vector
c
c                   bdd(1,3) := d^2/dxx b
c                   bdd(2,3) := d^2/dxy b
c                   bdd(3,3) := d^2/dyy b
c
c     The subroutine is only designed for solving 2-D locally
c     linear interpolation moving field problem. i.e.
c
c    l   p0 := (1,0,0); 
c
c    l   p(3) := (1,x,y);
c
c    There two different subroutines are called by this subroutine:
c
c    cubic2d.f ----- the cartesian product type window function;
c  
c    quintw1b.f ----- the circular disc type window function
c
c**************************************************************
      implicit double precision (a-h,o-z)
      include 'parameter.h'
c
      integer n_support,j,jpt
      dimension Lmap(mnsch)
      dimension cjp(2,nep),cpt(2),dcjp(2,nep),dwjp(nep)
      dimension gminv(3,3),gmdx(3,3),gmdy(3,3)
      dimension gmdxx(3,3),gmdxy(3,3),gmdyy(3,3)
      dimension adx(3,3),ady(3,3)
      dimension adxx(3,3),adxy(3,3),adyy(3,3)
      dimension b(*),bd(2,*),bdd(3,*)
c
c.....set the initial value for moment and its derivatives:
c
      zero  = 0.00
c
      am00  = 0.00
      am10  = 0.00
      am01  = 0.00
      am20  = 0.00
      am11  = 0.00
      am02  = 0.00
c
c
      am00dx = 0.00
      am10dx = 0.00
      am01dx = 0.00
      am20dx = 0.00
      am11dx = 0.00
      am02dx = 0.00
c
c
      am00dy = 0.00
      am10dy = 0.00
      am01dy = 0.00
      am20dy = 0.00
      am11dy = 0.00
      am02dy = 0.00
c
c
      am00dxx = 0.00
      am10dxx = 0.00
      am01dxx = 0.00
      am20dxx = 0.00
      am11dxx = 0.00
      am02dxx = 0.00
c
c
      am00dxy = 0.00
      am10dxy = 0.00
      am01dxy = 0.00
      am20dxy = 0.00
      am11dxy = 0.00
      am02dxy = 0.00
c
c
      am00dyy = 0.00
      am10dyy = 0.00
      am01dyy = 0.00
      am20dyy = 0.00
      am11dyy = 0.00
      am02dyy = 0.00

c
c.....set the initial value for all array: 
c
      do i = 1, 3
	 do  j = 1, 3
	    gminv(i,j)  = 0.00
	    gmdx(i,j)   = 0.00
	    gmdy(i,j)   = 0.00
	    gmdxx(i,j)  = 0.00
	    gmdxy(i,j)  = 0.00
	    gmdyy(i,j)  = 0.00
         end do
	 b(i)  = 0.00
      end do  
c
c.....input the value cpt 
c
      xpt = cpt(1)
      ypt = cpt(2)
c
c.....main loop: calculate moment by Trapezodial rule
c
      do 30 j = 1, n_support
c
         jpt = Lmap(j)	 
c
c........define intermediate variable
c
c........a1 and a2
c
	 a1 =  dcjp(1,jpt)
	 a2 =  dcjp(2,jpt)
c
         dx = -1.0/a1
	 dy = -1.0/a2
	 dxx = dx * dx
	 dxy = dx * dy
	 dyy = dy * dy
c
	 dsj = dwjp(jpt)

	 xj  = cjp(1,jpt)
	 yj  = cjp(2,jpt)
	 r10 = (xj - xpt)/a1
	 r01 = (yj - ypt)/a2
c
         r20 = r10*r10
	 r11 = r10*r01
	 r02 = r01*r01
c
         xx = dabs(r10)
	 yy = dabs(r01)

	 if(xx.ge.2.0.or.yy.ge.2.0) go to 30
c
         call window(aw,awdx,awdy,awdxx,awdxy,
     &               awdyy,xpt,ypt,xj,yj,a1,a2)
c
         aw    = aw*dsj
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
c.....the 1st derivatives
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
c........the second derivatives
c
	 awdxx = awdxx*dsj
	 awdxy = awdxy*dsj
	 awdyy = awdyy*dsj
c
c
	 am00dxx = am00dxx + awdxx
	 am10dxx = am10dxx + 2.0*dx*awdx + r10*awdxx
	 am01dxx = am01dxx + r01*awdxx
         am20dxx = am20dxx + 2.0*dxx*aw + 2.0*dx*r10*awdx + r20*awdxx
	 am11dxx = am11dxx + 2.0*dx*r01*aw + r11*awdxx
         am02dxx = am02dxx + r02*awdxx
c
	 am00dxy = am00dxy + awdxy
	 am10dxy = am10dxy + dx*awdy + r10*awdxy
	 am01dxy = am01dxy + dy*awdx + r01*awdxy
         am20dxy = am20dxy + 2.0*dx*r10*awdy + r20*awdxy
	 am11dxy = am11dxy + dxy*aw + r10*dy*awdx 
     &           + r01*dx*awdy + r11*awdxy
         am02dxy = am02dxy + 2.0*dy*r01*awdx+ r02*awdxy
c
	 am00dyy = am00dyy + awdyy
	 am10dyy = am10dyy + r10*awdyy
	 am01dyy = am01dyy + 2.0*dy*awdy + r01*awdyy
         am20dyy = am20dyy + r20*awdyy
	 am11dyy = am11dyy + 2.0*dy*r10*awdy + r11*awdyy
         am02dyy = am02dyy + 2.0*dyy*aw + 2.0*dy*r01*awdy + r02*awdyy
c
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
       a11 =  am20*am02 - am11*am11
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
      if(det .le. zero) then
	print *, 'det =', det
	print *, 'STOP! the determinat det < 0 '
	stop
      else
      end if
c
c.....assemble the gminv(i,j)
c
      cdet = 1.0d0/det
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
c.....calculate the second derivative of gm
c 
      gmdxx(1,1)  =  am00dxx
      gmdxx(1,2)  =  am10dxx
      gmdxx(2,1)  =  am10dxx
      gmdxx(1,3)  =  am01dxx
      gmdxx(3,1)  =  am01dxx
      gmdxx(2,2)  =  am20dxx
      gmdxx(2,3)  =  am11dxx
      gmdxx(3,2)  =  am11dxx
      gmdxx(3,3)  =  am02dxx
c
      gmdxy(1,1)  =  am00dxy
      gmdxy(1,2)  =  am10dxy
      gmdxy(2,1)  =  am10dxy
      gmdxy(1,3)  =  am01dxy
      gmdxy(3,1)  =  am01dxy
      gmdxy(2,2)  =  am20dxy
      gmdxy(2,3)  =  am11dxy
      gmdxy(3,2)  =  am11dxy
      gmdxy(3,3)  =  am02dxy
c
      gmdyy(1,1)  =  am00dyy
      gmdyy(1,2)  =  am10dyy
      gmdyy(2,1)  =  am10dyy
      gmdyy(1,3)  =  am01dyy
      gmdyy(3,1)  =  am01dyy
      gmdyy(2,2)  =  am20dyy
      gmdyy(2,3)  =  am11dyy
      gmdyy(3,2)  =  am11dyy
      gmdyy(3,3)  =  am02dyy
c
c
c.....assign the value for b vector
c
      b(1) = gminv(1,1)
      b(2) = gminv(1,2)
      b(3) = gminv(1,3)
c
c.....find the value for bd(2,3)
c
c.....compute M^(-1)* M_x and M^(-1) * M_y
c
      do i = 1,3
	 do j = 1,3
	    adx(i,j)  = 0.0    
	    ady(i,j)  = 0.0    
	    adxx(i,j) = 0.0
	    adxy(i,j) = 0.0
	    adyy(i,j) = 0.0
	    do k = 1,3
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
      do i = 1,3
	 bd(1,i) = 0.0
	 bd(2,i) = 0.0
         do j = 1,3
	    bd(1,i) = bd(1,i) - adx(i,j)*b(j)
	    bd(2,i) = bd(2,i) - ady(i,j)*b(j)
         end do
      end do
c
c....compute the 2nd derivatives bdd(3,3)
c
      do i = 1,3
         bdd(1,i) = 0.0
         bdd(2,i) = 0.0
         bdd(3,i) = 0.0
	 do j = 1,3
	    bdd(1,i) = bdd(1,i) - 2.0*adx(i,j)*bd(1,j)
     &               - adxx(i,j)*b(j)
	    bdd(2,i) = bdd(2,i) - adx(i,j)*bd(2,j)
     &               - ady(i,j)*bd(1,j) - adxy(i,j)*b(j)
	    bdd(3,i) = bdd(3,i) - 2.0*ady(i,j)*bd(2,j)
     &               - adyy(i,j)*b(j)
         end do
      end do
c
c
  999 continue
c
      return
      end
c
