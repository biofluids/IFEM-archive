      program main
c**********************************************
c
c  This is a testing program for the shape function
c  of reproducing kernel method.
c
c
c
c
c**********************************************
      implicit double precision (a-h,o-z)
      parameter (mnode=10000)
c
      dimension cjp(3,mnode),dcjp(3,mnode),dwjp(mnode)
      dimension b(7),bd(3,7),dhp(3)
      dimension cpt(3),cjt(3),dcjt(3)
      dimension shpd(3),cpx(51),cpy(51),cpz(51)
c
      open(10,file='cn.dat')
c
      open(11,file='cndx.dat')
c
      open(12,file='cndy.dat')
c
      open(13,file='cndz.dat')
c
c.....input data
c     
      numnpx = 21
      numnpy = 21
      numnpz = 21
      numnp2 = numnpx * numnpy   
      numnp  = numnp2 * numnpz
c
c.....domain definition
c     
      x0 = - 1.00
      xn =   1.00
      y0 = - 1.00
      yn =   1.00
      z0 = - 1.00
      zn =   1.00
c
      dx = (xn-x0)/(numnpx-1)
      dy = (yn-y0)/(numnpy-1)
      dz = (zn-z0)/(numnpz-1)
c
      ax = 1.0
      ay = 1.0
      az = 1.0
c
c.....calculate the coordinates x(i) and y(j)
c.....and re-assigne the value for dcjp
c
      do i= 1, numnpx
	 cjp(1,i)  = x0 + (i-1)*dx
	 dcjp(1,i) = ax*dx
      enddo
c
      do i= numnpx+1, numnp2
	 cjp(1,i)  = cjp(1,i - numnpx)
	 dcjp(1,i) = ax*dx
      enddo
c
      do  j= 1, numnpx
	 cjp(2,j)  = y0 
	 dcjp(2,j) = ay*dy
      enddo
c
      do  j= numnpx+1, numnp2
	 cjp(2,j)  = cjp(2,j-numnpx) + dy
	 dcjp(2,j) = ay * dy
      enddo
c
      do j = 1, numnp2
	 cjp(3,j) = z0
	 dcjp(3,j) = az * dz
      enddo
c
      do k = numnp2+1, numnp
	 cjp(3,k)  = cjp(3,k-numnp2) + dz
	 dcjp(3,k) = az * dz
c
	 cjp(1,k)  = cjp(1,k-numnp2)
	 dcjp(1,k) = ax*dx
c       
	 cjp(2,k)  = cjp(2,k-numnp2) 
	 dcjp(2,k) = ay * dy
      enddo

c
c.....calculate the weights
c
      do i = 1, numnp
	 dwjp(i) = dx*dy*dz
	 if (cjp(1,i).eq.x0) dwjp(i) = 0.50*dwjp(i)
	 if (cjp(1,i).eq.xn) dwjp(i) = 0.50*dwjp(i)
	 if (cjp(2,i).eq.y0) dwjp(i) = 0.50*dwjp(i)
	 if (cjp(2,i).eq.yn) dwjp(i) = 0.50*dwjp(i) 
	 if (cjp(3,i).eq.z0) dwjp(i) = 0.50*dwjp(i)
	 if (cjp(3,i).eq.zn) dwjp(i) = 0.50*dwjp(i) 
      enddo
c
c.....
c
      npt    = 21
c
      xp1 = -0.20
      xp2 =  0.20
      yp1 = -0.20
      yp2 =  0.20
      zp1 = -0.20
      zp2 =  0.20
c
      dpx = (xp2-xp1)/(npt - 1)
      dpy = (yp2-yp1)/(npt - 1)
      dpz = (zp2-zp1)/(npt - 1)
c
      do i = 1, npt
	 cpx(i) = xp1 + (i-1)*dpx
	 cpy(i) = yp1 + (i-1)*dpy
	 cpz(i) = zp1 + (i-1)*dpz
      enddo
c
      do 60 i = 1, npt
	 cpt(1) = cpx(i)
	 do 65 j =1,npt
	    cpt(2) = cpy(j)
	    do 68 k = 1, npt
	       cpt(3) = cpz(k)
               call correct3dtl(b,bd,cpt,cjp,dcjp,dwjp,numnp)
c
               jt     = 4631
	       cjt(1) = cjp(1,jt)
	       cjt(2) = cjp(2,jt)
	       cjt(3) = cjp(3,jt)
	       wjt    = dwjp(jt)
	       dhp(1) = dcjp(1,jt)
	       dhp(2) = dcjp(2,jt)
	       dhp(3) = dcjp(3,jt)
c
               call shape3dtl(shp,shpd,b,bd,cpt,cjt,dhp,wjt)
c
	       write(10,9110) cpt(1),cpt(2),cpt(3),shp
	       write(11,9110) cpt(1),cpt(2),cpt(3),shpd(1)
  68         continue
  65	 continue
  60  continue 
c
 9110 format(4(x,e12.5))
 9220 format(51(x,e12.5))
c       
       close(10)
       close(11)
       close(12)
c
       stop
       end
c
