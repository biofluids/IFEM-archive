c
      subroutine essbc(nebpx,nebpy,Lmebx,Lmeby,xm,dxm,dvm,
     &                 numnp,Lnp,Lmnp,disp,vel,acc,shpn,
     &                 time,dt)
c
c
      implicit double precision (a-h,o-z)
c
      include 'parameter.h'
c
      integer n_support,i,is
      dimension Lmap(mnsch)
      dimension gbcd(maxEssbc),fbc(maxEssbc)
      dimension Lmebx(maxEssbc),Lmeby(maxEssbc)
      dimension abc(maxEssbc,maxEssbc),atem(maxEssbc,maxEssbc)
      dimension Lnp(maxNumnp),Lmnp(maxNumnp,mnsch)
      dimension xm(2,maxNumnp),dxm(2,maxNumnp),dvm(maxNumnp)
      dimension disp(2,maxNumnp),vel(2,maxNumnp),acc(2,maxNumnp)
      dimension shpd(2),shpn(maxNumnp,mnsch)
      dimension b(3),bd(2,3),cjt(2),dcjt(2),cpt(2)
c
c     Essential Boundary Condition : Consistent BC Method
c
c-----x direction----------
c
      do ji = 1, nebpx
        gbcd(ji) = 0.0
      enddo
c
      do i = 1, nebpx
        ip = Lmebx(i)
        cpt(1) = xm(1,ip)
        cpt(2) = xm(2,ip)
c
	n_support = Lnp(ip)
	do is = 1, n_support
	   Lmap(is) = Lmnp(is,ip)
	enddo
c
        call correct(b,bd,cpt,xm,dxm,dvm,numnp,
     &               iInter,n_support,Lmap)
c
        do 210 j = 1, nebpx
          abc(i,j) = 0.0d0
          jp = Lmebx(j)
          cjt(1)  = xm(1,jp)
          cjt(2)  = xm(2,jp)
          dcjt(1) = dxm(1,jp)
          dcjt(2) = dxm(2,jp)
          wjt     = dvm(jp)
          xr = dabs((cjt(1)-cpt(1))/dcjt(1))
          yr = dabs((cjt(2)-cpt(2))/dcjt(2))
          if ( xr.ge.2.0.or.yr.ge.2.0) go to 210
          call RKPMshape(shp,shpd,b,bd,cpt,cjt,dcjt,wjt)
          abc(i,j) = shp
 210    continue
      enddo
c
      do ji = 1, nebpx
        ip = Lmebx(ji)
        disp(1,ip) = 0.0d0
      enddo
c
      do i = 1, nebpx
        iy = Lmebx(i)
        ftemp = 0.0d0
        iLoop = Lnp(iy)
        do k = 1, iLoop
          ky = Lmnp(iy,k)
          shbc = shpn(iy,k)
          ftemp = ftemp + shbc*disp(1,ky)
        enddo
        fbc(i) = gbcd(i) - ftemp
      enddo
c
      ifail = 0
      do i1 = 1, nebpx
      do j1 = 1, nebpx
        atem(i1,j1) = abc(i1,j1)
      enddo
      enddo
      call slnpd(nebpx,ifail,atem,fbc,mnpeb)
      do m12 = 1, nebpx
        m39 = Lmebx(m12)
        disp(1,m39) = fbc(m12)
        vel(1,m39) = 0.0d0
        acc(1,m39) = 0.0d0
      enddo
c      
c-----y direction---------
c
      do ji = 1, nebpy
        ip = Lmeby(ji)
        gbcd(ji) = 0.0
c        if (xm(2,ip).gt.9.9) gbcd(ji) = -5.0d2*(time+dt)
        if (xm(2,ip).gt.9.9*0.2) gbcd(ji) = -5.0d1*(time+dt)        
      enddo
c
      do i = 1, nebpy
        ip = Lmeby(i)
        cpt(1) = xm(1,ip)
        cpt(2) = xm(2,ip)
c
	n_support = Lnp(ip)
	do is = 1, n_support
	   Lmap(is) = Lmnp(is,ip)
	enddo
c
        call correct(b,bd,cpt,xm,dxm,dvm,numnp,
     &               n_support,Lmap)
c
        do 212 j = 1, nebpy
          abc(i,j) = 0.0d0
          jp = Lmeby(j)
          cjt(1)  = xm(1,jp)
          cjt(2)  = xm(2,jp)
          dcjt(1) = dxm(1,jp)
          dcjt(2) = dxm(2,jp)
          wjt     = dvm(jp)
          xr = dabs((cjt(1)-cpt(1))/dcjt(1))
          yr = dabs((cjt(2)-cpt(2))/dcjt(2))
          if ( xr.ge.2.0.or.yr.ge.2.0) go to 212
          call RKPMshape(shp,shpd,b,bd,cpt,cjt,dcjt,wjt)
          abc(i,j) = shp
 212    continue
      enddo
c
      do ji = 1, nebpy
        ip = Lmeby(ji)
        disp(2,ip) = 0.0d0
      enddo
c
      do i = 1, nebpy
        iy = Lmeby(i)
        ftemp = 0.0d0
        iLoop = Lnp(iy)
        do k = 1, iLoop
          ky = Lmnp(iy,k)
          shbc = shpn(iy,k)
          ftemp = ftemp + shbc*disp(2,ky)
        enddo
        fbc(i) = gbcd(i) - ftemp
      enddo
c
      ifail = 0
      do i1 = 1, nebpy
      do j1 = 1, nebpy
        atem(i1,j1) = abc(i1,j1)
      enddo
      enddo
      call slnpd(nebpy,ifail,atem,fbc,mnpeb)
c
      do m12 = 1, nebpy
        m39 = Lmeby(m12)
        disp(2,m39) = fbc(m12)
        vel(2,m39) = 0.0d0
        acc(2,m39) = 0.0d0
c        if (xm(2,m39).gt.9.9) vel(2,m39) = -5.0d2
        if (xm(2,m39).gt.9.9*0.2) vel(2,m39) = -5.0d1
        
      enddo
c      
c      print *,'ebc done'
c
      return
      end
c
