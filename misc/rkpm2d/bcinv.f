      subroutine bcinv(LmDispbcX,LmDispbcY,
     &                 LmVelbcX,LmVelbcY,
     &                 xm,dxm,dvm,
     &                 numnp)
c
c--------------------------------------------------------------
c
c
c   Note:
c
c   In DYME2D (Dynamic Meshless (FEM/RKPM)) Solid 2-D code:
c
c   Three types essential b.c. are specified:
c
c   (1) the purely displacement b.c., by which, we mean:
c
c       both the velocity and acceleration on the boundary
c       are zero.
c
c   (2) the purely velocity b.c. , by which, we mean:
c
c       the acceleration on the b.c. is zero.
c       and, the displacement on the b.c. is prescribed
c       in exact fashion
c
c   (3) the general velocity b.c., 
c
c         the acceleration on the b.c. is prescribed;
c         the velocity on the b.c. is prescribed;
c    and, the displacement on the b.c. is prescribed;
c       
c   To simplify the input process, only velocity b.c card
c   and displacement b.c. card are read.
c   
c   Thus, in bcdisp.f, the displacement b.c. on both Vbc
c   and Dbc are enforced.
c
c
c
c--------------------------------------------------------------
c
      implicit none
      include 'parameter.h'
c      
      integer i,n_support,Lmap(mnsch)
      integer nTract,nDispbcX,nDispbcY,
     &        nVelbcX,nVelbcY
      integer LmDispbcX(maxDispbcX),LmDispbcY(maxDispbcY),
     &        LmVelbcX(maxVelbcX),LmVelbcY(maxVelbcY)
c      
      real*8  xm(2,maxNumnp),dxm(2,maxNumnp),dvm(maxNumnp),
     &  shBCX(maxEssbcX,maxEssbcX),shBCY(maxEssbcY,maxEssbcY)
c      
c	! local vars
c
      integer LmEssbcX(maxEssbcX),LmEssbcY(maxEssbcY),
     &        ipvtX(maxEssbcX),ipvtY(maxEssbcY)
      real*8  shp,shpd(2),shpdd(3)
      real*8  workX(maxEssbcX),zX(maxEssbcX),
     &        workY(maxEssbcY),zY(maxEssbcY)
c           
      real*8 b(6),bd(2,6),bdd(3,6)
      real*8 cjt(2),dcjt(2),cpt(2),det(2)
      real*8 wjt,xr,yr,rcondX,rcondY,detM
      
      integer nEssbcX,nEssbcY,
     &        iEssbcX,jEssbcX,iEssbcY,jEssbcY,
     &        iDispbcX,iDispbcY,
     &        iVelbcX,iVelbcY,
     &        numnp,ip,jp
c
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      integer imethRKPM,imethFEM,imeth,iIntMeth,iInter
c
      common /shapeBC/shBCX,shBCY
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /meth/imethRKPM,imethFEM,imeth,iIntMeth,iInter
      common /bound/nTract,nDispbcX,nDispbcY,nVelbcX,nVelbcY
c
c
c Essential Boundary Condition : Consistent BC Method ( for RKPM )
c
c
c---------- below is for RKPM
c 
      iEssbcX=0
      iEssbcY=0
c      
      do iDispbcX = 1,nDispbcX
         iEssbcX               = iEssbcX+1
         LmEssbcX(iEssbcX)     = LmDispbcX(iDispbcX)
      enddo
c
      do iDispbcY = 1,nDispbcY
         iEssbcY               = iEssbcY+1
         LmEssbcY(iEssbcY)     = LmDispbcY(iDispbcY)
      enddo
c      
      do iVelbcX = 1,nVelbcX
         iEssbcX               = iEssbcX+1
         LmEssbcX(iEssbcX)     = LmVelbcX(iVelbcX)
      enddo
c
      do iVelbcY = 1,nVelbcY
         iEssbcY               = iEssbcY + 1
         LmEssbcY(iEssbcY)     = LmVelbcY(iVelbcY)
      enddo
c
      nEssbcX = nDispbcX + nVelbcX
      nEssbcY = nDispbcY + nVelbcY
c
c----- (1)  X-Direction ----------
c
      if (nEssbcX .gt. 0) then
c
         do iEssbcX=1,nEssbcX
            ip = LmEssbcX(iEssbcX)
            cpt(1) = xm(1,ip)
            cpt(2) = xm(2,ip)
c
	    n_support = Lnp(ip)
	    do i = 1, n_support
	       Lmap(i) = Lmnp(i,ip)
	    enddo
c
            call correct(b,bd,bdd,cpt,xm,dxm,dvm,numnp,iInter,
     &                   n_support,Lmap)
c
            do jEssbcX = 1, nEssbcX
               shBCX(iEssbcX,jEssbcX) = 0.0d0
               jp = LmEssbcX(jEssbcX)
               cjt(1)  = xm(1,jp)
               cjt(2)  = xm(2,jp)
               dcjt(1) = dxm(1,jp)
               dcjt(2) = dxm(2,jp)
               wjt     = dvm(jp)
               xr = abs((cjt(1)-cpt(1))/dcjt(1))
               yr = abs((cjt(2)-cpt(2))/dcjt(2))
               if (xr .ge. 2.0d0 .or. yr.ge.2.0d0) then
                           ! nothing
               else
                   call RKPMshape(shp,shpd,shpdd,
     &             b,bd,bdd,cpt,cjt,dcjt,wjt,iInter)
c
                   shBCX(iEssbcX,jEssbcX) = shp
               endif
            enddo       !jEssbcX
         enddo  ! iEssbcX
c
	 call dgeco(shBCX,maxEssbcX,nEssbcX,ipvtX,rcondX,zX)
         print *, 'rcondX-B.C.', rcondX
	 call dgedi(shBCX,maxEssbcX,nEssbcX,ipvtX,det,workX,11)
c
         detM = det(1)*10.0**det(2)
c
	 print *, 'detMX', detM
c
       endif
c      
c----- (2) Y- Direction ---------
c
c
      if (nEssbcY .gt. 0) then
c
         do iEssbcY=1,nEssbcY
            ip=LmEssbcY(iEssbcY)
            cpt(1) = xm(1,ip)
            cpt(2) = xm(2,ip)
c
            n_support = Lnp(ip)
	    do i = 1, n_support
	       Lmap(i) = Lmnp(i,ip)
	    enddo
c
            call correct(b,bd,bdd,cpt,xm,dxm,dvm,numnp,
     &                   iInter,n_support,Lmap)
c
            do jEssbcY = 1, nEssbcY
               shBCY(iEssbcY,jEssbcY) = 0.0d0
               jp = LmEssbcY(jEssbcY)
               cjt(1)  = xm(1,jp)
               cjt(2)  = xm(2,jp)
               dcjt(1) = dxm(1,jp)
               dcjt(2) = dxm(2,jp)
               wjt     = dvm(jp)
               xr = abs((cjt(1)-cpt(1))/dcjt(1))
               yr = abs((cjt(2)-cpt(2))/dcjt(2))
c
               if (xr.ge.2.0d0 .or. yr.ge.2.0d0) then
c.................... nothing .......................
               else
                   call RKPMshape(shp,shpd,shpdd,
     &             b,bd,bdd,cpt,cjt,dcjt,wjt,iInter)
c
                   shBCY(iEssbcY,jEssbcY) = shp
               endif
c
            enddo       !jEssbcY
         enddo  ! iEssbcY
c
c	 call GJinv2(abcY,shBCY,nEssbcY,maxEssbcY,flag)
c
	 call dgeco(shBCY,maxEssbcY,nEssbcY,ipvtY,rcondY,zY)
         print *, 'rcondY-B.C.', rcondY
	 call dgedi(shBCY,maxEssbcY,nEssbcY,ipvtY,det,workY,11)
c
         detM = det(1)*10.0**det(2)
c
	 print *, 'detMY', detM
c
      endif    ! nEssbcY >0
c
      return
      end
c
