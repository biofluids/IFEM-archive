      subroutine trabc(LmTract,vTract,xm,dxm,dvm,
     &           fext,nTract,numnp,
     &           istep,nstep,iextf_LoadingType,ifwLoad)
c
c
      implicit none
      include 'parameter.h'
c
c----------------------------------------------------------
c
c  This is without wavelet adaptivity 
c
c
c
      integer nTract,numnp,n_support,Lmap(mnsch)
      integer LmTract(2,maxTract)
      real*8  vTract(2,2,maxTract)
      real*8  xm(2,maxNumnp),dxm(2,maxNumnp),dvm(maxNumnp)
      real*8  fext(2,maxNumnp)
      integer istep,nstep
      integer ifwLoad
      integer iextf_Loadingtype
c      
      real*8 b(6),bd(2,6),bdd(3,6)
c
      real*8 cpt(2),cjt(2),dcjt(2)
      real*8 zs(10),zw(10),shpd(2),shpdd(3)
c
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      integer ngbdry
      integer i,ip,ig,iTract,idim
      integer ip1,ip2,nnc1,nnc2
      real*8  vlength,x1,y1,x2,y2,fx1,fy1,fx2,fy2
      real*8  shp1,shp2,fx_g,fy_g
      real*8  wjt,xxr,yyr,shpi
c      
      real*8 vloadpar
c
c###
c
      real*8  ax,ay,afact1,afact2
      integer iInter,iAdapt,num_apt,iad_stage,Istage
      common /rkpm/ax,ay,afact1,afact2,nnc1,nnc2
      common /adapt/iInter,iAdapt,num_apt,iad_stage,Istage
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      do ip = 1,numnp
         do idim=1,2
            fext(idim,ip) = 0.0d0
         enddo
      enddo
c
c  Gaussian integration on the boundary
c
      ngbdry = 4   ! the number of gauss point on each edge
      call gauss1D(zs,zw,ngbdry)
c
      do iTract=1,nTract
         ip1=LmTract(1,iTract)
         ip2=LmTract(2,iTract)
c        
         x1 = xm(1,ip1)
         y1 = xm(2,ip1)
         x2 = xm(1,ip2)
         y2 = xm(2,ip2)
c         
         fx1=vTract(1,1,iTract)
         fy1=vTract(2,1,iTract)

         fx2=vTract(1,2,iTract)
         fy2=vTract(2,2,iTract)
        
         vlength=sqrt((x2-x1)**2+(y2-y1)**2)
         
         do ig = 1, ngbdry
         
            shp1=(1.-zs(ig))/2.
            shp2=(1.+zs(ig))/2.

            cpt(1) = shp1*x1+shp2*x2
            cpt(2) = shp1*y1+shp2*y2
            
            fx_g=shp1*fx1+shp2*fx2
            fy_g=shp1*fy1+shp2*fy2
c
c
	    n_support = Lnp(ip)
	    do i = 1, n_support
	       Lmap(i) = Lmnp(i,ip)
	    enddo
c
           call correct(b,bd,bdd,cpt,xm,dxm,
     &          dvm,numnp,iInter,n_support,
     &          Lmap)
c
c##
c
            do 6324 ip = 1, numnp
               cjt(1) = xm(1,ip)
               cjt(2) = xm(2,ip)
               dcjt(1) = dxm(1,ip)
               dcjt(2) = dxm(2,ip)
               wjt     = dvm(ip)
               xxr = dabs(cpt(1)-cjt(1))/dcjt(1)
               yyr = dabs(cpt(2)-cjt(2))/dcjt(2)
c
               if (xxr .gt. 2.0d0 .or. yyr .gt. 2.0d0) goto 6324
c
                  call RKPMshape(shpi,shpd,shpdd,
     &                 b,bd,bdd,cpt,cjt,dcjt,wjt,iInter)
c#
c#Axsi symmetry
c                            ! ????????????
c                  fext(1,ip) = fext(1,ip)
c     &                + fx_g*shpi*dabs(cpt(1))*vlength/2.0*zw(ig)
c                                  ^^^^^^^^^^^^^^^???????????  
c                  fext(2,ip) = fext(2,ip)
c     &                + fy_g*shpi*dabs(cpt(1))*vlength/2.0*zw(ig)
c

                fext(1,ip) = fext(1,ip)+fx_g*shpi*vlength/2.*zw(ig)
                fext(2,ip) = fext(2,ip)+fy_g*shpi*vlength/2.*zw(ig)
c
c
 6324       continue

         enddo	! ig
      enddo	! iTract
c
c      
      call extf_loadpar(istep,nstep,iextf_LoadingType,vloadpar)
c
      if ( 1 .eq. 0) then
         write(ifwLoad,'(1x,2i8,2e15.6)') istep,nstep,
     &                 real(istep)/real(nstep),vloadpar
      endif
c      
      do ip = 1,numnp
         do idim=1,2
            fext(idim,ip)=fext(idim,ip)*vloadpar
         enddo
      enddo
c      
      end	!ends trabc
c
