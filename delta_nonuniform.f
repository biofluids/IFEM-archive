c Subroutine rkpm_delta
c Lucy Zhang
c 11/06/02
c Northwestern University

c This subroutine calculate the delta function using RKPM which is used for both
c interpolation and distribution of the velocities and forces respectively
c between fluids and solids domain.

      subroutine delta_nonuniform(shrknode,cnn,ncnn,nn_solids,x_solids,
     +     xna,ien,dwjp)

      include 'global.h'
	parameter (maxnn_solids=20000)

	!solids variables
	integer nn_solids
      real* 8 shrknode(maxconn,maxnn_solids)
      integer cnn(maxconn,maxnn_solids),ncnn(maxnn_solids)

      real* 8 x_solids(nsd,nn_solids)

	!fluids variables
	real* 8 xna(nsd,nn),xn(nsd,nen)
	integer ien(nen,ne),rng(neface,ne)
      real* 8 dwjp(nn), adist(nsd,nn)

	!local variables
      real* 8 x(3), y(3), a(3)
      real* 8 xr(nsd,nsd), cf(nsd,nsd) 
      real* 8 b(4), bd(3,4)
      real* 8 shp, shpd(3), det
	real* 8 xmax,ymax,zmax,vol
      real* 8 coef,avginf
      integer ncount, ierr
      integer maxinf,mininf,nmaxinf,nmininf,navginf,totinf
      integer ie,inl,isd,nnum
      integer inf(maxconn),ninf

C      coef = 0.5
      coef = 0.6d0
      maxinf = 0
      mininf = 9999
      avginf = 0
	cnn(:,:)=0
	ncnn(:)=0
	shrknode(:,:)=0.0d0

c  Calculate element coordinates
      call shape

c  Calculate nodal weights
      dwjp(:) = 0.0
      adist(:,:) = 0.0
      totinf = 0
      do ie = 1,ne
        do inl=1,nen
          do isd=1,nsd
            nnum = ien(inl,ie)
            xn(isd,inl) = xna(isd,nnum)
          enddo
        enddo

        xmax = coef*(maxval(xn(1,1:nen)) - minval(xn(1,1:nen)))
        ymax = coef*(maxval(xn(2,1:nen)) - minval(xn(2,1:nen)))
        zmax = coef*(maxval(xn(3,1:nen)) - minval(xn(3,1:nen)))
        
        do inl = 1,nen
          node = ien(inl,ie) 
          adist(1,node) = max(adist(1,node),xmax)
          adist(2,node) = max(adist(2,node),ymax)
          adist(3,node) = max(adist(3,node),zmax)
        enddo

c  Calculate volume
        if (nen.eq.4) then
          include "vol3d4n.h"
        else
          include "vol3d8n.h"
        endif
        
        do inl = 1,nen
          nnum = ien(inl,ie)
          dwjp(nnum) = dwjp(nnum) + vol/nen
        enddo
      enddo

c Calculate the RKPM shape function for the solids points
      do i = 1, nn_solids
         x(1:nsd)=x_solids(1:nsd,i) !get solids point coordinate
         ninf=0
	   inf(:)=0
c get a list of influence nodes from the fluids grid
         call getinf(inf,ninf,x,xna,adist,nn,nsd,maxconn)
         cnn(1:ninf,i)=inf(1:ninf)
         ncnn(i)=ninf

         if (ninf .gt. maxinf) maxinf = ninf
         if (ninf .lt. mininf) mininf = ninf
         totinf = totinf + ninf
c calculate the correction function
         call correct3d(b,bd,x,xna,adist,dwjp,nn,1,inf,ninf,maxconn)
         do n = 1, ninf
            nnum = inf(n)
            do isd = 1,nsd
               y(isd) = xna(isd,nnum)
               a(isd) = adist(isd,nnum)
            enddo
            call RKPMshape3d(shp,b,bd,x,y,a,dwjp(nnum))
            shrknode(n,i)=shp
         enddo
      enddo

      avginf = totinf/nn_solids
      write(6,'(" Maximum Influence Nodes = ",i7)') maxinf
      write(6,'(" Minimum Influence Nodes = ",i7)') mininf
      write(6,'(" Average Influence Nodes = ",f7.3)') avginf

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine finds the influence points of point x
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getinf(inf,ninf,x,xna,adist,nn,nsd,maxconn)

      real*8 x(3), xna(nsd,nn), adist(nsd,nn)
      real*8 r(nsd)
      integer inf(maxconn)
      integer ninf, i,nsd,nn,maxconn

c!!! MAKE SURE MAXCONN IS DEFINED IN COMMON.H
ccccccccccccccccccc
c   x = the coordinate of the point to be calculated for
c   xna = the coordinate of all points
c   r = the distance between point x and other points in the system
c   inf = a collection of all the influence points
c   ninf = total number of influence points
c   adist = the radial distance of the influence domain
ccccccccccccccccccc

      ninf = 0
      do i = 1,nn
        r(1:nsd) = x(1:nsd) - xna(1:nsd,i)
        if ((abs(r(1)).le.2*adist(1,i)).and.
     +       (abs(r(2)).le.2*adist(2,i)).and. 
     +       (abs(r(3)).le.2*adist(3,i))) then
          ninf = ninf + 1
          inf(ninf) = i
        endif
      enddo
	
      if (ninf > maxconn) then
        write (*,*) "Too many influence nodes!"
        write (*,*) ninf
      elseif (ninf.lt.4) then
       write (*,*) "Not enough influence nodes!"
       write (*,*) ninf
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine correct3d(b,bd,cpt,cjp,dcjp,
     &                   dwjp,nep,iInter,inf,ninf,maxconn)
c
c......3-D correct function......
c
      implicit none

	integer nep,iInter
	integer maxconn,ninf,inf(maxconn)
      real*8 b(*),bd(3,*),cpt(3)
      real*8 cjp(3,nep),dcjp(3,nep),dwjp(nep)
     

      if (iInter .eq. 1) then
	 call correct3dl(b,bd,cpt,cjp,dcjp,dwjp,nep,inf,ninf,maxconn)
      elseif(iInter .eq. 11)  then
	 call correct3dtl(b,bd,cpt,cjp,dcjp,dwjp,nep,inf,ninf,maxconn)
      else
	 print *, 'wrong iInter'
	 stop
      endif


      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

