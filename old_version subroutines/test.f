	subroutine test(shrknode,cnn,ncnn,nn_solids,x_solids,
     +	xna,ien)

c	real*8 shrknode(90,801)
c	integer cnn(90,801),ncnn(801)

      include 'global.h'
c	parameter (mno2=801)

	!solids variables
	integer nn_solids
	real* 8 shrknode(maxconn,nn_solids)
      integer cnn(maxconn,nn_solids),ncnn(nn_solids)
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

c	shrknode(:,:)=1.0
c	cnn(:,:)=2
c	ncnn(:)=3
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

