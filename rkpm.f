ccccccccccccccccccccccccccccccccccccccccccccccccccc
c rkpm.f
c G. Wagner 12/11/98
c - Driver to compute RKPM shape functions for
c   all local quadrature points
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rkpm(shrk,shrknode,cnn,ncnn,cnn2,ncnn2,ien,rng,xna)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c Inputs:
c
c   ien(nen,nec): finite element connectivies
c
c   xna(nsd,nn):  spatial locations for all nodes
c
c Outputs:
c
c   shrk(0:nsd,maxconn,nec*nquad): shape function values
c       and their first derivatives
c
c   cnn(maxconn,nqdc): quad point - node connectivities
c
c   ncnn(nqdc): number of influence nodes for each quad point
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      include "global.h"
      include "malloc.h"
      real* 8 shrk(0:nsd,maxconn,nec*nquad)
      real* 8 shrknode(maxconn,nnc)
      integer ien(nen,nec),rng(neface,nec),cnn(maxconn,nqdc),ncnn(nqdc)
      integer cnn2(maxconn,nnc),ncnn2(nnc)
      real* 8 xna(nsd,nn)
      real* 8 xqd(nsdpad,nquadpad), ssq(nenpad,nquadpad)
      real* 8 x(3), y(3), a(3), xn(nsdpad,nenpad)
      real* 8 xr(nsdpad,nsdpad), cf(nsdpad,nsdpad) 
      real* 8 b(4), bd(3,4)
      real* 8 shp, shpd(3), det, vol, coef
      real* 8 xmax, ymax, zmax
      integer nnum, qvcount, ncount, ierr, seed
      integer maxinf,mininf,avginf,nqvtot
      real* 8 starttime, endtime, inftime, comptime, nconntime
     
      real* 8 adist(nsd,nn), adistloc(nsd,nn)
      real* 8 dwjploc(nn), dwjp(nn)
      integer inf(maxconn),ninf
      pointer (adistptr,adist), (adistlocptr,adistloc), (dwjplocptr,dwjploc)
      pointer (dwjpptr,dwjp), (infptr,inf)

      integer nconn(maxconn,nn),nnconn(nn)
      pointer (nconnptr,nconn), (nnconnptr,nnconn)

      real *8 xnamin(nsdpad),xnamax(nsdpad),lcell(nsdpad)
      integer ncell(nsdpad),nx,ny,nz,AllocateStatus
      integer i,j,k

      TYPE List_Node
      INTEGER :: NodeNumber
      TYPE(List_Node), POINTER :: Next
      END TYPE List_Node

      TYPE List_Node_Ptr
      TYPE(List_Node), POINTER :: NodePtr
      END TYPE List_Node_Ptr

      TYPE(List_Node_Ptr), dimension(:,:,:), allocatable :: CellArray
      TYPE(List_Node), pointer :: TempPtr
      
      adistptr    = malloc(nn*nsd*fsize)
      adistlocptr = malloc(nn*nsd*fsize)
      dwjplocptr  = malloc(nn*fsize)
      dwjpptr     = malloc(nn*fsize)
      infptr      = malloc(maxconn*isize)

      nconnptr    = malloc(maxconn*nn*isize)
      nnconnptr   = malloc(nn*isize)

      coef = 0.5
      maxinf = 0
      mininf = 9999
      avginf = 0
      inftime = 0
      comptime = 0
      nconntime = 0


c  Calculate element coordinates
      call shape

c  Calculate nodal weights, in parallel
      dwjploc(:) = 0.0
      adistloc(:,:) = 0.0
      do ie = 1,nec
        do inl=1,nen
          do isd=1,nsd
            nnum = ien(inl,ie)
            xn(isd,inl) = xna(isd,nnum)
          end do
        end do

        xmax = coef*(maxval(xn(1,1:nen)) - minval(xn(1,1:nen)))
        ymax = coef*(maxval(xn(2,1:nen)) - minval(xn(2,1:nen)))
        zmax = coef*(maxval(xn(3,1:nen)) - minval(xn(3,1:nen)))
        
        do inl = 1,nen
          node = ien(inl,ie) 
          adistloc(1,node) = max(adistloc(1,node),xmax)
          adistloc(2,node) = max(adistloc(2,node),ymax)
          adistloc(3,node) = max(adistloc(3,node),zmax)
        end do

c  Calculate volume
        if (nen.eq.4) then
          include "vol3d4n.h"
        else
          include "vol3d8n.h"
        end if
        
        do inl = 1,nen
          nnum = ien(inl,ie)
          dwjploc(nnum) = dwjploc(nnum) + vol/nen
        end do
        
      end do
      
      call MPI_ALLREDUCE(dwjploc, dwjp, nn, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(adistloc, adist, nsd*nn, MPI_DOUBLE_PRECISION,
     &     MPI_MAX, MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
c  Compute nodal connectvities
      starttime = MPI_WTIME()

c  Calculate maximum range of geometry
      do isd = 1,nsd
        xnamin(isd) = minval(xna(isd,:))
        xnamax(isd) = maxval(xna(isd,:))
        lcell(isd) = 2.0*maxval(adist(isd,:))
        ncell(isd) = floor((xnamax(isd)-xnamin(isd))/lcell(isd))+1
      end do
        
c  Allocate memory for list
      allocate(CellArray(ncell(1),ncell(2),ncell(3)))

c  Loop over nodes
      do inl = 1,nn
        nx = max(ceiling((xna(1,inl)-xnamin(1))/lcell(1)),1)
        ny = max(ceiling((xna(2,inl)-xnamin(2))/lcell(2)),1)
        nz = max(ceiling((xna(3,inl)-xnamin(3))/lcell(3)),1)
        allocate(TempPtr, STAT=AllocateStatus)
        if (AllocateStatus /= 0) stop "*** Not enought memory ***"
        TempPtr%NodeNumber = inl
        TempPtr%Next => CellArray(nx,ny,nz)%NodePtr
        CellArray(nx,ny,nz)%NodePtr => TempPtr
      end do

      call getnconn(nconn,nnconn,ien)
      endtime = MPI_WTIME()
      nconntime = nconntime + (endtime-starttime)

      qvcount = 0
      ncount = 0
      cnn(:,:) = 0
      ncnn(:) = 0
      cnn2(:,:) = 0
      ncnn2(:) = 0

c Main element loop starts here
      do ie = 1,nec
        
        seed = ien(1,ie)

c  Calculate x value of quadrature points
        xqd(:,:) = 0
        do inl = 1,nen
          nnum = ien(inl,ie)
          do isd = 1,nsd
            xqd(isd,1:nquad)=xqd(isd,1:nquad)+xna(isd,nnum)*sq(0,inl,1:nquad)
          end do
        end do
        
c  Main quadrature point loop starts here
        do iq = 1,nquad
          qvcount = qvcount + 1
          x(1:nsd) = xqd(1:nsd,iq)
          
c  Get list of influence nodes
          starttime = MPI_WTIME()
          ninf=0
          inf=0
c          call newgetinf(inf,ninf,x,xna,adist,nconn,nnconn,seed)
          call getinf(inf,ninf,x,xna,adist,CellArray,xnamin,ncell,lcell)
          endtime = MPI_WTIME()
          inftime = inftime + (endtime - starttime)
          cnn(1:ninf,qvcount) = inf(1:ninf)
          ncnn(qvcount) = ninf
          if (ninf.gt.maxinf) maxinf = ninf
          if (ninf.lt.mininf) mininf = ninf
          avginf = avginf + ninf
          
          call correct3dl(b,bd,x,xna,adist,dwjp,nn,inf,ninf,maxconn)
          
c  loop over influence nodes and call RKPMshape
          do n = 1,ninf
            nnum = inf(n)
            do isd=1,nsd
              y(isd) = xna(isd,nnum)
              a(isd) = adist(isd,nnum)
            enddo
            starttime = MPI_WTIME()
            call RKPMshape3dl(shp,shpd,b,bd,x,y,a,dwjp(nnum))
            endtime = MPI_WTIME()
            comptime = comptime + (endtime - starttime)
            shrk(0,n,qvcount) = shp
            shrk(1:nsd,n,qvcount) = shpd(1:nsd)
          end do
        end do
        
c        if (myid.eq.0) then
c          if (mod(qvcount,10).eq.0) then
c            write (*,*) qvcount," qp's of ",nqdc
c          end if
c        end if
        
      end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


c Now calculate for nodes
      do inn = 1,nnc
        
        seed = myid*maxnnc + inn
        ncount = ncount + 1
        x(1:nsd) = xna(1:nsd, myid*maxnnc+inn)
          
c  Get list of influence nodes
        starttime = MPI_WTIME()
c        call newgetinf(inf,ninf,x,xna,adist,nconn,nnconn,seed)
        call getinf(inf,ninf,x,xna,adist,CellArray,xnamin,ncell,lcell)
        endtime = MPI_WTIME()
        inftime = inftime + (endtime - starttime)
        cnn2(1:ninf,ncount) = inf(1:ninf)
        ncnn2(ncount) = ninf
        if (ninf.gt.maxinf) maxinf = ninf
        if (ninf.lt.mininf) mininf = ninf

        call correct3dl(b,bd,x,xna,adist,dwjp,nn,inf,ninf,maxconn)
          
c  loop over influence nodes and call RKPMshape
        do n = 1,ninf
          nnum = inf(n)
          do isd=1,nsd
            y(isd) = xna(isd,nnum)
            a(isd) = adist(isd,nnum)
          enddo
          starttime = MPI_WTIME()
          call RKPMshape3dl(shp,shpd,b,bd,x,y,a,dwjp(nnum))
          endtime = MPI_WTIME()
          comptime = comptime + (endtime - starttime)
          shrknode(n,ncount) = shp
        end do

c        if (myid.eq.0) then
c          if (mod(ncount,10).eq.0) then
c            write (*,*) ncount," nodes of ",nnc
c          end if
c        end if
              
      end do
      
      
      call MPI_REDUCE(maxinf,maxinf,1,MPI_INTEGER,MPI_MAX,0,
     &     MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(mininf,mininf,1,MPI_INTEGER,MPI_MIN,0,
     &     MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(avginf,avginf,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(qvcount,nqvtot,1,MPI_INTEGER,MPI_SUM,0,
     &     MPI_COMM_WORLD,ierr)
      
      if (myid.eq.0) then
        avginf = avginf/nqvtot
        write(6,'(" Maximum Influence Nodes = ",i7)') maxinf
        write(6,'(" Minimum Influence Nodes = ",i7)') mininf
        write(6,'(" Average Influence Nodes = ",i7)') avginf
        write(6,'(" Influence nodes time    = ",f8.2)') inftime
        write(6,'(" Computation time        = ",f8.2)') comptime
        write(6,'(" Nodal connectivity time = ",f8.2)') nconntime
      end if
      
      call free(adistptr)
      call free(adistlocptr)
      call free(dwjpptr)
      call free(dwjplocptr)
      call free(infptr)
      call free(nconnptr)
      call free(nnconnptr)

c  Free memory used for Cell Array linked lists
      do i = 1,ncell(1)
        do j = 1,ncell(2)
          do k = 1,ncell(3)
            do
              if (.not.associated(CellArray(i,j,k)%NodePtr)) exit
              TempPtr => CellArray(i,j,k)%NodePtr
              CellArray(i,j,k)%NodePtr => TempPtr%Next
              deallocate(TempPtr)
            end do
          end do
        end do
      end do

      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c      subroutine getinf(inf,ninf,x,xna,adist)

c      include "global.h"

c      real*8 x(3), xna(nsd,nn), adist(nsd,nn)
c      real*8 r(nsdpad)
c      integer inf(maxconn)
c      integer ninf

c      ninf = 0
c      do i = 1,nn
c        r(1:nsd) = x(1:nsd) - xna(1:nsd,i)
c        if ((abs(r(1)).le.2*adist(1,i)).and.(abs(r(2)).le.2*adist(2,i)).and. 
c     &       (abs(r(3)).le.2*adist(3,i))) then
c          ninf = ninf + 1
c          inf(ninf) = i
c        end if
c      end do
c      if (ninf > maxconn) then
c        write (*,*) "Too many influence nodes!"
c        write (*,*) myid,ninf
c      else if (ninf.lt.4) then
c       write (*,*) "Not enough influence nodes!"
c       write (*,*) myid,ninf
c      end if


c      return
c      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine getinf(inf,ninf,x,xna,adist,CellArray,xnamin,ncell,lcell)

      include "global.h"

      real *8 x(3), xna(nsd,nn), adist(nsd,nn)
      real *8 r(nsdpad)
      integer inf(maxconn)
      integer ninf
      real *8 xnamin(nsdpad),lcell(nsdpad)
      integer ncell(nsdpad),nx,ny,nz
      integer i,j,k,node

      TYPE List_Node
      INTEGER :: NodeNumber
      TYPE(List_Node), POINTER :: Next
      END TYPE List_Node

      TYPE List_Node_Ptr
      TYPE(List_Node), POINTER :: NodePtr
      END TYPE List_Node_Ptr

      TYPE (List_Node_Ptr) CellArray(ncell(1),ncell(2),ncell(3))
      
      TYPE (List_Node), POINTER :: CurrPtr
      
c  Calculate nx, ny, nz
      nx = max(ceiling((x(1)-xnamin(1))/lcell(1)),1)
      ny = max(ceiling((x(2)-xnamin(2))/lcell(2)),1)
      nz = max(ceiling((x(3)-xnamin(3))/lcell(3)),1)
      
c  Loop over group of cells and test nodes in list
      ninf = 0
      do i = max(1,nx-1),min(ncell(1),nx+1)
        do j = max(1,ny-1),min(ncell(2),ny+1)
          do k = max(1,nz-1),min(ncell(3),nz+1)
            CurrPtr => CellArray(i,j,k)%NodePtr
            do
              if (.not.associated(CurrPtr)) exit
              node = CurrPtr%NodeNumber
              r(1:nsd) = x(1:nsd) - xna(1:nsd,node)
              if ((abs(r(1)).le.2*adist(1,node)).and.
     &             (abs(r(2)).le.2*adist(2,node)).and. 
     &             (abs(r(3)).le.2*adist(3,node))) then
                ninf = ninf + 1
                inf(ninf) = node
              end if
              CurrPtr => CurrPtr%Next
            end do
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine newgetinf(inf,ninf,x,xna,adist,nconn,nnconn,seed)

      include "global.h"
      include "malloc.h"

      real*8 x(3), xna(nsd,nn), adist(nsd,nn)
      real*8 r(nsdpad)
      integer inf(maxconn), ninf
      integer nconn(maxconn,nn), nnconn(nn), seed
      integer i,in,node,candidate

      integer checked(nn)
      pointer (checkedptr,checked)

      integer newinf1(maxconn),nnewinf1,newinf2(maxconn),nnewinf2
      pointer (newinf1ptr,newinf1), (newinf2ptr,newinf2)

      checkedptr = malloc(nn*isize)
      newinf1ptr = malloc(maxconn*isize)
      newinf2ptr = malloc(maxconn*isize)
      
      checked(:) = 0
      inf(:) = 0
      ninf = 0
      newinf1(:) = 0

      newinf1(1) = seed
      nnewinf1 = 1
      checked(seed) = 1

      do while (nnewinf1.gt.0)
        newinf2(:) = 0
        nnewinf2 = 0
        do i = 1,nnewinf1
          node = newinf1(i)
          do in = 1,nnconn(node)
            candidate = nconn(in,node)
            if (checked(candidate).eq.0) then
              r(1:nsd) = x(1:nsd) - xna(1:nsd,candidate)
              if ((abs(r(1)).le.2*adist(1,candidate)).and.
     &             (abs(r(2)).le.2*adist(2,candidate)).and. 
     &             (abs(r(3)).le.2*adist(3,candidate))) then
                nnewinf2 = nnewinf2 + 1
                newinf2(nnewinf2) = candidate
              end if
              checked(candidate) = 1
            end if
          end do
        end do
        inf(ninf+1:ninf+nnewinf1) = newinf1(1:nnewinf1)
        ninf = ninf + nnewinf1
        newinf1(:) = newinf2(:)
        nnewinf1 = nnewinf2
c        stop
      end do


      call free(checkedptr)
      call free(newinf1ptr)
      call free(newinf2ptr)

      return
      end

      

ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine testrkpm(shrk,cnn,ncnn)

      include "global.h"

      real*8 shrk(0:nsd,maxconn,nquad*nec)
      integer cnn(maxconn,nqdc), ncnn(nqdc)
      integer qvcount
      real*8 shp, shp1, sum, sum1

      qvcount = 0
      do ie = 1,nec
        do iq = 1,nquad
          qvcount = qvcount + 1
          sum = 0.0
          sum1 = 0.0
          do inl = 1,ncnn(iq)
            shp = shrk(0,inl,iq)
            shp1 = shrk(1,inl,iq)
            sum = sum + shp  
	    sum1 = sum1 + shp1
            if (shp1.gt.1E3) then
              write(6,*) qvcount, shp1
            end if
          end do
c	  write (6,*) myid,sum,sum1
        end do

      end do
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getnconn(nconn,nnconn,ien)

      include "global.h"
      include "malloc.h"

      integer nconn(maxconn,nn), nnconn(nn), ien(nen,nec)
      integer econn(2,12), neconn
      integer ie, icn, in, node1, node2
      integer ierr, status(MPI_STATUS_SIZE)
      integer iproc, numsend(nn), numrecv(nn), srce, dst
      logical unused

      pointer (numsendptr,numsend),(numrecvptr,numrecv)

      numsendptr = malloc(nn*isize)
      numrecvptr = malloc(nn*isize)

      nnconn(:) = 0

      if (nen.eq.4) then
        neconn = 6
        econn(1,1) = 1
        econn(2,1) = 2
        econn(1,2) = 1
        econn(2,2) = 3
        econn(1,3) = 1
        econn(2,3) = 4
        econn(1,4) = 2
        econn(2,4) = 3
        econn(1,5) = 2
        econn(2,5) = 4
        econn(1,6) = 3
        econn(2,6) = 4
      else if (nen.eq.8) then
        neconn = 12
        econn(1,1) = 1
        econn(2,1) = 2
        econn(1,2) = 2
        econn(2,2) = 3
        econn(1,3) = 3
        econn(2,3) = 4
        econn(1,4) = 4
        econn(2,4) = 1
        econn(1,5) = 5
        econn(2,5) = 6
        econn(1,6) = 6
        econn(2,6) = 7
        econn(1,1) = 7
        econn(2,1) = 8
        econn(1,2) = 8
        econn(2,2) = 5
        econn(1,3) = 1
        econn(2,3) = 5
        econn(1,4) = 2
        econn(2,4) = 6
        econn(1,5) = 3
        econn(2,5) = 7
        econn(1,6) = 4
        econn(2,6) = 8
      end if

      do ie = 1,nec
        do icn = 1,neconn
          node1 = ien(econn(1,icn),ie)
          node2 = ien(econn(2,icn),ie)
          unused = .true.
          do in = 1,nnconn(node1)
            if (nconn(in,node1).eq.node2) unused = .false.
          end do
          if (unused) then
            nnconn(node1) = nnconn(node1) + 1
            nconn(nnconn(node1),node1) = node2
          end if

          unused = .true.
          do in = 1,nnconn(node2)
            if (nconn(in,node2).eq.node1) unused = .false.
          end do
          if (unused) then
            nnconn(node2) = nnconn(node2) + 1
            nconn(nnconn(node2),node2) = node1
          end if
        end do
      end do

c Communication
      numsend(:) = nnconn(:)
      numrecv(:) = 0
      do iproc = 1,numproc-1
        srce = myid - iproc
        dst = myid + iproc
        if (srce.lt.0) srce = srce + numproc
        if (dst.gt.numproc-1) dst = dst - numproc
        call MPI_SEND(numsend,nn,MPI_INTEGER,dst,myid,
     &       MPI_COMM_WORLD,ierr)
        call MPI_RECV(numrecv,nn,MPI_INTEGER,srce,srce,
     &       MPI_COMM_WORLD,status,ierr)
        do in = 1,nn
          if (numsend(in).ne.0) then
            call MPI_SEND(nconn(1,in),numsend(in),MPI_INTEGER,dst,myid,
     &           MPI_COMM_WORLD,ierr)
          end if
          if (numrecv(in).ne.0) then
            call MPI_RECV(nconn(nnconn(in)+1,in),numrecv(in),MPI_INTEGER,
     &           srce,srce,MPI_COMM_WORLD,statur,ierr)
            nnconn(in) = nnconn(in) + numrecv(in)
          end if
        end do
      end do

      write(6,*) myid
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call free(numsendptr)
      call free(numrecvptr)

      return
      end


















