CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine commsetup(ien)
	  implicit none
	  include "global.h"
      include "malloc.h"

	  integer ien(nen,nec)
	  
      integer tmpw(1),tmpq(1),tmps(nen,1)
      pointer (tmpwptr,tmpw),(tmpqptr,tmpq),(tmpsptr,tmps)
	  integer ie,inl,i,j,k,sum,g_proc,l_proc    
	  integer ierr,isrc,ides,status(MPI_STATUS_SIZE)
      integer maxloc,maxglb

	  locptr = malloc(isize*(numproc-1))
	  locfromptr = malloc(isize*(numproc-1))
	  glbptr = malloc(isize*(numproc-1))
	  glbfromptr = malloc(isize*(numproc-1))
	  desptr = malloc(isize*(numproc-1))
	  srcptr = malloc(isize*(numproc-1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  PHASE 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  tmpsptr = malloc(isize*nec*nen)
	  tmpwptr = malloc(isize*nn)

	  do i =1,nn
		tmpw(i) = 0
	  enddo

	  do ie = 1,nec
		do inl = 1,nen
		  tmpw(ien(inl,ie)) =1
		enddo
	  enddo

	  nn_loc = 0
	  do i = 1,nn
		nn_loc = nn_loc + tmpw(i)
	  enddo

	  node_locptr = malloc(isize*nn_loc)

	  j = 0
	  do i = 1,nn
		if(tmpw(i).eq.1) then
		  j = j + 1
		  node_loc(j) = i
		  tmpw(i) = j  
		endif
	  enddo

	  do ie=1,nec
		do inl = 1,nen
		  tmps(inl,ie) = tmpw(ien(inl,ie))
		enddo
	  enddo 

	  do ie=1,nec
		do inl = 1,nen
		  ien(inl,ie) = tmps(inl,ie)
		enddo
	  enddo

	  call free(tmpwptr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  PHASE 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  tmpwptr = malloc(isize*nn_loc)
	  tmpqptr = malloc(isize*nn_loc)

	  do i=1,numproc-1
		loc(i) = 0
	  enddo

	  nnon = 0   

	  do i=1,nn_loc
		j = node_loc(i)			!! glb node_loc number
		g_proc = (j-1)/maxnnc
		l_proc = g_proc - myid
		if(l_proc.lt.0) l_proc = l_proc + numproc
		if(l_proc.eq.0) then
		  nnon = nnon + 1
		else
		  loc(l_proc) = loc(l_proc)+1
		endif
		tmpw(i) = l_proc
	  enddo

	  k = 0
	  do i=0,numproc-1
        do j=1,nn_loc
		  if(tmpw(j).eq.i) then
			k = k +1
			tmpq(j) = k
		  endif
        enddo
	  enddo

	  do ie = 1,nec
        do inl = 1,nen
		  tmps(inl,ie) = tmpq(ien(inl,ie))
        enddo
	  enddo

	  do ie=1,nec
        do inl = 1,nen
		  ien(inl,ie) = tmps(inl,ie)
        enddo
	  enddo

	  do i=1,nn_loc
		tmpw(tmpq(i)) = node_loc(i)
	  enddo 

	  do i=1,nn_loc
		j = tmpw(i)
		node_loc(i) = mod(j-1,maxnnc) + 1
	  enddo

	  sum = nnon + 1
	  do i=1,numproc-1
        locfrom(i) = sum
        sum = sum + loc(i)
	  enddo

	  call free(tmpwptr)
	  call free(tmpqptr)
	  call free(tmpsptr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  PHASE 3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  do i = 1,numproc-1
		isrc = myid + i
		if(isrc.gt.numproc-1) isrc = isrc - numproc
		src(i) = isrc
		ides = myid-i
		if(ides.lt.0) ides = numproc + ides
		des(i) = ides
	  enddo

	  maxglb = 0
	  maxloc = 0

	  do i = 1,numproc-1
		ides = src(i)
		isrc = des(i)
        call MPI_SEND(loc(i),1,MPI_INTEGER,ides,101,
	1		 MPI_COMM_WORLD,ierr)
        call MPI_RECV(glb(i),1,MPI_INTEGER,isrc,101,
	1		 MPI_COMM_WORLD,status,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
		maxloc = max(maxloc,loc(i))
		maxglb = max(maxglb,glb(i))
	  enddo

	  if(maxloc.ne.0) buflocptr = malloc(fsize*maxloc)      
	  if(maxglb.ne.0) bufglbptr = malloc(fsize*maxglb)      

	  sum = nnon + 1
	  do i=1,numproc-1
		glbfrom(i) = sum
		sum = sum + glb(i)
	  enddo
	  nn_glb = sum - 1          

	  node_glbptr = malloc(isize*nn_glb)

	  do i=1,nnon
		node_glb(i) = node_loc(i)
	  enddo
	  
	  do i = 1,numproc-1
		ides = src(i)
		isrc = des(i)
		j = locfrom(i)
		k = glbfrom(i)
        call MPI_SEND(node_loc(j),loc(i),MPI_INTEGER,ides,101,
	1		 MPI_COMM_WORLD,ierr)
        call MPI_RECV(node_glb(k),glb(i),MPI_INTEGER,isrc,101,
	1		 MPI_COMM_WORLD,status,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	  enddo

	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine gather_1(x,xn)
	  implicit none
	  include "global.h"
      include "malloc.h"
      integer ierr,status(MPI_STATUS_SIZE)

	  integer iproc,iptr,jptr,ides,isrc,i,j,node
	  real*8 xn(nnc),x(nn_loc)

	  do i = 1, nn_loc
		x(i) = 0.0
	  enddo

C...OFF PN GATHER (STEP 1)
	  do iproc = 1,numproc-1

		ides = des(iproc)
		isrc = src(iproc)
		iptr = locfrom(iproc)
		jptr = glbfrom(iproc)

		if(glb(iproc).gt.0) then
		  do j = 1,glb(iproc)
			node = node_glb(jptr+j-1)
			bufglb(j) = xn(node)
		  enddo
		  call MPI_SEND(bufglb,glb(iproc),MPI_DOUBLE_PRECISION,ides,
	1		   101,MPI_COMM_WORLD,ierr)
		endif

		if(loc(iproc).gt.0) then
		  call MPI_RECV(bufloc,loc(iproc),MPI_DOUBLE_PRECISION,isrc,
	1		   101,MPI_COMM_WORLD,status,ierr)
		  do i = 1,loc(iproc)
			x(iptr+i-1) = bufloc(i)
		  enddo
		endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  enddo

C...ON PN GATHER(STEP 1) 

	  do i = 1,nnon
		node = node_loc(i)
		x(i) = xn(node)
	  enddo

	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine scatter_1(x,xn,assemble)
	  implicit none
	  include "global.h"
      include "malloc.h"
      integer ierr,status(MPI_STATUS_SIZE)
	  integer iproc,iptr,jptr,ides,isrc,i,j,node
	  real*8 xn(nnc),x(nn_loc),flag
      logical assemble

	  flag = 0.0
	  if(assemble) flag = 1.0

	  do i = 1,nnc    
		xn(i) = 0.0
	  enddo

C...ON PN SCATTER(STEP1)
	  do i = 1,nnon 
		node = node_glb(i)
        xn(node) = x(i)
	  enddo

	  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C...OFF PN SCATTER(STEP2)

	  do iproc = 1,numproc-1

		ides = src(iproc)
		isrc = des(iproc)
		iptr = locfrom(iproc)
		jptr = glbfrom(iproc)

		if(loc(iproc).gt.0) then
		  do i = 1,loc(iproc)
			bufloc(i) = x(iptr+i-1)
		  enddo
		  call MPI_SEND(bufloc,loc(iproc),MPI_DOUBLE_PRECISION,ides,
	1		   101,MPI_COMM_WORLD,ierr)
		endif

		if(glb(iproc).gt.0) then
		  call MPI_RECV(bufglb,glb(iproc),MPI_DOUBLE_PRECISION,isrc,
	1		   101,MPI_COMM_WORLD,status,ierr)
		  do j = 1,glb(iproc)
			node = node_glb(jptr+j-1)
			if (flag.eq.0) then
			  xn(node) = max(xn(node),bufglb(j))
			else
			  xn(node) = xn(node)*flag + bufglb(j)
			end if
		  enddo
		endif
		
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  enddo

	  return
	  end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine gather(x,xn,m,yn,y)
	  include "global.h"
	  
      integer m,i,j
	  real*8 xn(m,nnc),x(m,nn_loc)
	  real*8 yn(nnc),y(nn_loc)

	  do i = 1,m      
		do j=1,nnc
		  yn(j) = xn(i,j)
		enddo
        call gather_1 (y,yn)
		do j=1,nn_loc
		  x(i,j) = y(j)
		enddo
	  enddo
	  
	  return
	  end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine scatter(x,xn,m,assemble,yn,y)
	  include "global.h"

      integer m,i,j
	  real*8 xn(m,nnc),x(m,nn_loc)
	  real*8 yn(nnc),y(nn_loc)
      logical assemble

	  do i = 1,m      
		do j=1,nn_loc
		  y(j) = x(i,j)
		enddo
        call scatter_1 (y,yn,assemble)
		do j=1,nnc
		  xn(i,j) = yn(j)
		enddo
	  enddo
	  
      return
      end







