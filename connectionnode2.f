CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine nodesetup2(cnn2, ncnn2, mn) 
	  implicit none
	  include "global.h"
      include "malloc.h"
	  
	  integer cnn2(mn,nnc), ncnn2(nnc), mn
	  
      integer tmpw(1),tmpq(1),tmps(mn,1)
      pointer (tmpwptr,tmpw),(tmpqptr,tmpq),(tmpsptr,tmps)
	  integer i,j,k,sum,g_proc,l_proc    
	  integer isource2,idestin2,status(MPI_STATUS_SIZE),ierr
      integer maxon,maxal

	  on2ptr = malloc(isize*(numproc-1))
	  onfrom2ptr = malloc(isize*(numproc-1))
	  al2ptr = malloc(isize*(numproc-1))
	  alfrom2ptr = malloc(isize*(numproc-1))
	  destin2ptr = malloc(isize*(numproc-1))
	  source2ptr = malloc(isize*(numproc-1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Phase 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  tmpsptr = malloc(isize*nnc*mn)
	  tmpwptr = malloc(isize*nn)
	  
	  do i =1,nn
		tmpw(i)=0
	  enddo
	  
	  do i =1,nnc
		do j = 1,ncnn2(i)
		  tmpw(cnn2(j,i))=1
		enddo
	  enddo

	  nn_on2 = 0
	  do i = 1,nn
		nn_on2 = nn_on2 + tmpw(i)
	  enddo

      if(nn_on2.gt.0) node_on2ptr = malloc(isize*nn_on2)
	  
      j = 0
      do i = 1,nn
		if(tmpw(i).eq.1) then
		  j = j + 1
		  node_on2(j) = i
		  tmpw(i) = j  
		endif
      enddo
	  
c      write(6,*) nn_on2, nnc, myid
	  
      do i=1,nnc 
		do j = 1,ncnn2(i)
		  tmps(j,i) = 0
		enddo
	  enddo
	  
	  do i=1,nnc
		do j = 1,ncnn2(i)
		  tmps(j,i) = tmpw(cnn2(j,i))
		enddo
	  enddo
	  
      do i=1,nnc
		do j=1,ncnn2(i) 
		  cnn2(j,i) = tmps(j,i)   
		enddo
      enddo

      call free(tmpwptr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  PHASE 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(nn_on2.gt.0) then
		tmpwptr = malloc(isize*nn_on2)
		tmpqptr = malloc(isize*nn_on2)
		almap2ptr = malloc(isize*nn_on2)
	  endif

	  do i=1,numproc-1
		on2(i) = 0
	  enddo

	  mmon2 = 0   
	  
	  do i=1,nn_on2
		j = node_on2(i)			!! al node_on2 number
		g_proc = (j-1)/maxnnc
		l_proc = g_proc - myid
		if(l_proc.lt.0) l_proc = l_proc + numproc
		if(l_proc.eq.0) then
		  mmon2 = mmon2 + 1
		else
		  on2(l_proc) = on2(l_proc)+1
		endif
		tmpw(i) = l_proc
	  enddo

      k = 0
      do i=0,numproc-1
        do j=1,nn_on2
		  if(tmpw(j).eq.i) then
			k = k +1
			tmpq(j) = k
		  endif
        enddo
      enddo
	  
      do i = 1,nnc
		do j = 1,ncnn2(i) 
		  tmps(j,i) = tmpq(cnn2(j,i))
		enddo
      enddo

	  do i=1,nnc
		do j = 1,ncnn2(i)
		  cnn2(j,i) = tmps(j,i)
		enddo
	  enddo

      do i=1,nn_on2
		tmpw(tmpq(i)) = node_on2(i)
		almap2(tmpq(i)) = node_on2(i)
      enddo 
	  
	  do i=1,nn_on2
		j = tmpw(i)
		node_on2(i) = mod(j-1,maxnnc) + 1
	  enddo

	  sum = mmon2 + 1
	  do i=1,numproc-1
        onfrom2(i) = sum
        sum = sum + on2(i)
	  enddo
	  
	  call free(tmpwptr)
	  call free(tmpqptr)
	  call free(tmpsptr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  PHASE 3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  do i = 1,numproc-1
		isource2 = myid + i
		if(isource2.gt.numproc-1) isource2 = isource2 - numproc
		source2(i) = isource2
		idestin2 = myid-i
		if(idestin2.lt.0) idestin2 = numproc + idestin2
		destin2(i) = idestin2
	  enddo

	  maxal = 0
	  maxon = 0

	  do i = 1,numproc-1
		idestin2 = source2(i)
		isource2 = destin2(i)
        call MPI_SEND(on2(i),1,MPI_INTEGER,idestin2,101,
	1		 MPI_COMM_WORLD,ierr)
        call MPI_RECV(al2(i),1,MPI_INTEGER,isource2,101,
	1		 MPI_COMM_WORLD,status,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
		maxon = max(maxon,on2(i))
		maxal = max(maxal,al2(i))
	  enddo

	  if(maxon.gt.0) bufon2ptr = malloc(fsize*maxon)      
	  if(maxal.gt.0) bufal2ptr = malloc(fsize*maxal)      

	  sum = mmon2 + 1
	  do i=1,numproc-1
		alfrom2(i) = sum
		sum = sum + al2(i)
	  enddo
	  nn_al2 = sum - 1          
	  
	  if(nn_al2.gt.0) node_al2ptr = malloc(isize*nn_al2)
	  
	  do i=1,mmon2
		node_al2(i) = node_on2(i)
	  enddo
	  
	  do i = 1,numproc-1
		idestin2 = source2(i)
		isource2 = destin2(i)
		j = onfrom2(i)
		k = alfrom2(i)
        call MPI_SEND(node_on2(j),on2(i),MPI_INTEGER,idestin2,101,
	1		 MPI_COMM_WORLD,ierr)
        call MPI_RECV(node_al2(k),al2(i),MPI_INTEGER,isource2,101,
	1		 MPI_COMM_WORLD,status,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	  enddo

	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine grab_one2(x,xn)
	  implicit none
	  include "global.h"
      include "malloc.h"
      integer status(MPI_STATUS_SIZE),ierr

	  integer iproc,iptr,jptr,idestin2,isource2,i,j,node
	  real*8 xn(nnc),x(nn_on2)

	  do i = 1,nn_on2
		x(i) = 0.0
	  enddo
          
C...OFF PN GATHER (STEP 1)
	  do iproc = 1,numproc-1
		
		idestin2 = destin2(iproc)
		isource2 = source2(iproc)
		iptr = onfrom2(iproc)
		jptr = alfrom2(iproc)
		
		if(al2(iproc).gt.0) then
		  do j = 1,al2(iproc)
			node = node_al2(jptr+j-1)
			bufal2(j) = xn(node)
		  enddo
		  call MPI_SEND(bufal2,al2(iproc),MPI_DOUBLE_PRECISION,idestin2,
	1		   101,MPI_COMM_WORLD,ierr)
		endif
		
		if(on2(iproc).gt.0) then
		  call MPI_RECV(bufon2,on2(iproc),MPI_DOUBLE_PRECISION,isource2,
	1		   101,MPI_COMM_WORLD,status,ierr)
		  do i = 1,on2(iproc)
			x(iptr+i-1) = bufon2(i)
		  enddo
		endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
	  enddo

C...ON PN GATHER(STEP 1) 

	  do i = 1,mmon2
		node = node_on2(i)
		x(i) = xn(node)
	  enddo

	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine send_one2(x,xn,assemble)
	  implicit none
	  include "global.h"
      include "malloc.h"
      integer ierr,status(MPI_STATUS_SIZE)
	  integer iproc,iptr,jptr,idestin2,isource2,i,j,node
	  real*8 xn(nnc),x(nn_on2),flag
      logical assemble

	  flag = 0.0
	  if(assemble) flag = 1.0

	  do i = 1,nnc    
		xn(i) = 0.0
	  enddo

C...ON PN SCATTER(STEP1)
	  do i = 1,mmon2 
		node = node_al2(i)
        xn(node) = x(i)
	  enddo
	  
	  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C...OFF PN SCATTER(STEP2)

	  do iproc = 1,numproc-1

		idestin2 = source2(iproc)
		isource2 = destin2(iproc)
		iptr = onfrom2(iproc)
		jptr = alfrom2(iproc)

		if(on2(iproc).gt.0) then
		  do i = 1,on2(iproc)
			bufon2(i) = x(iptr+i-1)
		  enddo
		  call MPI_SEND(bufon2,on2(iproc),MPI_DOUBLE_PRECISION,idestin2,
	1		   101,MPI_COMM_WORLD,ierr)
		endif

		if(al2(iproc).gt.0) then
		  call MPI_RECV(bufal2,al2(iproc),MPI_DOUBLE_PRECISION,isource2,
	1		   101,MPI_COMM_WORLD,status,ierr)
		  do j = 1,al2(iproc)
			node = node_al2(jptr+j-1)
			xn(node) = xn(node)*flag+bufal2(j)
		  enddo
		endif
		
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	  enddo

	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine grab_all2(x,xn,m,yn,y)
	  include "global.h"
	  
      integer m,i
	  real*8 xn(m,nnc),x(m,nn_on2)
	  real*8 yn(nnc),y(nn_on2)

	  do i = 1,m      
        yn(:) = xn(i,:)
        call grab_one2 (y,yn)
        x(i,:) = y(:)
	  enddo
	  
	  return
	  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  subroutine send_all2(x,xn,m,assemble,yn,y)
	  include "global.h"

      integer m,i,j
	  real*8 xn(m,nnc),x(m,nn_on2)
	  real*8 yn(nnc),y(nn_on2)
      logical assemble

	  do i = 1,m      
		do j=1,nn_on2
		  y(j) = x(i,j)
		enddo
        call send_one2 (y,yn,assemble)
		do j=1,nnc
		  xn(i,j) = yn(j)
		enddo
	  enddo
	  
      return
      end

c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  updated.fcm                                                          c
c  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  subroutine example2(cnn2,ncnn2,d,mdf)
	  
	  include "global.h"

	  real* 8 dg(mdf,nn_on2), d(mdf,nnc)
	  integer cnn2(maxconn,nnc), ncnn2(nnc)
	  
	  dg(:,:) = 0.0
	  call grab_all(dg,d,mdf)    
	  
	  do i=1,nnc
		do j=1,ncnn2(i)
		  node= cnn2(j,i)
		  do k=1,mdf
			known function = dg(k,node)
		  enddo
		enddo
	  enddo
	  
	  return
	  end


