CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine nodesetup(cnn, mn) 
	implicit none
	include "global.h"
      include "malloc.h"

	integer cnn(mn,nqdc), mn
        
      integer tmpw(1),tmpq(1),tmps(mn,1)
      pointer (tmpwptr,tmpw),(tmpqptr,tmpq),(tmpsptr,tmps)
	integer i,j,k,sum,g_proc,l_proc    
	integer isource,idestin,status(MPI_STATUS_SIZE),ierr
      integer maxon,maxal

           onptr = malloc(isize*(numproc-1))
       onfromptr = malloc(isize*(numproc-1))
           alptr = malloc(isize*(numproc-1))
       alfromptr = malloc(isize*(numproc-1))
       destinptr = malloc(isize*(numproc-1))
       sourceptr = malloc(isize*(numproc-1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Phase 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 tmpsptr = malloc(isize*nqdc*mn)
	 tmpwptr = malloc(isize*nn)

	 do i =1,nn
		tmpw(i)=0
	 enddo

	do i =1,nqdc
		 do j = 1, mn
			tmpw(cnn(j,i))=1
		 enddo
	enddo

	nn_on = 0
	do i = 1,nn
	        nn_on = nn_on + tmpw(i)
	enddo

      if(nn_on.gt.0) node_onptr = malloc(isize*nn_on)

      j = 0
      do i = 1,nn
          if(tmpw(i).eq.1) then
              j = j + 1
              node_on(j) = i
              tmpw(i) = j  
          endif
      enddo

      write(6,*) nn_on, myid

      do i=1,nqdc 
           do j = 1,mn 
                tmps(j,i) = 0
           enddo
	enddo

	do i=1,nqdc
		  do j = 1, mn
			 tmps(j,i) = tmpw(cnn(j,i))
		  enddo
	enddo

      do i=1,nqdc
            do j=1,mn 
                 cnn(j,i) = tmps(j,i)   
            enddo
      enddo

      call free(tmpwptr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       PHASE 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(nn_on.gt.0) then
	tmpwptr = malloc(isize*nn_on)
	tmpqptr = malloc(isize*nn_on)
	endif

        do i=1,numproc-1
              on(i) = 0
        enddo

        mmon = 0   

        do i=1,nn_on
              j = node_on(i)        !! al node_on number
              g_proc = (j-1)/maxnqdc
              l_proc = g_proc - myid
              if(l_proc.lt.0) l_proc = l_proc + numproc
              if(l_proc.eq.0) then
                       mmon = mmon + 1
              else
                       on(l_proc) = on(l_proc)+1
              endif
              tmpw(i) = l_proc
        enddo

      k = 0
      do i=0,numproc-1
        do j=1,nn_on
               if(tmpw(j).eq.i) then
                      k = k +1
                      tmpq(j) = k
               endif
        enddo
      enddo

      do i = 1,nqdc
             do j = 1,mn 
                tmps(j,i) = tmpq(cnn(j,i))
	       enddo
      enddo

	do i=1,nqdc
		  do j = 1, mn
			 cnn(j,i) = tmps(j,i)
		  enddo
	enddo

      do i=1,nn_on
             tmpw(tmpq(i)) = node_on(i)
      enddo 

        do i=1,nn_on
               j = tmpw(i)
               node_on(i) = mod(j-1,maxnqdc) + 1
        enddo

        sum = mmon + 1
        do i=1,numproc-1
        onfrom(i) = sum
        sum = sum + on(i)
        enddo

        call free(tmpwptr)
        call free(tmpqptr)
        call free(tmpsptr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       PHASE 3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i = 1,numproc-1
                isource = myid + i
                if(isource.gt.numproc-1) isource = isource - numproc
                source(i) = isource
                idestin = myid-i
                if(idestin.lt.0) idestin = numproc + idestin
                destin(i) = idestin
        enddo

	maxal = 0
	maxon = 0

	do i = 1,numproc-1
	idestin = source(i)
	isource = destin(i)
        call MPI_SEND(on(i),1,MPI_INTEGER,idestin,101,
     &                MPI_COMM_WORLD,ierr)
        call MPI_RECV(al(i),1,MPI_INTEGER,isource,101,
     &                MPI_COMM_WORLD,status,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	maxon = max(maxon,on(i))
	maxal = max(maxal,al(i))
	enddo

        if(maxon.gt.0) bufonptr = malloc(fsize*maxon)      
        if(maxal.gt.0) bufalptr = malloc(fsize*maxal)      

        sum = mmon + 1
        do i=1,numproc-1
           alfrom(i) = sum
           sum = sum + al(i)
        enddo
        nn_al = sum - 1          

	if(nn_al.gt.0) node_alptr = malloc(isize*nn_al)

        do i=1,mmon
              node_al(i) = node_on(i)
        enddo
	
	do i = 1,numproc-1
	idestin = source(i)
	isource = destin(i)
	j = onfrom(i)
	k = alfrom(i)
        call MPI_SEND(node_on(j),on(i),MPI_INTEGER,idestin,101,
     &                MPI_COMM_WORLD,ierr)
        call MPI_RECV(node_al(k),al(i),MPI_INTEGER,isource,101,
     &                MPI_COMM_WORLD,status,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
	enddo

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine grab_one(x,xn)
	implicit none
	include "global.h"
      include "malloc.h"
      integer status(MPI_STATUS_SIZE),ierr

	integer iproc,iptr,jptr,idestin,isource,i,j,node
	real*8 xn(nnc),x(nn_on)

	do i = 1,nn_on
	x(i) = 0.0
	enddo

C...OFF PN GATHER (STEP 1)
	do iproc = 1,numproc-1

		idestin = destin(iproc)
		isource = source(iproc)
		iptr = onfrom(iproc)
		jptr = alfrom(iproc)

	if(al(iproc).gt.0) then
        do j = 1,al(iproc)
        node = node_al(jptr+j-1)
        bufal(j) = xn(node)
	enddo
        call MPI_SEND(bufal,al(iproc),MPI_DOUBLE_PRECISION,idestin,
     &                101,MPI_COMM_WORLD,ierr)
	endif

	if(on(iproc).gt.0) then
        call MPI_RECV(bufon,on(iproc),MPI_DOUBLE_PRECISION,isource,
     &                101,MPI_COMM_WORLD,status,ierr)
	do i = 1,on(iproc)
	x(iptr+i-1) = bufon(i)
	enddo
	endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	enddo

C...ON PN GATHER(STEP 1) 

	do i = 1,mmon
	node = node_on(i)
	x(i) = xn(node)
	enddo

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine grab_all(x,xn,m)
	include "global.h"
  
      integer m,i
	real*8 xn(m,nnc),x(m,nn_on)
	real*8 yn(nnc),y(nn_on)

	do i = 1,m      
        yn(:) = xn(i,:)
        call grab_one (y,yn)
        x(i,:) = y(:)
	enddo
       
        return
        end


c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	updated.fcm                                                          c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine example(cnn,d,mdf)

	include "global.h"

	real* 8 dg(mdf,nn_on), d(mdf,nnc)
	integer cnn(maxconn,nqdc)
	
	dg(:,:) = 0.0
	call grab_all(dg,d,mdf)    

        do i=1,nqdc
	  do j=1,maxconn
             if(cnn(j,i).gt.0) then
	         node= cnn(j,i)
	         do k=1,mdf
		      known function = dg(k,node)
		 enddo
		endif
	   enddo
	enddo

	return
	end
