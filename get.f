	  subroutine getnodebc(nodebc,nodebcon2,nodebcv,bnodesall,
	1	   ien,rng,hn,hm2,hloc)

	  implicit none
	  include "global.h"
	  include "malloc.h"

	  real* 8 nodebc(ndf,nnc),nodebcon2(ndf,nn_on2),nodebcv(ndf,nnc)
	  integer bnodesall(0:numproc-1,ndf)
	  integer ien(nen,nec), rng(neface,nec)
	  real* 8 hn(nnc),hm2(nn_on2),hloc(nn_loc)

	  real* 8 bcloc(ndf,nn_loc),bcvloc(ndf,nn_loc)
	  pointer (bclocptr,bcloc),(bcvlocptr,bcvloc)
	  integer inn,ie,ieface,irng,idf,inface,node,iproc
	  integer ierr
	  logical assemble
	  
	  bclocptr = malloc(ndf*nn_loc*fsize)
	  bcvlocptr = malloc(ndf*nn_loc*fsize)

	  bnodes(:) = 0
	  bnodesall(:,:) = 0
	  bcloc(:,:) = 0
	  bcvloc(:,:) = 0

	  do ie = 1,nec
		do ieface = 1,neface
		  irng = rng(ieface,ie)
		  if (irng.gt.0) then
			do idf = 1,ndf
			  if (bc(idf,irng).eq.1) then
				do inface = 1,nnface
				  node = ien(map(ieface,inface,etype),ie)
				  bcloc(idf,node) = 1
				  bcvloc(idf,node) = bv(idf,irng)
				end do
			  end if
			end do
		  end if
		end do
	  end do

	  assemble = .false.
	  call scatter(bcloc,nodebc,ndf,assemble,hn,hloc)
	  call scatter(bcvloc,nodebcv,ndf,assemble,hn,hloc)
	  call grab_all2(nodebcon2,nodebc,ndf,hn,hm2)

	  do idf = 1,ndf
		do inn = 1,nnc
		  if (nodebc(idf,inn).eq.1) bnodes(idf) = bnodes(idf) + 1
		end do
		bnodesall(myid,idf) = bnodes(idf)
		do iproc = 0,numproc-1
		  call MPI_BCAST(bnodesall(iproc,idf),1,MPI_INTEGER,
	1		   iproc,MPI_COMM_WORLD,ierr)
		end do
	  end do

		

	  call MPI_ALLREDUCE(bnodes,bnodestot,ndf,MPI_INTEGER,
	1	   MPI_SUM, MPI_COMM_WORLD, ierr)
	  
	  call free(bclocptr)
	  call free(bcvlocptr)

	  end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getun(un,dn,d2,shrknode,cnn2,ncnn2,hn,hm2)

      implicit none
      include "global.h"

      real* 8 un(ndf,nnc),dn(ndf,nnc),d2(ndf,nn_on2)
      real* 8 shrknode(maxconn,nnc)
      integer cnn2(maxconn,nnc),ncnn2(nnc)
      real* 8 hn(nnc),hm2(nn_on2)
      
      integer inn,inl,isd,node

      call grab_all2(d2,dn,ndf,hn,hm2)

      un(:,:) = 0

      do inn = 1,nnc
        do inl = 1,ncnn2(inn)
          node = cnn2(inl,inn)
          un(:,inn) = un(:,inn) + shrknode(inl,inn) * d2(:,node)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  subroutine dbartod(dbarn,dn,d,d2,blist,bvlist,bnodesall,
	1	   nodebcon2,homog,maxb,
	2	   shrknode,cnn2,ncnn2,hn,hm,hm2)

	  implicit none
	  include "global.h"

	  integer maxb
	  real* 8 dbarn(ndf,nnc), dn(ndf,nnc), d(ndf,nn_on), d2(ndf,nn_on2)
	  integer blist(maxb,ndf),bnodesall(0:numproc-1,ndf)
	  logical homog
	  real* 8 bvlist(maxb,ndf),shrknode(maxconn,nnc),nodebcon2(ndf,nn_on2)
	  integer cnn2(maxconn,nnc),ncnn2(nnc)
	  real* 8 hn(nnc),hm(nn_on),hm2(nn_on2)

	  integer idf,ib,node,inl,neighbor,iproc,counter,mtype,size,ub,lb,extent
	  integer ierr,status(MPI_STATUS_SIZE),tag


c  Communicate data needed by processors
	  dn(:,:) = dbarn(:,:)
      call grab_all2 (d2,dn,ndf,hn,hm2)
	  
c  Ammend nodes of dn corresponding to boundary conditions
	  do idf = 1,ndf
c		write (*,*) myid
c		stop

		if (bnodes(idf).gt.0) then
		  do ib = 1,bnodes(idf)
			node = blist(ib,idf)
			if (homog) then
			  dn(idf,node) = 0
			else
			  dn(idf,node) = bvlist(ib,idf)
			end if
			do inl = 1,ncnn2(node)
			  neighbor = cnn2(inl,node)
			  if (nodebcon2(idf,neighbor).eq.0) then
				dn(idf,node) = dn(idf,node)-shrknode(inl,node)*d2(idf,neighbor)
			  end if
			end do
		  end do

c  Send these values to root processor for each idf
		  tag = myid*10 + idf
		  call MPI_SEND(dn(idf,1),1,bctype(idf),idf-1,tag,
	1		   MPI_COMM_WORLD,ierr)
		end if

		if (myid.eq.idf-1) then
		  if (bnodestot(idf).gt.0) then
			counter = 0
			do iproc = 0,numproc-1
			  if (bnodesall(iproc,idf).gt.0) then
				tag = iproc*10 + idf
				call MPI_RECV(rhs(1+counter),bnodesall(iproc,idf),
	1				 MPI_DOUBLE_PRECISION,
	2				 iproc,tag,MPI_COMM_WORLD,status,ierr)
				counter = counter + bnodesall(iproc,idf)
			  end if
			end do
		  end if
		end if
	  end do

c  Solve on each root processor
	  if (myid.lt.ndf) then
		idf = myid+1
		if (bnodestot(idf).gt.0) then
		  mtype = 1
		  call ma28cd(order, amat, licn, icoln, ikeep, rhs, wgt, mtype)
c  Return values to original processors
		  counter = 0
		  do iproc = 0,numproc-1
			if (bnodesall(iproc,idf).gt.0) then
			  call MPI_SEND(rhs(1+counter),bnodesall(iproc,idf),
	1			   MPI_DOUBLE_PRECISION,iproc,idf,
	2			   MPI_COMM_WORLD,ierr)
			  counter = counter + bnodesall(iproc,idf)
			end if
		  end do
		end if
	  end if

c  Receive on each processor
	  do idf = 1,ndf
		if (bnodes(idf).gt.0) then
		  call MPI_RECV(dn(idf,1),1,bctype(idf),idf-1,idf,
	1		   MPI_COMM_WORLD,status,ierr)
		end if
	  end do

c  Recommunicate data needed by each processor
	  call grab_all(d,dn,ndf,hn,hm)

	  return
	  end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  subroutine rtorbar(rbar,rbarn,rbar2,rinc,blist,bvlist,bnodesall,
	1	   nodebcon2,maxb,shrknode,cnn2,ncnn2,hn,hm,hm2)

	  implicit none
	  include "global.h"

	  integer maxb
	  real* 8 rbarn(ndf,nnc), rbar(ndf,nn_on), rbar2(ndf,nn_on2)
	  real* 8 rinc(ndf,nnc)
	  integer blist(maxb,ndf),bnodesall(0:numproc-1,ndf)
	  real* 8 bvlist(maxb,ndf),nodebcon2(ndf,nn_on2),shrknode(maxconn,nnc)
	  integer cnn2(maxconn,nnc),ncnn2(nnc)
	  real* 8 hn(nnc),hm(nn_on),hm2(nn_on2)

	  integer idf,ib,node,inl,neighbor,iproc,counter,mtype
	  integer ierr,status(MPI_STATUS_SIZE)
	  logical assemble

c  Communicate data needed by processors
	  assemble = .true.
	  call send_all (rbar, rbarn, ndf, assemble, hn, hm)

c  Send values from each processor to root
	  do idf = 1,ndf
		if (bnodes(idf).gt.0) then
		  call MPI_SEND(rbarn(idf,1),1,bctype(idf),idf-1,idf,
	1		   MPI_COMM_WORLD,ierr)
		end if
	  end do

c  Receive on root processor and solve
	  if (myid.lt.ndf) then
		idf = myid+1
		if (bnodestot(idf).gt.0) then
		  counter = 0
		  do iproc = 0,numproc-1
			if (bnodesall(iproc,idf).gt.0) then
			  call MPI_RECV(rhs(1+counter),bnodesall(iproc,idf),
	1			   MPI_DOUBLE_PRECISION,
	2			   iproc,idf,MPI_COMM_WORLD,status,ierr)
			  counter = counter + bnodesall(iproc,idf)
			end if
		  end do
c  Solve (transpose)
		  mtype = 0
		  call ma28cd(order, amat, licn, icoln, ikeep, rhs, wgt, mtype)
c  Return values to original processors
		  counter = 0
		  do iproc = 0,numproc-1
			if (bnodesall(iproc,idf).gt.0) then
			  call MPI_SEND(rhs(1+counter),bnodesall(iproc,idf),
	1			   MPI_DOUBLE_PRECISION,iproc,idf,
	2			   MPI_COMM_WORLD,ierr)
			  counter = counter + bnodesall(iproc,idf)
			end if
		  end do
		end if
	  end if

c  Receive values of solution
	  rbar2(:,:) = 0
	  do idf = 1,ndf
		if (bnodes(idf).gt.0) then
		  call MPI_RECV(rbarn(idf,1),1,bctype(idf),idf-1,idf,
	1		   MPI_COMM_WORLD,status,ierr)
c  Ammend values of r for non-boundary neighbors
		  do ib = 1,bnodes(idf)
			node = blist(ib,idf)
			do inl = 1,ncnn2(node)
			  neighbor = cnn2(inl,node)
			  if (nodebcon2(idf,neighbor).eq.0) then
				rbar2(idf,neighbor) = rbar2(idf,neighbor)
	1				 - shrknode(inl,node)*rbarn(idf,node)
			  end if
			end do
		  end do
		end if
	  end do

c  Increment rbarn and set equal to 0 on boundaries
	  rinc(:,:) = 0
	  assemble = .true.
	  call send_all2(rbar2,rinc,ndf,assemble,hn,hm2)
	  rbarn(:,:) = rbarn(:,:) + rinc(:,:)
	  do idf = 1,ndf
		if (bnodes(idf).gt.0) then
		  do ib = 1,bnodes(idf)
			node = blist(ib,idf)
			rbarn(idf,node) = 0
		  end do
		end if
	  end do
	  
	  return
	  end





