
      subroutine bcsetup(bmap,blist,bvlist,ac,irnc,icnc,maxb,nodebc,
     &     nodebcon2,nodebcv,cnn2,ncnn2,shrknode)

      implicit none
      include "global.h"
      include "malloc.h"

      integer maxb
      integer bmap(nn,ndf),blist(maxb,ndf)
      real* 8 bvlist(maxb,ndf)
      real* 8 ac(maxb*maxconn,ndf)
      integer irnc(maxb*maxconn,ndf),icnc(maxb*maxconn,ndf)
      real* 8 nodebc(ndf,nnc),nodebcon2(ndf,nn_on2),nodebcv(ndf,nnc)
      integer cnn2(maxconn,nnc),ncnn2(nnc)
      real* 8 shrknode(maxconn,nnc)

      integer bcount,nbcount,idf,iproc,inn,getnnc
      integer ib,node,inl,row,neighbor
      integer counter
      real* 8 u, timer
      integer ierr,status(MPI_STATUS_SIZE),rowelemtype
      real* 8 bctemp(ndf,maxnnc)
      integer block(maxb),disp(maxb)
      pointer (bctempptr,bctemp),(blockptr,block),(dispptr,disp)

      bctempptr = malloc(ndf*maxnnc*fsize)
      if (maxb.gt.0) then
        blockptr  = malloc(maxb*isize)
        dispptr   = malloc(maxb*isize)
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Phase 1 -- set up node number map in bmap
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      do idf = 1,ndf
        if (myid.eq.idf-1) then
          bcount = 0
          nbcount = 0
          do iproc = 0,numproc-1            
            if (iproc.eq.myid) then
              do inn = 1,nnc
                if (nodebc(idf,inn).gt.0) then
                  bcount = bcount + 1
c                  bmap(bcount,idf) = myid*maxnnc + inn
                  bmap(myid*maxnnc + inn,idf) = bcount
                else
                  nbcount = nbcount + 1
c                  bmap(bnodestot(idf)+nbcount,idf) = myid*maxnnc + inn
                  bmap(myid*maxnnc + inn,idf) = bnodestot(idf)+nbcount
                end if
              end do
            else
              getnnc = (nn - 1) / numproc + 1
              if (iproc.eq.numproc-1) getnnc = nn - (numproc - 1) * maxnnc
              call MPI_RECV(bctemp(1,1),ndf*getnnc,MPI_DOUBLE_PRECISION,
     &             iproc,101,MPI_COMM_WORLD,status,ierr)
              do inn = 1,getnnc
                if (bctemp(idf,inn).gt.0) then
                  bcount = bcount + 1
c                  bmap(bcount,idf) = iproc*maxnnc + inn
                  bmap(iproc*maxnnc + inn,idf) = bcount
               else
                  nbcount = nbcount + 1
c                  bmap(bnodestot(idf)+nbcount,idf) = iproc*maxnnc + inn
                  bmap(iproc*maxnnc + inn,idf) = bnodestot(idf)+nbcount
                end if
              end do
            end if
          end do
        else
          call MPI_SEND(nodebc(1,1),ndf*nnc,MPI_DOUBLE_PRECISION,idf-1,101,
     &         MPI_COMM_WORLD,ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(bmap(1,idf),nn,MPI_INTEGER,idf-1,MPI_COMM_WORLD,ierr)

      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Phase 2 -- Form matrix A
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      do idf = 1,ndf
        bcount = 0
        do inn = 1,nnc
          if (nodebc(idf,inn).gt.0) then
            bcount = bcount + 1
            blist(bcount,idf) = inn
            bvlist(bcount,idf) = 0
            bvlist(bcount,idf) = nodebcv(idf,inn)
          end if
        end do
      end do

c  Look at nodebcon
c      do idf = 1,ndf
c        counter = 0
c        do inl = 1,nn_on2
c          if (nodebcon(idf,inl).gt.0) counter = counter + 1
c        end do
c        write (*,*) myid,idf,counter
c      end do
c      stop

      do idf = 1,ndf
        ecount(idf) = 0
        if (bnodes(idf).gt.0) then
          do ib = 1,bnodes(idf)
            node = blist(ib,idf)
            row = bmap(myid*maxnnc + node, idf)
            do inl = 1,ncnn2(node)
              neighbor = cnn2(inl,node)
              if (nodebcon2(idf,neighbor).gt.0) then
                ecount(idf) = ecount(idf) + 1
                ac(ecount(idf),idf)   = shrknode(inl,node)
                irnc(ecount(idf),idf) = row
                icnc(ecount(idf),idf) = bmap(almap2(neighbor), idf)
              end if
            end do
            
c            if (counter.eq.0) write (*,*) myid,idf,"flag!",node,ib,row

          end do
        end if
      end do

	  call MPI_ALLREDUCE(ecount,ecounttot,ndf,MPI_INTEGER,
     &     MPI_SUM, MPI_COMM_WORLD, ierr)

      if (myid.lt.ndf) then
        if (ecounttot(myid+1).gt.0) then
          order    = bnodestot(myid+1)
          nz       = ecounttot(myid+1)
          licn     = nz * 15
          lirn     = nz * 15
          amatptr  = malloc(licn*fsize)
          rhsptr   = malloc(order*fsize)
          icolnptr = malloc(licn*isize)
          irownptr = malloc(lirn*isize)
          ikeepptr = malloc(5*order*isize)
          iwptr    = malloc(8*order*isize)
          wgtptr   = malloc(order*isize)
        end if
      end if
      
      do idf = 1,ndf
        if (myid.eq.idf-1) then
          counter = 0
          do iproc = 0,numproc-1
            if (iproc.eq.myid) then
              do inn = 1,ecount(idf)
                counter = counter + 1
                amat(counter)  = ac(inn,idf)
                irown(counter) = irnc(inn,idf)
                icoln(counter) = icnc(inn,idf)
              end do
            else
              call MPI_RECV(getnnc,1,MPI_INTEGER,iproc,101,
     &             MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(amat(counter+1),getnnc,
     &             MPI_DOUBLE_PRECISION,
     &             iproc,102,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(irown(counter+1),getnnc,
     &             MPI_INTEGER,
     &             iproc,103,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(icoln(counter+1),getnnc,
     &             MPI_INTEGER,
     &             iproc,104,MPI_COMM_WORLD,status,ierr)

c              do inl = counter+1, counter+getnnc
c                if (irown(inl).eq.icoln(inl)) then
c                  if (amat(inl).eq.0) then
c                    write (*,*) myid,irown(inl),amat(inl)
c                  end if
c                end if
c              end do

              counter = counter + getnnc

              
            end if
          end do
          
c          write (*,*) myid, counter, ecounttot(myid+1), nz

        else          
          call MPI_SEND(ecount(idf),1,MPI_INTEGER,idf-1,101,
     &         MPI_COMM_WORLD,ierr)
          call MPI_SEND(ac(1,idf),ecount(idf),MPI_DOUBLE_PRECISION,idf-1,102,
     &         MPI_COMM_WORLD,ierr)
          call MPI_SEND(irnc(1,idf),ecount(idf),MPI_INTEGER,idf-1,103,
     &         MPI_COMM_WORLD,ierr)
          call MPI_SEND(icnc(1,idf),ecount(idf),MPI_INTEGER,idf-1,104,
     &         MPI_COMM_WORLD,ierr)
 
c          do inl = 1,ecount(idf)
c            if (irnc(inl,idf).eq.icnc(inl,idf)) then
c              if (ac(inl,idf).eq.0) then
c                write (*,*) myid,inl,idf,irnc(inl,idf),ac(inl,idf)
c              end if
c            end if
c          end do

       end if

      end do

 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Phase 3 -- LU Decomposition of matrix A
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (myid.lt.ndf) then
        idf = myid+1
c  Compute LU decomposition or read from file
        if (readdecomp) then
          call readlu
        else
          if (bnodestot(idf).gt.0) then
            u = 1.0
            timer = MPI_WTIME()
            call ma28ad(order,nz,amat,licn,irown,lirn,icoln,
     &           u,ikeep,iw,wgt,ierr)
            timer = MPI_WTIME() - timer
            write (*,*) myid, timer
          end if
          call writelu
        end if
      end if

c  Create new MPI types for easy communication
      call MPI_TYPE_VECTOR(1,1,ndf,MPI_DOUBLE_PRECISION,rowelemtype,ierr)
      do idf = 1,ndf
        if (bnodes(idf).gt.0) then
          do ib = 1,bnodes(idf)
            block(ib) = 1
c            disp(ib)  = blist(ib,idf)-1
            disp(ib)  = (blist(ib,idf)-1)*ndf
          end do
c          call MPI_TYPE_INDEXED(bnodes(idf),block,disp,rowelemtype,
c     &         bctype(idf),ierr)
          call MPI_TYPE_INDEXED(bnodes(idf),block,disp,MPI_DOUBLE_PRECISION,
     &         bctype(idf),ierr)
          call MPI_TYPE_COMMIT(bctype(idf),ierr)
        end if
      end do
      
      call free(bctempptr)
      call free(blockptr)
      call free(dispptr)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... for mesh update 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bcsetupd(dbmap,dblist,dbvlist,dac,dirnc,
     & dicnc,maxb,dnodebc,dnodebcon2,dnodebcv,cnn2,ncnn2,shrknode)

      implicit none
      include "global.h"
      include "malloc.h"

      integer maxb
      integer dbmap(nn,nsd),dblist(maxb,nsd)
      real* 8 dbvlist(maxb,nsd)
      real* 8 dac(maxb*maxconn,nsd)
      integer dirnc(maxb*maxconn,nsd),dicnc(maxb*maxconn,nsd)
      real* 8 dnodebc(nsd,nnc),dnodebcon2(nsd,nn_on2),dnodebcv(nsd,nnc)
      integer cnn2(maxconn,nnc),ncnn2(nnc)
      real* 8 shrknode(maxconn,nnc)

      integer bcount,nbcount,isd,iproc,inn,getnnc
      integer ib,node,inl,row,neighbor
      integer counter
      real* 8 u, timer
      integer ierr,status(MPI_STATUS_SIZE),rowelemtype
      real* 8 bctemp(nsd,maxnnc)
      integer block(maxb),disp(maxb)
      pointer (bctempptr,bctemp),(blockptr,block),(dispptr,disp)

      bctempptr = malloc(nsd*maxnnc*fsize)
      if (maxb.gt.0) then
        blockptr  = malloc(maxb*isize)
        dispptr   = malloc(maxb*isize)
      end if
      do isd = 1,nsd
        if (myid.eq.isd-1) then
          bcount = 0
          nbcount = 0
          do iproc = 0,numproc-1            
            if (iproc.eq.myid) then
              do inn = 1,nnc
                if (dnodebc(isd,inn).gt.0) then
                  bcount = bcount + 1
                  dbmap(myid*maxnnc + inn,isd) = bcount
                else
                  nbcount = nbcount + 1
                  dbmap(myid*maxnnc + inn,isd) = dbnodestot(isd)+nbcount
                end if
              end do
            else
              getnnc = (nn - 1) / numproc + 1
              if (iproc.eq.numproc-1) getnnc = nn - (numproc - 1) * maxnnc
              call MPI_RECV(bctemp(1,1),nsd*getnnc,MPI_DOUBLE_PRECISION,
     &             iproc,101,MPI_COMM_WORLD,status,ierr)
              do inn = 1,getnnc
                if (bctemp(isd,inn).gt.0) then
                  bcount = bcount + 1
                  dbmap(iproc*maxnnc + inn,isd) = bcount
               else
                  nbcount = nbcount + 1
                  dbmap(iproc*maxnnc + inn,isd) = dbnodestot(isd)+nbcount
                end if
              end do
            end if
          end do
        else
          call MPI_SEND(dnodebc(1,1),nsd*nnc,MPI_DOUBLE_PRECISION,isd-1,101,
     &         MPI_COMM_WORLD,ierr)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dbmap(1,isd),nn,MPI_INTEGER,isd-1,MPI_COMM_WORLD,ierr)

      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Phase 2 -- Form matrix A
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      do isd = 1,nsd
        bcount = 0
        do inn = 1,nnc
          if (dnodebc(isd,inn).gt.0) then
            bcount = bcount + 1
            dblist(bcount,isd) = inn
            dbvlist(bcount,isd) = 0
            dbvlist(bcount,isd) = dnodebcv(isd,inn)
          end if
        end do
      end do

      do isd = 1,nsd
        decount(isd) = 0
        if (dbnodes(isd).gt.0) then
          do ib = 1,dbnodes(isd)
            node = dblist(ib,isd)
            row = dbmap(myid*maxnnc + node, isd)
            do inl = 1,ncnn2(node)
              neighbor = cnn2(inl,node)
              if (dnodebcon2(isd,neighbor).gt.0) then
                decount(isd) = decount(isd) + 1
                dac(decount(isd),isd)   = shrknode(inl,node)
                dirnc(decount(isd),isd) = row
                dicnc(decount(isd),isd) = dbmap(almap2(neighbor), isd)
              end if
            end do
          end do
        end if
      end do

	  call MPI_ALLREDUCE(decount,decounttot,nsd,MPI_INTEGER,
     &     MPI_SUM, MPI_COMM_WORLD, ierr)

      if (myid.lt.nsd) then
        if (decounttot(myid+1).gt.0) then
          dorder    = dbnodestot(myid+1)
          dnz       = decounttot(myid+1)
          dlicn     = dnz * 15
          dlirn     = dnz * 15
          damatptr  = malloc(dlicn*fsize)
          drhsptr   = malloc(dorder*fsize)
          dicolnptr = malloc(dlicn*isize)
          dirownptr = malloc(dlirn*isize)
          dikeepptr = malloc(5*dorder*isize)
          diwptr    = malloc(8*dorder*isize)
          dwgtptr   = malloc(dorder*isize)
        end if
      end if
      
      do isd = 1,nsd
        if (myid.eq.isd-1) then
          counter = 0
          do iproc = 0,numproc-1
            if (iproc.eq.myid) then
              do inn = 1,decount(isd)
                counter = counter + 1
                damat(counter)  = dac(inn,isd)
                dirown(counter) = dirnc(inn,isd)
                dicoln(counter) = dicnc(inn,isd)
              end do
            else
              call MPI_RECV(getnnc,1,MPI_INTEGER,iproc,101,
     &             MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(damat(counter+1),getnnc,
     &             MPI_DOUBLE_PRECISION,
     &             iproc,102,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(dirown(counter+1),getnnc,
     &             MPI_INTEGER,
     &             iproc,103,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(dicoln(counter+1),getnnc,
     &             MPI_INTEGER,
     &             iproc,104,MPI_COMM_WORLD,status,ierr)

              counter = counter + getnnc

              
            end if
          end do

        else          
          call MPI_SEND(decount(isd),1,MPI_INTEGER,isd-1,101,
     &         MPI_COMM_WORLD,ierr)
          call MPI_SEND(dac(1,isd),decount(isd),MPI_DOUBLE_PRECISION,isd-1,102,
     &         MPI_COMM_WORLD,ierr)
          call MPI_SEND(dirnc(1,isd),decount(isd),MPI_INTEGER,isd-1,103,
     &         MPI_COMM_WORLD,ierr)
          call MPI_SEND(dicnc(1,isd),decount(isd),MPI_INTEGER,isd-1,104,
     &         MPI_COMM_WORLD,ierr)
 
       end if

      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Phase 3 -- LU Decomposition of matrix A
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (myid.lt.nsd) then
        isd = myid+1
c  Compute LU decomposition or read from file
        if (readdecomp) then
          call readlud
        else
          if (dbnodestot(isd).gt.0) then
            u = 1.0
            timer = MPI_WTIME()
            call ma28ad(dorder,dnz,damat,dlicn,dirown,dlirn,dicoln,
     &           u,dikeep,diw,dwgt,ierr)
            timer = MPI_WTIME() - timer
            write (*,*) myid, timer
          end if
          call writelud
        end if
      end if

c  Create new MPI types for easy communication
      call MPI_TYPE_VECTOR(1,1,nsd,MPI_DOUBLE_PRECISION,rowelemtype,ierr)
      do isd = 1,nsd
        if (dbnodes(isd).gt.0) then
          do ib = 1,dbnodes(isd)
            block(ib) = 1
            disp(ib)  = (dblist(ib,isd)-1)*nsd
          end do
          call MPI_TYPE_INDEXED(dbnodes(isd),block,disp,MPI_DOUBLE_PRECISION,
     &         dbctype(isd),ierr)
          call MPI_TYPE_COMMIT(dbctype(isd),ierr)
        end if
      end do
      
      call free(bctempptr)
      call free(blockptr)
      call free(dispptr)

      return
      end
