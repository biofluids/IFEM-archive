      subroutine hypo

      include "global.h"
      include "malloc.h"

      integer ien(nen,nec),cnn(maxconn,nqdc),rng(neface,nec)
      integer ncnn(nqdc),cnn2(maxconn,nnc),ncnn2(nnc)
      real* 8 xn(nsd,nnc), x(nsd,nn_on), xloc(nsd,nn_loc)
      real* 8 hg(nec), xna(nsd,nn)
      real* 8 shrk(0:nsd,maxconn,nquad*nec)
      real* 8 shrknode(maxconn,nnc)
      pointer (rngptr,rng),(ienptr,ien),(cnnptr,cnn),(ncnnptr,ncnn)
      pointer (cnn2ptr,cnn2),(ncnn2ptr,ncnn2)
      pointer (shrkptr,shrk)
      pointer (shrknodeptr,shrknode)
      pointer (xnptr,xn),(xptr,x),(xlocptr,xloc)
      pointer (hgptr,hg),(xnaptr,xna)

      real* 8 nodebc(ndf,nnc),nodebcon2(ndf,nn_on2)
      real* 8 nodebcv(ndf,nnc)
      integer bmap(nn,ndf)
c      integer blist(ndf,1)
c      real* 8 ac(1,ndf)
c      integer irnc(1,ndf),icnc(1,ndf)
      real* 8, dimension(:,:), allocatable :: ac
      integer, dimension(:,:), allocatable :: blist, bnodesall, irnc, icnc
      real* 8, dimension(:,:), allocatable :: bvlist
      pointer (nodebcptr,nodebc)
      pointer (nodebcon2ptr,nodebcon2)
      pointer (nodebcvptr,nodebcv),(bmapptr,bmap)
c      pointer (blistptr,blist)
c      pointer (acptr,ac),(irncptr,irnc),(icncptr,icnc)

      integer id(ndf,nnc)
      real* 8 dbarn(ndf,nnc),dbar(ndf,nn_on),dbar2(ndf,nn_on2)
      real* 8 dn(ndf,nnc),un(ndf,nnc)
      real* 8 d(ndf,nn_on),do(ndf,nn_on),u(ndf,nn_on)
      real* 8 d2(ndf,nn_on2)
      real* 8 dd(ndf,nnc),dd2(ndf+1,nnc)
      real* 8 bg(ndf,nnc),dg(ndf,nnc),p(ndf,nn_on)
      real* 8 wg(ndf,nnc), w(ndf,nn_on)
      real* 8 p2(ndf,nn_on2),w2(ndf,nn_on2)
      pointer (idptr,id),(fnptr,fn)
      pointer (dbarnptr,dbarn),(dbarptr,dbar),(dbar2ptr,dbar2)
      pointer (dnptr,dn),(unptr,un)
      pointer (dptr,d),(doptr,do),(uptr,u)
      pointer (d2ptr,d2)
      pointer (ddptr,dd),(dd2ptr,dd2)
      pointer (bgptr,bg),(dgptr,dg),(pptr,p)
      pointer (wgptr,wg),(wptr,w)
      pointer (p2ptr,p2),(w2ptr,w2)

      real* 8  hn(nnc),hm(nn_on),hm2(nn_on2),hloc(nn_loc)
      pointer (hnptr,hn),(hmptr,hm),(hm2ptr,hm2),(hlocptr,hloc)

      real* 8 z(ndf*nnc,inner)
      real* 8 v(ndf*nnc,inner+1)
      real* 8 zg(ndf*nnc), avg(ndf*nnc), sm(ndf*nnc)
      real* 8 von(ndf,nn_on),avon(ndf,nn_on)
      real* 8 vn(ndf,nnc),v2(ndf,nn_on2)
      pointer (zptr,z),(vptr,v)
      pointer (avgptr,avg),(zgptr,zg),(smptr,sm)
      pointer (vonptr,von),(avonptr,avon)
      pointer (vnptr,vn),(v2ptr,v2)

      real* 8 h(inner+1,inner)
      real* 8 y(inner+1)
      real* 8 cc(inner), ss(inner)
      pointer (hptr,h),(yptr,y),(ccptr,cc),(ssptr,ss)

      logical assemble,homog
      integer ierr
      real* 8 starttime,endtime,totaltime

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("initialization",-999,.false.)
      tt = t_start
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("read ien",-999,.false.)
      ienptr = malloc(nen*nec*isize)
      rngptr = malloc(nec*neface*isize)
      call readien(ien)
      call readrng(rng)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      cnnptr = malloc(nqdc*maxconn*isize)      
      cnn2ptr = malloc(nnc*maxconn*isize)
      ncnnptr = malloc(nqdc*isize)
      ncnn2ptr = malloc(nnc*isize)
      shrkptr = malloc(nec*nquad*maxconn*(nsd+1)*fsize)
      shrknodeptr = malloc(nnc*maxconn*fsize)
      if(readshape) then
        call error("read tmp data",-999,.false.)
        call readtmp(shrk,shrknode,cnn,ncnn,cnn2,ncnn2)
      else
        xnaptr = malloc(nn*nsd*fsize)
        call readallx(xna)
        call error("rkpm shape function",-999,.false.)
        call rkpm(shrk,shrknode,cnn,ncnn,cnn2,ncnn2,ien,rng,xna)
        call error("write tmp data",-999,.false.)
        call writetmp(shrk,shrknode,cnn,ncnn,cnn2,ncnn2)
        call free(xnaptr)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("communications",-999,.false.)
      call commsetup(ien)
      call nodesetup(cnn,ncnn,maxconn)
      call nodesetup2(cnn2,ncnn2,maxconn)

      xnptr = malloc(nnc*nsd*fsize)
      xptr = malloc(nn_on*nsd*fsize)
      xlocptr = malloc(nn_loc*nsd*fsize)

      hgptr = malloc(nec*fsize)

      nodebcptr = malloc(nnc*ndf*fsize)
      nodebcon2ptr = malloc(nn_on2*ndf*fsize)
      nodebcvptr = malloc(nnc*ndf*fsize)
      bmapptr = malloc(nn*ndf*isize)
      idptr = malloc(nnc*ndf*isize)
      dbarnptr = malloc(nnc*ndf*fsize)
      dnptr = malloc(nnc*ndf*fsize)
      fnptr = malloc(nnc*ndf*fsize)
      unptr = malloc(nnc*ndf*fsize)
      bgptr = malloc(nnc*ndf*fsize)
      p2ptr = malloc(nn_on2*ndf*fsize)
      dgptr = malloc(nnc*ndf*fsize)
      wgptr = malloc(nnc*ndf*fsize)
      w2ptr = malloc(nn_on2*ndf*fsize)
      ddptr = malloc(nnc*ndf*fsize)
      dd2ptr = malloc(nnc*(ndf+1)*fsize)
      dptr = malloc(nn_on*ndf*fsize)
      d2ptr = malloc(nn_on2*ndf*fsize)
      doptr = malloc(nn_on*ndf*fsize)
      pptr = malloc(nn_on*ndf*fsize)
      wptr = malloc(nn_on*ndf*fsize)
      dbarptr = malloc(nn_on*ndf*fsize)
      dbar2ptr = malloc(nn_on2*ndf*fsize)
      uptr = malloc(nn_on*ndf*fsize)

c      idfptr = malloc(nnc*isize)
c      dnfptr = malloc(nnc*fsize)
c      fnfptr = malloc(nnc*fsize)
c      bgfptr = malloc(nnc*fsize)
c      dgfptr = malloc(nnc*fsize)
c      wgfptr = malloc(nnc*fsize)
c      dfptr = malloc(nn_on*fsize)
c      dofptr = malloc(nn_on*fsize)
c      pfptr = malloc(nn_on*fsize)
c      wfptr = malloc(nn_on*fsize)

      hnptr = malloc(nnc*fsize)
      hmptr = malloc(nn_on*fsize)
      hm2ptr = malloc(nn_on2*fsize)
      hlocptr = malloc(nn_loc*fsize)

      hptr = malloc((inner+1)*inner*fsize)
      yptr = malloc((inner+1)*fsize)
      ccptr = malloc(inner*fsize)
      ssptr = malloc(inner*fsize)

      zptr = malloc(ndf*nnc*inner*fsize)
      vptr = malloc(ndf*nnc*(inner+1)*fsize)
      zgptr = malloc(ndf*nnc*fsize)
      avgptr = malloc(ndf*nnc*fsize)
      smptr = malloc(ndf*nnc*fsize)
 	  vonptr = malloc(ndf*nn_on*fsize)
      avonptr = malloc(ndf*nn_on*fsize)
      vnptr = malloc(ndf*nnc*fsize)
      v2ptr = malloc(ndf*nn_on2*fsize)

c      zfptr = malloc(nnc*inner*fsize)
c      vfptr = malloc(nnc*(inner+1)*fsize)
c      zgfptr = malloc(nnc*fsize)
c	   avgfptr = malloc(nnc*fsize)
c      smfptr = malloc(nnc*fsize)
c      vonfptr = malloc(nn_on*fsize)
c      avonfptr = malloc(nn_on*fsize)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("read x ",-999,.false.)
      call readx(xn)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("boundary condition setup",-999,.false.)
      allocate(bnodesall(0:numproc-1,ndf))
      call getnodebc(nodebc,nodebcon2,nodebcv,bnodesall,
     &     ien,rng,hn,hm2,hloc)
      if (maxval(bnodes).gt.0) then
        allocate(blist(maxval(bnodes),ndf),
     &       bvlist(maxval(bnodes),ndf),
     &       ac(maxval(bnodes)*maxconn,ndf),
     &       irnc(maxval(bnodes)*maxconn,ndf),
     &       icnc(maxval(bnodes)*maxconn,ndf))
      end if
      call bcsetup(bmap,blist,bvlist,bnodesall,ac,irnc,icnc,maxval(bnodes),
     &     nodebc,nodebcon2,nodebcv,cnn2,ncnn2,shrknode)
      call free(nodebcvptr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("form d and id",-999,.false.)

      call formd (dbarn)

      if(restart) then
        call error("restart",-999,.false.)
        call diskin(dbarn,dd)
      endif
      
      call error("diskout",-999,.false.)
      
      starttime = MPI_WTIME()
      homog = .false.
      call dbartod(dbarn,dn,d,d2,blist,bvlist,bnodesall,nodebcon2,homog,
     &     maxval(bnodes),shrknode,cnn2,ncnn2,hn,hm,hm2)
      endtime = MPI_WTIME()
      totaltime = endtime - starttime
      write (*,*) myid,"time",totaltime
      
      call getun(un,dn,d2,shrknode,cnn2,ncnn2,hn,hm2)

      call diskout(un,dd)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("local shape functions",-999,.false.)
      call shape
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("localization",-999,.false.)
      call gather (xloc,xn,nsd,hn,hloc)
      call grab_all (x,xn,nsd,hn,hm)
      call grab_all (d,dn,ndf,hn,hm)
      call grab_all2 (d2,dn,ndf,hn,hm2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("mesh information",-999,.false.)
	  call lenght(xloc,ien,hg)
	  if (myid.eq.0) then
        write(7,'(" Minimum element lenght... = ",e15.8)') hmin
        write(7,'(" Maximum element lenght... = ",e15.8)') hmax
        write(7,'(" Minimum element volume... = ",e15.8)') vmin
        write(7,'(" Maximum element volume... = ",e15.8)') vmax
	  endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      nq = 0
      if (myid.eq.0) write (7,101) nq
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do its=1,nts
        if (myid.eq.0) write (6,*) ' '
        if (myid.eq.0) write (6,*) ' TIME STEP = ', its
        if (myid.eq.0) write (7,*) ' '
        if (myid.eq.0) write (7,*) ' TIME STEP = ', its
        tt = tt + dt
        call equal(d,do,ndf*nn_on)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do iit=1,nit
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          call fclear(p,ndf*nn_on)
          call fclear(w,ndf*nn_on)
          call block(xloc,shrk,d,do,p,w,hg,ien,rng,cnn,ncnn)
c          assemble=.true.
          call send_all(p,bg,ndf,assemble,hn,hm)
          call send_all(w,wg,ndf,assemble,hn,hm)
          call rtorbar(p,bg,p2,dd,blist,bvlist,bnodesall,nodebcon2,
     &         maxval(bnodes),shrknode,cnn2,ncnn2,hn,hm,hm2)
          call rtorbar(w,wg,w2,dd,blist,bvlist,bnodesall,nodebcon2,
     &         maxval(bnodes),shrknode,cnn2,ncnn2,hn,hm,hm2)
          call getnorm(bg,bg,ndf*nnc,res_l)
          res_l= sqrt(res_l/nq)
          call fclear(dg,ndf*nnc)
          call gmres(xloc,shrk,shrknode,d,do,wg,bg,p2,dd,dg,hg,ien,rng,
     &         blist,bvlist,maxval(bnodes),
     &         bnodesall,nodebcon2,cnn,ncnn,cnn2,ncnn2,hn,hm,hm2,
     &         z,v,zg,vn,v2,avg,sm,von,avon,h,y,cc,ss)
          call getnorm(dg,dg,ndf*nnc,del_l)
          del_l = sqrt(del_l/nq)
          call update(dbarn,dg)
          call dbartod(dbarn,dn,d,d2,blist,bvlist,bnodesall,nodebcon2,.false.,
     &         maxval(bnodes),shrknode,cnn2,ncnn2,hn,hm,hm2)
          if (myid.eq.0) write(6,102)iit,res_l,del_l
          if (myid.eq.0) write(7,102)iit,res_l,del_l
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(mod(its,ntsbout).eq.0) then
          idisk = idisk+1
          assemble=.false.
          call send_all(d,dn,ndf,assemble,hn,hm)
          call getun(un,dn,d2,shrknode,cnn2,ncnn2,hn,hm2)
          call diskout(un,dd)
          call restartout(dbarn,dd)
        endif
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (myid.eq.0) then
        write(7,*)'Processors...............',numproc 
      endif


      return

 101  format(/"Number of equations for Flow.........(nq) = ",i10)
 102  format("Iteration",i3,':  ',2e14.7)

      end



