      subroutine hypo

      include "global.h"
      include "malloc.h"

      integer ien(nen,nec),rng(neface,nec)
      real* 8 xn(nsd,nnc), x(nsd,nn_loc), hg(nec)
      pointer (rngptr,rng),(ienptr,ien)
      pointer (xnptr,xn),(xptr,x),(hgptr,hg)

      integer id(ndf,nnc) 
      real* 8 dn(ndf,nnc),fn(ndf,nnc)
      real* 8 d(ndf,nn_loc),do(ndf,nn_loc)
      real* 8 bg(ndf,nnc),dg(ndf,nnc),p(ndf,nn_loc)
      real* 8 wg(ndf,nnc), w(ndf,nn_loc)
      pointer (idptr,id),(fnptr,fn),(dnptr,dn)
      pointer (dptr,d),(doptr,do),(uptr,u)
      pointer (bgptr,bg),(dgptr,dg),(pptr,p)
      pointer (wgptr,wg),(wptr,w)

      real* 8 dd(ndf,nnc),hn(nnc),hm(nn_loc)
      pointer (ddptr,dd),(hnptr,hn),(hmptr,hm)

      real* 8 z(ndf*nnc,inner)
      real* 8 v(ndf*nnc,inner+1)
      real* 8 zg(ndf*nnc), avg(ndf*nnc), sm(ndf*nnc)
      real* 8 vloc(ndf,nn_loc),avloc(ndf,nn_loc)
      pointer (zptr,z),(vptr,v)
      pointer (avgptr,avg),(zgptr,zg),(smptr,sm)
      pointer (vlocptr,vloc),(avlocptr,avloc)

      real* 8 h(inner+1,inner)
      real* 8 y(inner+1)
      real* 8 cc(inner), ss(inner)
      pointer (hptr,h),(yptr,y),(ccptr,cc),(ssptr,ss)

      logical assemble
      integer ierr

      real* 8 fdrag(nsdpad)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("initialization",-999,.false.)
      tt = t_start
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("read ien",-999,.false.)
      ienptr = malloc(nen*nec*isize)
      call readien(ien)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("communications",-999,.false.)
      call commsetup(ien)
c  call nodesetup (cnn,mn) !!!!change me

      xnptr = malloc(nnc*nsd*fsize)
      xptr = malloc(nn_loc*nsd*fsize)
      rngptr = malloc(nec*neface*isize)
      hgptr = malloc(nec*fsize)

      idptr = malloc(nnc*ndf*isize)
      dnptr = malloc(nnc*ndf*fsize)
      fnptr = malloc(nnc*ndf*fsize)
      bgptr = malloc(nnc*ndf*fsize)
      dgptr = malloc(nnc*ndf*fsize)
      wgptr = malloc(nnc*ndf*fsize)
      dptr = malloc(nn_loc*ndf*fsize)
      uptr = malloc(nn_loc*nsd*fsize)
      doptr = malloc(nn_loc*ndf*fsize)
      pptr = malloc(nn_loc*ndf*fsize)
      wptr = malloc(nn_loc*ndf*fsize)

      ddptr = malloc(ndf*nnc*fsize)
      hnptr = malloc(nnc*fsize)
      hmptr = malloc(nn_loc*fsize)

      hptr = malloc((inner+1)*inner*fsize)
      yptr = malloc((inner+1)*fsize)
      ccptr = malloc(inner*fsize)
      ssptr = malloc(inner*fsize)

      zptr = malloc(ndf*nnc*inner*fsize)
      vptr = malloc(ndf*nnc*(inner+1)*fsize)
      zgptr = malloc(ndf*nnc*fsize)
      avgptr = malloc(ndf*nnc*fsize)
      smptr = malloc(ndf*nnc*fsize)
 	  vlocptr = malloc(ndf*nn_loc*fsize)
      avlocptr = malloc(ndf*nn_loc*fsize)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("read x and rng",-999,.false.)
      call readx(xn)
      call readrng(rng)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("form d and id",-999,.false.)
      call formid(id,rng,ien,hn,hm)
      call formd (xn,fn,rng,ien,hn,hm)
      call equal(fn ,dn ,ndf*nnc)
      if(restart) then
        call error("restart",-999,.false.)
        call diskin(dn,dd)
	  endif
      call setd(dn ,fn ,id ,ndf)
      call error("diskout",-999,.false.)

      call diskout(dn,dd)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("local shape functions",-999,.false.)
      call shape
      call shape2d
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("localization",-999,.false.)
      call gather (x,xn,nsd,hn,hm)
      call gather (d,dn,ndf,hn,hm)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (calcforce) then
        call error("lift and drag force",-999,.false.)
        call getfdrag(fdrag,d,x,ien,rng,hn,hm)
        call fdragout(tt,fdrag,1,3)
      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call error("mesh information",-999,.false.)
      call vol(x,ien)
	  if (myid.eq.0) then
        write(7,*)' '
        write(7,'(" Volume of liquid......... = ",e15.8)') liq
	  endif
	  call lenght(x,ien,hg)
	  if (myid.eq.0) then
        write(7,'(" Minimum element lenght... = ",e15.8)') hmin
        write(7,'(" Maximum element lenght... = ",e15.8)') hmax
        write(7,'(" Minimum element volume... = ",e15.8)') vmin
        write(7,'(" Maximum element volume... = ",e15.8)') vmax
	  endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
        call equal(d,do,ndf*nn_loc)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do iit=1,nit
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          call fclear(p,ndf*nn_loc)
          call fclear(w,ndf*nn_loc)
          call block(x,d,do,p,w,hg,ien)
          assemble=.true.
          call scatter(p,bg,ndf,assemble,hn,hm)
          call scatter(w,wg,ndf,assemble,hn,hm)
          call setid(bg,id,ndf)
          call getnorm(bg,bg,ndf*nnc,res_l)
          res_l= sqrt(res_l/nq)
          call fclear(dg,ndf*nnc)
          call gmres(x,d,do,id,wg,bg,dg,hg,ien,hn,hm,
     &         z,v,zg,avg,sm,vloc,avloc,h,y,cc,ss)
          call getnorm(dg,dg,ndf*nnc,del_l)
          del_l = sqrt(del_l/nq)
          call update(p, d, dg, ndf,hn,hm)
          if (myid.eq.0) write(6,102)iit,res_l,del_l
          if (myid.eq.0) write(7,102)iit,res_l,del_l
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (calcforce) then
          call error("lift and drag force",-999,.false.)
          call getfdrag(fdrag,d,x,ien,rng,hn,hm)
          call fdragout(tt,fdrag,1,3)
        end if
        if(mod(its,ntsbout).eq.0) then
          if(surf(0).gt.0) call dyno(x,d,rng,ien)
          idisk = idisk+1
          assemble=.false.
          call scatter(d,dn,ndf,assemble,hn,hm)
          call diskout(dn,dd)
        endif
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (myid.eq.0) then
        write(7,*)'Processors...............',numproc 
      endif

      return

 101  format(/"Number of equations for Flow.........(nq) = ",i10)
 102  format("Iteration",i3,':  ',4e14.7)

      end
